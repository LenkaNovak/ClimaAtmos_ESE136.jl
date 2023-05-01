import ClimaAtmos as CA

(s, parsed_args_defaults) = CA.parse_commandline()
if !(@isdefined parsed_args)
    parsed_args = parsed_args_defaults
end

# TODO: can we move this into src/?
using NVTX
using ClimaComms
if NVTX.isactive()
    NVTX.enable_gc_hooks()
    # makes output on buildkite a bit nicer
    if ClimaComms.iamroot(comms_ctx)
        atexit() do
            println("--- Saving profiler information")
        end
    end
end

parse_arg(pa, key, default) = isnothing(pa[key]) ? default : pa[key]

const FT = parsed_args["FLOAT_TYPE"] == "Float64" ? Float64 : Float32

include("parameter_set.jl")
params, parsed_args = create_parameter_set(FT, parsed_args, CA.cli_defaults(s))

include("comms.jl")
fps = parsed_args["fps"]
import ClimaAtmos.RRTMGPInterface as RRTMGPI
import ClimaAtmos.InitialConditions as ICs

include(joinpath(pkgdir(CA), "artifacts", "artifact_funcs.jl"))

import ClimaAtmos.TurbulenceConvection as TC

atmos = CA.get_atmos(FT, parsed_args, params.turbconv_params)
@info "AtmosModel: \n$(summary(atmos))"
numerics = CA.get_numerics(parsed_args)
simulation = CA.get_simulation(FT, parsed_args)
initial_condition = CA.get_initial_condition(parsed_args)

# TODO: use import instead of using
using Colors
using OrdinaryDiffEq
using PrettyTables
using DiffEqCallbacks
using JLD2
using ClimaCore.DataLayouts
using NCDatasets
using ClimaCore
using ClimaTimeSteppers

import Random
Random.seed!(1234)

using OrdinaryDiffEq
using DiffEqCallbacks
using JLD2

import ClimaCore
if parsed_args["trunc_stack_traces"]
    ClimaCore.Fields.truncate_printing_field_types() = true
end

using Statistics: mean
import SurfaceFluxes as SF
using CloudMicrophysics
const CCG = ClimaCore.Geometry
import ClimaAtmos.TurbulenceConvection as TC
import ClimaCore.Operators as CCO
const CM = CloudMicrophysics
import ClimaAtmos.Parameters as CAP

import ClimaCore: enable_threading
const enable_clima_core_threading = parsed_args["enable_threading"]
enable_threading() = enable_clima_core_threading

@time "Allocating Y" if simulation.restart
    (Y, t_start) = CA.get_state_restart(comms_ctx)
    spaces = CA.get_spaces_restart(Y)
else
    spaces = CA.get_spaces(parsed_args, params, comms_ctx)
    Y = ICs.atmos_state(
        initial_condition(params),
        atmos,
        spaces.center_space,
        spaces.face_space,
    )
    t_start = FT(0)
end

@time "Allocating cache (p)" begin
    p = CA.get_cache(
        Y,
        parsed_args,
        params,
        spaces,
        atmos,
        numerics,
        simulation,
        initial_condition,
        comms_ctx,
    )
end

if parsed_args["discrete_hydrostatic_balance"]
    CA.set_discrete_hydrostatic_balanced_state!(Y, p)
end

@time "ode_configuration" ode_algo = CA.ode_configuration(Y, parsed_args, atmos)

@time "get_callbacks" callback =
    CA.get_callbacks(parsed_args, simulation, atmos, params)
tspan = (t_start, simulation.t_end)
@time "args_integrator" integrator_args, integrator_kwargs =
    CA.args_integrator(parsed_args, Y, p, tspan, ode_algo, callback)

if haskey(ENV, "CI_PERF_SKIP_INIT") # for performance analysis
    throw(:exit_profile_init)
end

@time "get_integrator" integrator =
    CA.get_integrator(integrator_args, integrator_kwargs)

if haskey(ENV, "CI_PERF_SKIP_RUN") # for performance analysis
    throw(:exit_profile)
end

@info "Running" job_id = simulation.job_id output_dir = simulation.output_dir tspan

struct SimulationResults{S, RT, WT}
    sol::S
    ret_code::RT
    walltime::WT
end
function perform_solve!(integrator, simulation, comms_ctx)
    try
        if simulation.is_distributed
            OrdinaryDiffEq.step!(integrator)
            # GC.enable(false) # disabling GC causes a memory leak
            GC.gc()
            ClimaComms.barrier(comms_ctx)
            if ClimaComms.iamroot(comms_ctx)
                @timev begin
                    walltime = @elapsed sol = OrdinaryDiffEq.solve!(integrator)
                end
            else
                walltime = @elapsed sol = OrdinaryDiffEq.solve!(integrator)
            end
            ClimaComms.barrier(comms_ctx)
            GC.enable(true)
            return SimulationResults(sol, :success, walltime)
        else
            sol = @timev OrdinaryDiffEq.solve!(integrator)
            return SimulationResults(sol, :success, nothing)
        end
    catch ret_code
        @error "ClimaAtmos simulation crashed. Stacktrace for failed simulation" exception =
            (ret_code, catch_backtrace())
        return SimulationResults(nothing, :simulation_crashed, nothing)
    end
end

sol_res = perform_solve!(integrator, simulation, comms_ctx)

