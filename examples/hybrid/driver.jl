import ClimaAtmos as CA

# load in default CLI arguments
(s, parsed_args) = CA.parse_commandline()

# overwrite default CLI arguments
parsed_args["dt"] = "580secs"
parsed_args["moist"] = "equil"
parsed_args["precip_model"] = "0M"
parsed_args["apply_limiter"] = false
parsed_args["surface_scheme"] = "bulk"
parsed_args["vert_diff"] = true
parsed_args["rad"] = "gray"
parsed_args["C_E"] = Float64(0.0044)
parsed_args["kappa_4"] = Float64(2e17)
parsed_args["rayleigh_sponge"] = true # rayleigh spponge at the top
parsed_args["alpha_rayleigh_uh"] = 0
parsed_args["zd_rayleigh"] = 30e3
const FT = parsed_args["FLOAT_TYPE"] == "Float64" ? Float64 : Float32

# give your job a name
job_id = parsed_args["job_id"] = "my_first_run"


# set the length of the simulation, and the appropriate output parameters
debug_mode = false
if debug_mode
    parsed_args["dt_save_to_sol"] = "1days" # frequency at which the output is saved to cache
    parsed_args["dt_save_to_disk"] = "1days" # frequency at which the output is saved to hdf5 files
    parsed_args["t_end"] = "2days" # end of simulation
    diagnostics_day_range = collect(0:1:2) # range of days for averaging post-processed disgnostics output
else
    parsed_args["dt_save_to_sol"] = "2days" # frequency at which the output is saved to cache
    parsed_args["dt_save_to_disk"] = "2days" # frequency at which the output is saved to hdf5 files
    parsed_args["t_end"] = "100days" # end of simulation
    diagnostics_day_range = collect(0:50:100) # range of days for averaging post-processed disgnostics output
end

# load parameters
include("parameter_set.jl")
params, parsed_args = create_parameter_set(FT, parsed_args, CA.cli_defaults(s))

# setup communication for parallel computing
include("comms.jl")
fps = parsed_args["fps"]

# import local modules
import ClimaAtmos.RRTMGPInterface as RRTMGPI
import ClimaAtmos.InitialConditions as ICs
import ClimaAtmos.TurbulenceConvection as TC

# instantiate the needded packages
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

# include any external artifacts
include(joinpath(pkgdir(CA), "artifacts", "artifact_funcs.jl"))

# define the model mathematical formulation (with boundary conditions), numerics, timestepping and initial condition specification
atmos = CA.get_atmos(FT, parsed_args, params.turbconv_params)
@info "AtmosModel: \n$(summary(atmos))"
numerics = CA.get_numerics(parsed_args)
simulation = CA.get_simulation(FT, parsed_args)
initial_condition = CA.get_initial_condition(parsed_args)

# initialize the model's prognostic state vector
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

# initialize the model's cached auxiliary information and variables
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

# initialize the timestepping method
@time "ode_configuration" ode_algo = CA.ode_configuration(Y, parsed_args, atmos)

@time "get_callbacks" callback =
    CA.get_callbacks(parsed_args, simulation, atmos, params)
tspan = (t_start, simulation.t_end)

# initialize the entire simulation
@time "args_integrator" integrator_args, integrator_kwargs =
    CA.args_integrator(parsed_args, Y, p, tspan, ode_algo, callback)
@time "get_integrator" integrator =
    CA.get_integrator(integrator_args, integrator_kwargs)

# logger struct
struct SimulationResults{S, RT, WT}
    sol::S
    ret_code::RT
    walltime::WT
end

# print to screen
@info "Running" job_id = simulation.job_id output_dir = simulation.output_dir tspan

# main stepping function
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

# run the model
sol_res = perform_solve!(integrator, simulation, comms_ctx);

# # inspect your data
cwd = @__DIR__
include("$cwd/../../post_processing/plot/plot_model_data.jl")
# end step plots
png_end(
    integrator.sol,
    :meridional_wind,
    level = 5,
    output_dir = "$cwd/../output/",
)
png_end(integrator.sol, :zonal_wind, level = 5, output_dir = "$cwd/../output/") # end step
# animation
mp4(
    integrator.sol,
    :total_energy,
    level = 2,
    clims = (0, 1e5),
    output_dir = "$cwd/../output/",
)

# more serious post-processing
more_serious_post_processing = true
time_mean_post_processing = true
if more_serious_post_processing
    cmd = `julia --color=yes --project $cwd/../../post_processing/remap/remap_pipeline.jl --data_dir $cwd/../output/$job_id --out_dir $cwd/../output/nc_output_$job_id`
    run(cmd)
    cmd = `julia --color=yes --project $cwd/../../post_processing/plot/plot_pipeline.jl --nc_dir $cwd/../output/nc_output_$job_id --fig_dir $cwd/../output/nc_plots_$job_id --case_name moist_baroclinic_wave`
    run(cmd)
    # time-mean plots
    if time_mean_post_processing
        include("$cwd/../../post_processing/plot/time_mean_plots.jl")
        plot_time_mean(
            :temperature,
            nc_dir = "$cwd/../output/nc_output_$job_id",
            fig_dir = "$cwd/../output/nc_plots_$job_id",
            day_range = diagnostics_day_range,
        )
        plot_time_mean(
            :zonal_wind,
            nc_dir = "$cwd/../output/nc_output_$job_id",
            fig_dir = "$cwd/../output/nc_plots_$job_id",
            day_range = diagnostics_day_range,
        )
    end
end
