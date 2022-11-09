import ClimaCore.DataLayouts as DL
import ClimaAtmos.RRTMGPInterface as RRTMGPI
import Thermodynamics as TD
import LinearAlgebra
import ClimaCore.Fields
import ClimaComms
import ClimaCore as CC
import ClimaCore.Operators as CCO
import ClimaCore.Spaces
import OrdinaryDiffEq as ODE
import ClimaAtmos.Parameters as CAP
import DiffEqCallbacks as DEQ
import ClimaCore: InputOutput
import Dates
using Insolation: instantaneous_zenith_angle

function get_callbacks(parsed_args, simulation, atmos, params)
    FT = eltype(params)
    (; dt) = simulation

    callback_filters = call_every_n_steps(affect_filter!; skip_first = true)
    tc_callbacks =
        call_every_n_steps(turb_conv_affect_filter!; skip_first = true)

    additional_callbacks =
        if atmos.radiation_mode isa RRTMGPI.AbstractRRTMGPMode
            # TODO: better if-else criteria?
            dt_rad = if parsed_args["config"] == "column"
                dt
            else
                FT(time_to_seconds(parsed_args["dt_rad"]))
            end
            (call_every_dt(rrtmgp_model_callback!, dt_rad),)
        else
            ()
        end

    if !isnothing(atmos.turbconv_model)
        additional_callbacks = (additional_callbacks..., tc_callbacks)
    end
    if atmos.moisture_model isa EquilMoistModel &&
       parsed_args["apply_moisture_filter"]
        additional_callbacks = (additional_callbacks..., callback_filters)
    end

    dt_save_to_disk = time_to_seconds(parsed_args["dt_save_to_disk"])
    dt_save_restart = time_to_seconds(parsed_args["dt_save_restart"])

    dss_cb = if startswith(parsed_args["ode_algo"], "ODE.")
        call_every_n_steps(dss_callback)
    else
        nothing
    end
    save_to_disk_callback = if dt_save_to_disk == Inf
        nothing
    else
        call_every_dt(save_to_disk_func, dt_save_to_disk)
    end

    save_restart_callback = if dt_save_restart == Inf
        nothing
    else
        call_every_dt(save_restart_func, dt_save_restart)
    end

    gc_callback = if simulation.is_distributed
        call_every_n_steps(
            gc_func,
            parse(Int, get(ENV, "CLIMAATMOS_GC_NSTEPS", "1000")),
            skip_first = true,
        )
    else
        nothing
    end

    return ODE.CallbackSet(
        dss_cb,
        save_to_disk_callback,
        save_restart_callback,
        gc_callback,
        additional_callbacks...,
    )
end

function call_every_n_steps(f!, n = 1; skip_first = false, call_at_end = false)
    previous_step = Ref(0)
    return ODE.DiscreteCallback(
        (u, t, integrator) ->
            (previous_step[] += 1) % n == 0 ||
                (call_at_end && t == integrator.sol.prob.tspan[2]),
        f!;
        initialize = (cb, u, t, integrator) -> skip_first || f!(integrator),
        save_positions = (false, false),
    )
end

function call_every_dt(f!, dt; skip_first = false, call_at_end = false)
    next_t = Ref{typeof(dt)}()
    affect! = function (integrator)
        f!(integrator)

        t = integrator.t
        t_end = integrator.sol.prob.tspan[2]
        next_t[] = max(t, next_t[] + dt)
        if call_at_end
            next_t[] = min(next_t[], t_end)
        end
    end
    return ODE.DiscreteCallback(
        (u, t, integrator) -> t >= next_t[],
        affect!;
        initialize = (cb, u, t, integrator) -> begin
            skip_first || f!(integrator)
            t_end = integrator.sol.prob.tspan[2]
            next_t[] =
                (call_at_end && t < t_end) ? min(t_end, t + dt) : t + dt
        end,
        save_positions = (false, false),
    )
end

function affect_filter!(Y::Fields.FieldVector)
    @. Y.c.ρq_tot = max(Y.c.ρq_tot, 0)
    return nothing
end

function dss_callback(integrator)
    Y = integrator.u
    ghost_buffer = integrator.p.ghost_buffer
    @nvtx "dss callback" color = colorant"yellow" begin
        Spaces.weighted_dss_start!(Y.c, ghost_buffer.c)
        Spaces.weighted_dss_start!(Y.f, ghost_buffer.f)
        Spaces.weighted_dss_internal!(Y.c, ghost_buffer.c)
        Spaces.weighted_dss_internal!(Y.f, ghost_buffer.f)
        Spaces.weighted_dss_ghost!(Y.c, ghost_buffer.c)
        Spaces.weighted_dss_ghost!(Y.f, ghost_buffer.f)
    end
    # ODE.u_modified!(integrator, false) # TODO: try this
end


function affect_filter!(integrator)
    (; apply_moisture_filter) = integrator.p
    affect_filter!(integrator.u)
    # We're lying to OrdinaryDiffEq.jl, in order to avoid
    # paying for an additional tendency call, which is required
    # to support supplying a continuous representation of the solution.
    ODE.u_modified!(integrator, false)
end

function turb_conv_affect_filter!(integrator)
    (; edmf_cache, Δt) = p
    (; edmf, param_set, aux, case, surf_params) = edmf_cache
    t = integrator.t
    Y = integrator.u
    tc_params = CAP.turbconv_params(param_set)

    Fields.bycolumn(axes(Y.c)) do colidx
        state = TC.tc_column_state(Y, p, nothing, colidx)
        grid = TC.Grid(state)
        surf = TCU.get_surface(surf_params, grid, state, t, tc_params)
        TC.affect_filter!(edmf, grid, state, tc_params, surf, t)
    end

    # We're lying to OrdinaryDiffEq.jl, in order to avoid
    # paying for an additional `∑tendencies!` call, which is required
    # to support supplying a continuous representation of the
    # solution.
    ODE.u_modified!(integrator, false)
end

function rrtmgp_model_callback!(integrator)
    Y = integrator.u
    p = integrator.p
    t = integrator.t

    (; ᶜK, ᶜts, T_sfc, params, thermo_dispatcher) = p
    (; idealized_insolation, idealized_h2o, idealized_clouds) = p
    (; insolation_tuple, ᶠradiation_flux, radiation_model) = p
    (; ᶜinterp) = p.operators
    thermo_params = CAP.thermodynamics_params(params)
    insolation_params = CAP.insolation_params(params)

    radiation_model.surface_temperature .= RRTMGPI.field2array(T_sfc)

    C123 = Geometry.Covariant123Vector
    ᶜp = RRTMGPI.array2field(radiation_model.center_pressure, axes(Y.c))
    ᶜT = RRTMGPI.array2field(radiation_model.center_temperature, axes(Y.c))
    @. ᶜK = LinearAlgebra.norm_sqr(C123(Y.c.uₕ) + C123(ᶜinterp(Y.f.w))) / 2
    CA.thermo_state!(Y, p, ᶜinterp)
    @. ᶜp = TD.air_pressure(thermo_params, ᶜts)
    @. ᶜT = TD.air_temperature(thermo_params, ᶜts)

    if !(radiation_model.radiation_mode isa RRTMGPI.GrayRadiation)
        ᶜvmr_h2o = RRTMGPI.array2field(
            radiation_model.center_volume_mixing_ratio_h2o,
            axes(Y.c),
        )
        if idealized_h2o
            # slowly increase the relative humidity from 0 to 0.6 to account for
            # the fact that we have a very unrealistic initial condition
            max_relative_humidity = FT(0.6)
            t_increasing_humidity = FT(60 * 60 * 24 * 30)
            if t < t_increasing_humidity
                max_relative_humidity *= t / t_increasing_humidity
            end

            # temporarily store ᶜq_tot in ᶜvmr_h2o
            ᶜq_tot = ᶜvmr_h2o
            @. ᶜq_tot =
                max_relative_humidity * TD.q_vap_saturation(thermo_params, ᶜts)

            # filter ᶜq_tot so that it is monotonically decreasing with z
            for i in 2:Spaces.nlevels(axes(ᶜq_tot))
                level = Fields.field_values(Spaces.level(ᶜq_tot, i))
                prev_level = Fields.field_values(Spaces.level(ᶜq_tot, i - 1))
                @. level = min(level, prev_level)
            end

            # assume that ᶜq_vap = ᶜq_tot when computing ᶜvmr_h2o
            @. ᶜvmr_h2o = TD.shum_to_mixing_ratio(ᶜq_tot, ᶜq_tot)
        else
            @. ᶜvmr_h2o = TD.vol_vapor_mixing_ratio(
                thermo_params,
                TD.PhasePartition(thermo_params, ᶜts),
            )
        end
    end

    if !idealized_insolation
        current_datetime = p.simulation.start_date + Dates.Second(round(Int, t)) # current time
        max_zenith_angle = FT(π) / 2 - eps(FT)
        irradiance = FT(CAP.tot_solar_irrad(params))
        au = FT(CAP.astro_unit(params))

        bottom_coords = Fields.coordinate_field(Spaces.level(Y.c, 1))
        if eltype(bottom_coords) <: Geometry.LatLongZPoint
            solar_zenith_angle = RRTMGPI.array2field(
                radiation_model.solar_zenith_angle,
                axes(bottom_coords),
            )
            weighted_irradiance = RRTMGPI.array2field(
                radiation_model.weighted_irradiance,
                axes(bottom_coords),
            )
            ref_insolation_params = Ref(insolation_params)
            @. insolation_tuple = instantaneous_zenith_angle(
                current_datetime,
                Float64(bottom_coords.long),
                Float64(bottom_coords.lat),
                ref_insolation_params,
            ) # the tuple is (zenith angle, azimuthal angle, earth-sun distance)
            @. solar_zenith_angle =
                min(first(insolation_tuple), max_zenith_angle)
            @. weighted_irradiance =
                irradiance * (au / last(insolation_tuple))^2
        else
            # assume that the latitude and longitude are both 0 for flat space
            insolation_tuple = instantaneous_zenith_angle(
                current_datetime,
                0.0,
                0.0,
                insolation_params,
            )
            radiation_model.solar_zenith_angle .=
                min(first(insolation_tuple), max_zenith_angle)
            radiation_model.weighted_irradiance .=
                irradiance * (au / last(insolation_tuple))^2
        end
    end

    if !idealized_clouds && !(
        radiation_model.radiation_mode isa RRTMGPI.GrayRadiation ||
        radiation_model.radiation_mode isa RRTMGPI.ClearSkyRadiation
    )
        ᶜΔz = Fields.local_geometry_field(Y.c).∂x∂ξ.components.data.:9
        ᶜlwp = RRTMGPI.array2field(
            radiation_model.center_cloud_liquid_water_path,
            axes(Y.c),
        )
        ᶜiwp = RRTMGPI.array2field(
            radiation_model.center_cloud_ice_water_path,
            axes(Y.c),
        )
        ᶜfrac = RRTMGPI.array2field(
            radiation_model.center_cloud_fraction,
            axes(Y.c),
        )
        # multiply by 1000 to convert from kg/m^2 to g/m^2
        @. ᶜlwp =
            1000 * Y.c.ρ * TD.liquid_specific_humidity(thermo_params, ᶜts) * ᶜΔz
        @. ᶜiwp =
            1000 * Y.c.ρ * TD.ice_specific_humidity(thermo_params, ᶜts) * ᶜΔz
        @. ᶜfrac =
            ifelse(TD.has_condensate(thermo_params, ᶜts), FT(1), FT(0) * ᶜΔz)
    end

    RRTMGPI.update_fluxes!(radiation_model)
    RRTMGPI.field2array(ᶠradiation_flux) .= radiation_model.face_flux
end

function save_to_disk_func(integrator)

    (; t, u, p) = integrator
    (; output_dir) = p.simulation
    Y = u

    if :ᶜS_ρq_tot in propertynames(p)
        (;
            ᶜts,
            ᶜp,
            ᶜS_ρq_tot,
            ᶜ3d_rain,
            ᶜ3d_snow,
            params,
            ᶜK,
            col_integrated_rain,
            col_integrated_snow,
            T_sfc,
            q_sfc,
        ) = p
    else
        (; ᶜts, ᶜp, params, ᶜK, T_sfc, q_sfc) = p
    end

    thermo_params = CAP.thermodynamics_params(params)
    cm_params = CAP.microphysics_params(params)

    ᶜuₕ = Y.c.uₕ
    ᶠw = Y.f.w
    # kinetic
    @. ᶜK = norm_sqr(C123(ᶜuₕ) + C123(ᶜinterp(ᶠw))) / 2

    # thermo state
    CA.thermo_state!(Y, p, ᶜinterp)
    @. ᶜp = TD.air_pressure(thermo_params, ᶜts)
    ᶜT = @. TD.air_temperature(thermo_params, ᶜts)
    ᶜθ = @. TD.dry_pottemp(thermo_params, ᶜts)

    # vorticity
    curl_uh = @. curlₕ(Y.c.uₕ)
    ᶜvort = Geometry.WVector.(curl_uh)
    Spaces.weighted_dss!(ᶜvort)

    dry_diagnostic = (;
        pressure = ᶜp,
        temperature = ᶜT,
        potential_temperature = ᶜθ,
        kinetic_energy = ᶜK,
        vorticity = ᶜvort,
        sfc_temperature = T_sfc,
        sfc_qt = q_sfc,
    )

    # cloudwater (liquid and ice), watervapor and RH for moist simulation
    if :ρq_tot in propertynames(Y.c)

        ᶜq = @. TD.PhasePartition(thermo_params, ᶜts)
        ᶜcloud_liquid = @. ᶜq.liq
        ᶜcloud_ice = @. ᶜq.ice
        ᶜwatervapor = @. TD.vapor_specific_humidity(ᶜq)
        ᶜRH = @. TD.relative_humidity(thermo_params, ᶜts)

        moist_diagnostic = (;
            cloud_liquid = ᶜcloud_liquid,
            cloud_ice = ᶜcloud_ice,
            water_vapor = ᶜwatervapor,
            relative_humidity = ᶜRH,
        )
        # precipitation
        if :ᶜS_ρq_tot in propertynames(p)

            @. ᶜS_ρq_tot =
                Y.c.ρ * CM.Microphysics0M.remove_precipitation(
                    cm_params,
                    TD.PhasePartition(thermo_params, ᶜts),
                )

            # rain vs snow
            @. ᶜ3d_rain = ifelse(ᶜT >= FT(273.15), ᶜS_ρq_tot, FT(0))
            @. ᶜ3d_snow = ifelse(ᶜT < FT(273.15), ᶜS_ρq_tot, FT(0))

            CCO.column_integral_definite!(col_integrated_rain, ᶜ3d_rain)
            CCO.column_integral_definite!(col_integrated_snow, ᶜ3d_snow)

            @. col_integrated_rain /= CAP.ρ_cloud_liq(params)
            @. col_integrated_snow /= CAP.ρ_cloud_liq(params)


            moist_diagnostic = (
                moist_diagnostic...,
                precipitation_removal = ᶜS_ρq_tot,
                column_integrated_rain = col_integrated_rain,
                column_integrated_snow = col_integrated_snow,
            )
        end
    else
        moist_diagnostic = NamedTuple()
    end

    if :edmf_cache in propertynames(p) && p.simulation.is_debugging_tc

        tc_cent(p) = p.edmf_cache.aux.cent.turbconv
        tc_face(p) = p.edmf_cache.aux.face.turbconv
        turbulence_convection_diagnostic = (;
            bulk_up_area = tc_cent(p).bulk.area,
            bulk_up_h_tot = tc_cent(p).bulk.h_tot,
            bulk_up_buoyancy = tc_cent(p).bulk.buoy,
            bulk_up_q_tot = tc_cent(p).bulk.q_tot,
            bulk_up_q_liq = tc_cent(p).bulk.q_liq,
            bulk_up_q_ice = tc_cent(p).bulk.q_ice,
            bulk_up_temperature = tc_cent(p).bulk.T,
            bulk_up_cloud_fraction = tc_cent(p).bulk.cloud_fraction,
            bulk_up_e_tot_tendency_precip_formation = tc_cent(
                p,
            ).bulk.e_tot_tendency_precip_formation,
            bulk_up_qt_tendency_precip_formation = tc_cent(
                p,
            ).bulk.qt_tendency_precip_formation,
            env_area = tc_cent(p).en.area,
            env_q_tot = tc_cent(p).en.q_tot,
            env_q_liq = tc_cent(p).en.q_liq,
            env_q_ice = tc_cent(p).en.q_ice,
            env_theta_liq_ice = tc_cent(p).en.θ_liq_ice,
            env_theta_virt = tc_cent(p).en.θ_virt,
            env_theta_dry = tc_cent(p).en.θ_dry,
            env_e_tot = tc_cent(p).en.e_tot,
            env_e_kin = tc_cent(p).en.e_kin,
            env_h_tot = tc_cent(p).en.h_tot,
            env_RH = tc_cent(p).en.RH,
            env_temperature = tc_cent(p).en.T,
            env_buoyancy = tc_cent(p).en.buoy,
            env_cloud_fraction = tc_cent(p).en.cloud_fraction,
            env_TKE = tc_cent(p).en.tke,
            env_Hvar = tc_cent(p).en.Hvar,
            env_QTvar = tc_cent(p).en.QTvar,
            env_HQTcov = tc_cent(p).en.HQTcov,
            env_e_tot_tendency_precip_formation = tc_cent(
                p,
            ).en.e_tot_tendency_precip_formation,
            env_qt_tendency_precip_formation = tc_cent(
                p,
            ).en.qt_tendency_precip_formation,
            env_Hvar_rain_dt = tc_cent(p).en.Hvar_rain_dt,
            env_QTvar_rain_dt = tc_cent(p).en.QTvar_rain_dt,
            env_HQTcov_rain_dt = tc_cent(p).en.HQTcov_rain_dt,
            face_bulk_w = tc_face(p).bulk.w,
            face_env_w = tc_face(p).en.w,
        )
    else
        turbulence_convection_diagnostic = NamedTuple()
    end

    if vert_diff
        (; dif_flux_uₕ, dif_flux_energy, dif_flux_ρq_tot) = p
        vert_diff_diagnostic = (;
            sfc_flux_momentum = dif_flux_uₕ,
            sfc_flux_energy = dif_flux_energy,
            sfc_evaporation = dif_flux_ρq_tot,
        )
    else
        vert_diff_diagnostic = NamedTuple()
    end

    if atmos.radiation_mode isa RRTMGPI.AbstractRRTMGPMode
        (; face_lw_flux_dn, face_lw_flux_up, face_sw_flux_dn, face_sw_flux_up) =
            p.radiation_model
        rad_diagnostic = (;
            lw_flux_down = RRTMGPI.array2field(FT.(face_lw_flux_dn), axes(Y.f)),
            lw_flux_up = RRTMGPI.array2field(FT.(face_lw_flux_up), axes(Y.f)),
            sw_flux_down = RRTMGPI.array2field(FT.(face_sw_flux_dn), axes(Y.f)),
            sw_flux_up = RRTMGPI.array2field(FT.(face_sw_flux_up), axes(Y.f)),
        )
        if atmos.radiation_mode isa
           RRTMGPI.AllSkyRadiationWithClearSkyDiagnostics
            (;
                face_clear_lw_flux_dn,
                face_clear_lw_flux_up,
                face_clear_sw_flux_dn,
                face_clear_sw_flux_up,
            ) = p.radiation_model
            rad_clear_diagnostic = (;
                clear_lw_flux_down = RRTMGPI.array2field(
                    FT.(face_clear_lw_flux_dn),
                    axes(Y.f),
                ),
                clear_lw_flux_up = RRTMGPI.array2field(
                    FT.(face_clear_lw_flux_up),
                    axes(Y.f),
                ),
                clear_sw_flux_down = RRTMGPI.array2field(
                    FT.(face_clear_sw_flux_dn),
                    axes(Y.f),
                ),
                clear_sw_flux_up = RRTMGPI.array2field(
                    FT.(face_clear_sw_flux_up),
                    axes(Y.f),
                ),
            )
        else
            rad_clear_diagnostic = NamedTuple()
        end
    elseif atmos.radiation_mode isa CA.RadiationDYCOMS_RF01
        # TODO: add radiation diagnostics
        rad_diagnostic = NamedTuple()
        rad_clear_diagnostic = NamedTuple()
    elseif atmos.radiation_mode isa CA.RadiationTRMM_LBA
        # TODO: add radiation diagnostics
        rad_diagnostic = NamedTuple()
        rad_clear_diagnostic = NamedTuple()
    else
        rad_diagnostic = NamedTuple()
        rad_clear_diagnostic = NamedTuple()
    end

    diagnostic = merge(
        dry_diagnostic,
        moist_diagnostic,
        vert_diff_diagnostic,
        rad_diagnostic,
        rad_clear_diagnostic,
        turbulence_convection_diagnostic,
    )

    day = floor(Int, t / (60 * 60 * 24))
    sec = floor(Int, t % (60 * 60 * 24))
    @info "Saving diagnostics to HDF5 file on day $day second $sec"
    output_file = joinpath(output_dir, "day$day.$sec.hdf5")
    hdfwriter = InputOutput.HDF5Writer(output_file, comms_ctx)
    InputOutput.HDF5.write_attribute(hdfwriter.file, "time", t) # TODO: a better way to write metadata
    InputOutput.write!(hdfwriter, Y, "Y")
    InputOutput.write!(
        hdfwriter,
        Fields.FieldVector(; pairs(diagnostic)...),
        "diagnostics",
    )
    Base.close(hdfwriter)
    return nothing
end

function save_restart_func(integrator)
    (; t, u, p) = integrator
    (; output_dir) = p.simulation
    Y = u
    day = floor(Int, t / (60 * 60 * 24))
    sec = floor(Int, t % (60 * 60 * 24))
    @info "Saving restart file to HDF5 file on day $day second $sec"
    mkpath(joinpath(output_dir, "restart"))
    output_file = joinpath(output_dir, "restart", "day$day.$sec.hdf5")
    hdfwriter = InputOutput.HDF5Writer(output_file, comms_ctx)
    InputOutput.HDF5.write_attribute(hdfwriter.file, "time", t) # TODO: a better way to write metadata
    InputOutput.write!(hdfwriter, Y, "Y")
    Base.close(hdfwriter)
    return nothing
end

function gc_func(integrator)
    full = true # whether to do a full GC
    num_pre = Base.gc_num()
    alloc_since_last = (num_pre.allocd + num_pre.deferred_alloc) / 2^20
    live_pre = Base.gc_live_bytes() / 2^20
    GC.gc(full)
    live_post = Base.gc_live_bytes() / 2^20
    num_post = Base.gc_num()
    gc_time = (num_post.total_time - num_pre.total_time) / 10^9 # count in ns
    @debug(
        "GC",
        t = integrator.t,
        "alloc since last GC (MB)" = alloc_since_last,
        "live mem pre (MB)" = live_pre,
        "live mem post (MB)" = live_post,
        "GC time (s)" = gc_time,
        "# pause" = num_post.pause,
        "# full_sweep" = num_post.full_sweep,
    )
end
