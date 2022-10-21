import ClimaAtmos.Parameters as CAP
using Statistics: mean
using Dierckx: Spline1D
using Dates: Second
import Pkg
import NCDatasets as NC
import RRTMGP
using Insolation: instantaneous_zenith_angle

#####
##### RRTMGP
#####

"""
    rrtmgp_artifact(file_name)

Returns the filename of an artifact stored in
`RRTMGPReferenceData/<file_name>`.
"""
function rrtmgp_artifact(file_name)
    artifact_name = "RRTMGPReferenceData"
    artifacts_file = joinpath(pkgdir(RRTMGP), "test", "Artifacts.toml")
    data_folder = joinpath(
        Pkg.Artifacts.ensure_artifact_installed(artifact_name, artifacts_file),
        artifact_name,
    )
    return joinpath(data_folder, file_name)
end

"""
    data_loader(fn, file_name)

Loads data from an `NCDataset` from the `RRTMGP.jl` artifact stored in
`RRTMGPReferenceData/<file_name>`, and calls a function, `fn` on the
dataset.
"""
function data_loader(fn, file_name)
    NC.Dataset(rrtmgp_artifact(file_name), "r") do ds
        fn(ds)
    end
end

function radiation_model_cache(
    Y,
    params,
    radiation_mode::RRTMGPI.AbstractRRTMGPMode = RRTMGPI.ClearSkyRadiation();
    interpolation = RRTMGPI.BestFit(),
    bottom_extrapolation = RRTMGPI.SameAsInterpolation(),
    idealized_insolation = true,
    idealized_h2o = false,
    idealized_clouds = false,
    thermo_dispatcher,
)
    rrtmgp_params = CAP.rrtmgp_params(params)
    thermo_params = CAP.thermodynamics_params(params)
    if idealized_h2o && radiation_mode isa RRTMGPI.GrayRadiation
        error("idealized_h2o can't be used with $radiation_mode")
    end
    if idealized_clouds && (
        radiation_mode isa RRTMGPI.GrayRadiation ||
        radiation_mode isa RRTMGPI.ClearSkyRadiation
    )
        error("idealized_clouds can't be used with $radiation_mode")
    end

    bottom_coords = Fields.coordinate_field(Spaces.level(Y.c, 1))
    if eltype(bottom_coords) <: Geometry.LatLongZPoint
        latitude = RRTMGPI.field2array(bottom_coords.lat)
    else
        latitude = RRTMGPI.field2array(zero(bottom_coords.z)) # flat space is on Equator
    end
    input_data =
        NC.Dataset(rrtmgp_artifact(joinpath("atmos_state", "clearsky_as.nc")))
    if radiation_mode isa RRTMGPI.GrayRadiation
        kwargs = (;
            lapse_rate = 3.5,
            optical_thickness_parameter = (@. 7.2 +
                                              (1.8 - 7.2) * sind(latitude)^2),
        )
    else
        # the pressure and ozone concentrations are provided for each of 100
        # sites, which we average across
        n = input_data.dim["layer"]
        input_center_pressure =
            vec(mean(reshape(input_data["pres_layer"][:, :], n, :); dims = 2))
        # the first values along the third dimension of the ozone concentration
        # data are the present-day values
        input_center_volume_mixing_ratio_o3 =
            vec(mean(reshape(input_data["ozone"][:, :, 1], n, :); dims = 2))

        # interpolate the ozone concentrations to our initial pressures (set the
        # kinetic energy to 0 when computing the pressure using total energy)
        pressure2ozone =
            Spline1D(input_center_pressure, input_center_volume_mixing_ratio_o3)
        ᶜts = CA.thermo_state(Y, thermo_params, thermo_dispatcher, ᶜinterp, 0)
        ᶜp = @. TD.air_pressure(thermo_params, ᶜts)
        center_volume_mixing_ratio_o3 =
            RRTMGPI.field2array(@. FT(pressure2ozone(ᶜp)))

        # the first value for each global mean volume mixing ratio is the
        # present-day value
        input_vmr(name) =
            input_data[name][1] * parse(FT, input_data[name].attrib["units"])
        kwargs = (;
            use_global_means_for_well_mixed_gases = true,
            center_volume_mixing_ratio_h2o = NaN, # initialize in tendency
            center_volume_mixing_ratio_o3,
            volume_mixing_ratio_co2 = input_vmr("carbon_dioxide_GM"),
            volume_mixing_ratio_n2o = input_vmr("nitrous_oxide_GM"),
            volume_mixing_ratio_co = input_vmr("carbon_monoxide_GM"),
            volume_mixing_ratio_ch4 = input_vmr("methane_GM"),
            volume_mixing_ratio_o2 = input_vmr("oxygen_GM"),
            volume_mixing_ratio_n2 = input_vmr("nitrogen_GM"),
            volume_mixing_ratio_ccl4 = input_vmr("carbon_tetrachloride_GM"),
            volume_mixing_ratio_cfc11 = input_vmr("cfc11_GM"),
            volume_mixing_ratio_cfc12 = input_vmr("cfc12_GM"),
            volume_mixing_ratio_cfc22 = input_vmr("hcfc22_GM"),
            volume_mixing_ratio_hfc143a = input_vmr("hfc143a_GM"),
            volume_mixing_ratio_hfc125 = input_vmr("hfc125_GM"),
            volume_mixing_ratio_hfc23 = input_vmr("hfc23_GM"),
            volume_mixing_ratio_hfc32 = input_vmr("hfc32_GM"),
            volume_mixing_ratio_hfc134a = input_vmr("hfc134a_GM"),
            volume_mixing_ratio_cf4 = input_vmr("cf4_GM"),
            volume_mixing_ratio_no2 = 1e-8, # not available in input_data
            latitude,
        )
        if !(radiation_mode isa RRTMGPI.ClearSkyRadiation)
            kwargs = (;
                kwargs...,
                center_cloud_liquid_effective_radius = 12,
                center_cloud_ice_effective_radius = 95,
                ice_roughness = 2,
            )
            ᶜz = Fields.coordinate_field(Y.c).z
            ᶜΔz = Fields.local_geometry_field(Y.c).∂x∂ξ.components.data.:9
            if idealized_clouds # icy cloud on top and wet cloud on bottom
                ᶜis_bottom_cloud = Fields.Field(
                    DataLayouts.replace_basetype(Fields.field_values(ᶜz), Bool),
                    axes(Y.c),
                ) # need to fix several ClimaCore bugs in order to simplify this
                ᶜis_top_cloud = similar(ᶜis_bottom_cloud)
                @. ᶜis_bottom_cloud = ᶜz > 1e3 && ᶜz < 1.5e3
                @. ᶜis_top_cloud = ᶜz > 4e3 && ᶜz < 5e3
                kwargs = (;
                    kwargs...,
                    center_cloud_liquid_water_path = RRTMGPI.field2array(
                        @. ifelse(ᶜis_bottom_cloud, FT(0.002) * ᶜΔz, FT(0))
                    ),
                    center_cloud_ice_water_path = RRTMGPI.field2array(
                        @. ifelse(ᶜis_top_cloud, FT(0.001) * ᶜΔz, FT(0))
                    ),
                    center_cloud_fraction = RRTMGPI.field2array(
                        @. ifelse(
                            ᶜis_bottom_cloud || ᶜis_top_cloud,
                            FT(1),
                            FT(0) * ᶜΔz,
                        )
                    ),
                )
            else
                kwargs = (;
                    kwargs...,
                    center_cloud_liquid_water_path = NaN, # initialized in callback
                    center_cloud_ice_water_path = NaN, # initialized in callback
                    center_cloud_fraction = NaN, # initialized in callback
                )
            end
        end
    end

    if RRTMGPI.requires_z(interpolation) ||
       RRTMGPI.requires_z(bottom_extrapolation)
        kwargs = (;
            kwargs...,
            center_z = RRTMGPI.field2array(Fields.coordinate_field(Y.c).z),
            face_z = RRTMGPI.field2array(Fields.coordinate_field(Y.f).z),
        )
    end

    if idealized_insolation # perpetual equinox with no diurnal cycle
        solar_zenith_angle = FT(π) / 3
        weighted_irradiance =
            @. 1360 * (1 + FT(1.2) / 4 * (1 - 3 * sind(latitude)^2)) /
               (4 * cos(solar_zenith_angle))
    else
        solar_zenith_angle = weighted_irradiance = NaN # initialized in callback
    end

    radiation_model = RRTMGPI.RRTMGPModel(
        rrtmgp_params,
        data_loader,
        Float64,
        RRTMGP.Device.array_type();
        ncol = length(Spaces.all_nodes(axes(Spaces.level(Y.c, 1)))),
        domain_nlay = Spaces.nlevels(axes(Y.c)),
        radiation_mode,
        interpolation,
        bottom_extrapolation,
        add_isothermal_boundary_layer = true,
        center_pressure = NaN, # initialized in callback
        center_temperature = NaN, # initialized in callback
        surface_temperature = NaN, # initialized in callback
        surface_emissivity = 1,
        direct_sw_surface_albedo = 0.38,
        diffuse_sw_surface_albedo = 0.38,
        solar_zenith_angle,
        weighted_irradiance,
        kwargs...,
    )
    close(input_data)
    return (;
        idealized_insolation,
        idealized_h2o,
        idealized_clouds,
        insolation_tuple = similar(Spaces.level(Y.c, 1), Tuple{FT, FT, FT}),
        ᶠradiation_flux = similar(Y.f, Geometry.WVector{FT}),
        radiation_model,
    )
end
function radiation_tendency!(Yₜ, Y, p, t, colidx, ::RRTMGPI.RRTMGPModel)
    (; ᶠradiation_flux) = p
    ᶜdivᵥ = Operators.DivergenceF2C()
    if :ρθ in propertynames(Y.c)
        error("radiation_tendency! not implemented for ρθ")
    elseif :ρe_tot in propertynames(Y.c)
        @. Yₜ.c.ρe_tot[colidx] -= ᶜdivᵥ(ᶠradiation_flux[colidx])
    elseif :ρe_int in propertynames(Y.c)
        @. Yₜ.c.ρe_int[colidx] -= ᶜdivᵥ(ᶠradiation_flux[colidx])
    end
    return nothing
end
function rrtmgp_model_callback!(integrator)
    Y = integrator.u
    p = integrator.p
    t = integrator.t

    (; ᶜK, ᶜts, T_sfc, params, thermo_dispatcher) = p
    (; idealized_insolation, idealized_h2o, idealized_clouds) = p
    (; insolation_tuple, ᶠradiation_flux, radiation_model) = p
    thermo_params = CAP.thermodynamics_params(params)
    insolation_params = CAP.insolation_params(params)

    radiation_model.surface_temperature .= RRTMGPI.field2array(T_sfc)

    ᶜp = RRTMGPI.array2field(radiation_model.center_pressure, axes(Y.c))
    ᶜT = RRTMGPI.array2field(radiation_model.center_temperature, axes(Y.c))
    @. ᶜK = norm_sqr(C123(Y.c.uₕ) + C123(ᶜinterp(Y.f.w))) / 2
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
        current_datetime = p.simulation.start_date + Second(round(Int, t)) # current time
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

#####
##### DYCOMS_RF01 radiation
#####

function radiation_model_cache(
    Y,
    params,
    radiation_mode::CA.RadiationDYCOMS_RF01;
)
    FT = Spaces.undertype(axes(Y.c))
    # TODO: should we distinguish between radiation mode vs model?
    return (;
        ᶠρ = similar(Y.f, FT),
        ᶠf_rad = similar(Y.f, FT),
        ᶠq_tot = similar(Y.f, FT),
        ᶜq_tot = similar(Y.c, FT),
        ᶜq_liq = similar(Y.c, FT),
        ᶜdTdt_rad = similar(Y.c, FT),
        radiation_model = radiation_mode,
    )
end

"""
    radiation_tendency!(Yₜ, Y, p, t, colidx, self::CA.RadiationDYCOMS_RF01)

see eq. 3 in Stevens et. al. 2005 DYCOMS paper
"""
function radiation_tendency!(Yₜ, Y, p, t, colidx, self::CA.RadiationDYCOMS_RF01)
    (; params) = p
    FT = Spaces.undertype(axes(Y.c))
    thermo_params = CAP.thermodynamics_params(params)
    cp_d = CAP.cp_d(params)
    # Unpack
    ᶜρ = Y.c.ρ[colidx]
    ᶜdTdt_rad = p.ᶜdTdt_rad[colidx]
    ᶜts_gm = p.ᶜts[colidx]
    ᶠf_rad = p.ᶠf_rad[colidx]
    ᶠq_tot = p.ᶠq_tot[colidx]
    ᶠρ = p.edmf_cache.aux.face.ρ[colidx]
    ᶜq_tot = p.edmf_cache.aux.cent.q_tot[colidx]
    ᶜq_liq = p.edmf_cache.aux.cent.q_liq[colidx]
    # TODO: unify (ᶜq_tot, ᶜq_liq) with grid-mean
    @. ᶜq_tot = TD.total_specific_humidity(thermo_params, ᶜts_gm)
    @. ᶜq_liq = TD.liquid_specific_humidity(thermo_params, ᶜts_gm)

    ᶜq_tot_surf = Spaces.level(ᶜq_tot, 1)

    zc = Fields.coordinate_field(axes(ᶜρ)).z
    zf = Fields.coordinate_field(axes(ᶠρ)).z
    # find zi (level of 8.0 g/kg isoline of qt)
    # TODO: report bug: zi and ρ_i are not initialized
    zi = FT(0)
    ρ_i = FT(0)
    Ifq = Operators.InterpolateC2F(;
        bottom = Operators.SetValue(ᶜq_tot_surf),
        top = Operators.Extrapolate(),
    )
    @. ᶠq_tot .= Ifq(ᶜq_tot)
    n = length(parent(ᶠρ))
    for k in CCO.PlusHalf(1):CCO.PlusHalf(n)
        ᶠq_tot_lev = Spaces.level(ᶠq_tot, k)
        ᶠq_tot_lev = parent(ᶠq_tot_lev)[1]
        if (ᶠq_tot_lev < 8.0 / 1000)
            # will be used at cell faces
            zi = Spaces.level(zf, k)
            zi = parent(zi)[1]
            ρ_i = Spaces.level(ᶠρ, k)
            ρ_i = parent(ρ_i)[1]
            break
        end
    end

    ρ_z = Dierckx.Spline1D(vec(zc), vec(ᶜρ); k = 1)
    q_liq_z = Dierckx.Spline1D(vec(zc), vec(ᶜq_liq); k = 1)

    integrand(ρq_l, params, z) = params.κ * ρ_z(z) * q_liq_z(z)
    rintegrand(ρq_l, params, z) = -integrand(ρq_l, params, z)

    z_span = (minimum(zf), maximum(zf))
    rz_span = (z_span[2], z_span[1])
    params = (; κ = self.kappa)

    Δz = parent(Fields.dz_field(axes(ᶜρ)))[1]
    rprob = ODE.ODEProblem(rintegrand, 0.0, rz_span, params; dt = Δz)
    rsol = ODE.solve(rprob, ODE.Tsit5(), reltol = 1e-12, abstol = 1e-12)
    q_0 = rsol.(vec(zf))

    prob = ODE.ODEProblem(integrand, 0.0, z_span, params; dt = Δz)
    sol = ODE.solve(prob, ODE.Tsit5(), reltol = 1e-12, abstol = 1e-12)
    q_1 = sol.(vec(zf))
    parent(ᶠf_rad) .= self.F0 .* exp.(-q_0)
    parent(ᶠf_rad) .+= self.F1 .* exp.(-q_1)

    # cooling in free troposphere
    for k in CCO.PlusHalf(1):CCO.PlusHalf(n)
        zf_lev = Spaces.level(zf, k)
        zf_lev = parent(zf_lev)[1]
        ᶠf_rad_lev = Spaces.level(ᶠf_rad, k)
        if zf_lev > zi
            cbrt_z = cbrt(zf_lev - zi)
            @. ᶠf_rad_lev +=
                ρ_i *
                cp_d *
                self.divergence *
                self.alpha_z *
                (cbrt_z^4 / 4 + zi * cbrt_z)
        end
    end

    ∇c = Operators.DivergenceF2C()
    wvec = Geometry.WVector
    @. ᶜdTdt_rad = -∇c(wvec(ᶠf_rad)) / ᶜρ / cp_d

    @. Yₜ.c.ρe_tot[colidx] += ᶜρ * TD.cv_m(thermo_params, ᶜts_gm) * ᶜdTdt_rad

    return
end

#####
##### TRMM_LBA radiation
#####

function radiation_model_cache(Y, params, radiation_mode::CA.RadiationTRMM_LBA;)
    FT = Spaces.undertype(axes(Y.c))
    # TODO: should we distinguish between radiation mode vs model?
    return (; ᶜdTdt_rad = similar(Y.c, FT), radiation_model = radiation_mode)
end
function radiation_tendency!(
    Yₜ,
    Y,
    p,
    t,
    colidx,
    radiation_mode::CA.RadiationTRMM_LBA,
)
    # TODO: get working (need to add cache / function)
    rad = radiation_mode.rad_profile
    thermo_params = CAP.thermodynamics_params(params)
    ᶜdTdt_rad = p.ᶜdTdt_rad[colidx]
    ᶜρ = Y.c.ρ[colidx]
    ᶜts_gm = p.ᶜts[colidx]
    zc = Fields.coordinate_field(axes(ᶜρ)).z
    @. ᶜdTdt_rad = rad(t, zc)
    @. Yₜ.c.ρe_tot[colidx] += ᶜρ * TD.cv_m(thermo_params, ᶜts_gm) * ᶜdTdt_rad
    return nothing
end
