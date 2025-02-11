
function remaining_tendency!(Yₜ, Y, p, t)
    p.test_dycore_consistency && fill_with_nans!(p)
    @nvtx "remaining tendency" color = colorant"yellow" begin
        Yₜ .= zero(eltype(Yₜ))
        @nvtx "precomputed quantities" color = colorant"orange" begin
            set_precomputed_quantities!(Y, p, t)
        end
        @nvtx "horizontal" color = colorant"orange" begin
            horizontal_advection_tendency!(Yₜ, Y, p, t)
            @nvtx "hyperdiffusion tendency" color = colorant"yellow" begin
                hyperdiffusion_tendency!(Yₜ, Y, p, t)
            end
        end
        @nvtx "vertical" color = colorant"orange" begin
            explicit_vertical_advection_tendency!(Yₜ, Y, p, t)
        end
        @nvtx "additional_tendency!" color = colorant"orange" begin
            additional_tendency!(Yₜ, Y, p, t)
        end
        @nvtx "dss_remaining_tendency" color = colorant"blue" begin
            dss!(Yₜ, p, t)
        end
    end
    return Yₜ
end

function additional_tendency!(Yₜ, Y, p, t)
    viscous_sponge_tendency!(Yₜ, Y, p, t, p.atmos.viscous_sponge)

    # Vertical tendencies
    Fields.bycolumn(axes(Y.c)) do colidx
        rayleigh_sponge_tendency!(Yₜ, Y, p, t, colidx, p.atmos.rayleigh_sponge)
        forcing_tendency!(Yₜ, Y, p, t, colidx, p.forcing_type)
        subsidence_tendency!(Yₜ, Y, p, t, colidx, p.subsidence)
        edmf_coriolis_tendency!(Yₜ, Y, p, t, colidx, p.edmf_coriolis)
        large_scale_advection_tendency!(Yₜ, Y, p, t, colidx, p.ls_adv)

        (; vert_diff) = p.atmos
        if p.atmos.coupling isa Decoupled
            get_surface_fluxes!(Y, p, t, colidx, vert_diff)
        end
        vertical_diffusion_boundary_layer_tendency!(
            Yₜ,
            Y,
            p,
            t,
            colidx,
            vert_diff,
        )

        radiation_tendency!(Yₜ, Y, p, t, colidx, p.radiation_model)
        explicit_sgs_flux_tendency!(Yₜ, Y, p, t, colidx, p.turbconv_model)
        precipitation_tendency!(Yₜ, Y, p, t, colidx, p.precip_model)
    end
    # TODO: make bycolumn-able
    non_orographic_gravity_wave_tendency!(
        Yₜ,
        Y,
        p,
        t,
        p.atmos.non_orographic_gravity_wave,
    )
    orographic_gravity_wave_tendency!(
        Yₜ,
        Y,
        p,
        t,
        p.atmos.orographic_gravity_wave,
    )
end
