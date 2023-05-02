# plot raw model data (uses a low-order interpolation)
using Plots
using ClimaCorePlots
using ClimaCore: Fields, Geometry
using Statistics: std

solution =  integrator.sol;

function png_end(sol, var_name::Symbol; output_dir = ".", clims = nothing, level = 1)
    var = get(sol.u[end], Val(var_name))
    clims = isnothing(clims) ? map(x -> x * std( parent(var)) * 2 , (-1,1) )  : clims
    Plots.plot(var, clims = clims, title = string(var_name), level = level)
    Plots.png(joinpath(output_dir, "plot_"*string(var_name)*".png"))
end

function mp4(sol, var_name::Symbol; output_dir = ".", clims = nothing, level = 1)
    clims = isnothing(clims) ? map(x -> x * std( parent(get(sol.u[end], Val(var_name)))) * 2 , (-1,1) )  : clims
    anim = Plots.@animate for (i, u) in enumerate(sol.u) # note: sol.u = solution vector here, not zonal wind
        Plots.plot(get(u, Val(var_name)), clims = clims, title = string(var_name) *" "* string(sol.t[i]/86400), level = level, right_margin = 10Plots.mm)
    end
    Plots.mp4(anim, joinpath(output_dir, "anim_"*string(var_name)*".mp4"), fps = 10)
end
get(u, ::Val{:zonal_wind}) = Geometry.UVVector.(u.c.uₕ).components.data.:1
get(u, ::Val{:meridional_wind}) = Geometry.UVVector.(u.c.uₕ).components.data.:2
get(u, ::Val{:total_energy}) = u.c.ρe_tot

get_coordinates(integrator) = Fields.local_geometry_field(integrator.u.c).coordinates
