using Plots
using Statistics

val2string(vv::Val) = string(vv)[6:end-3]

function plot_time_mean(var_name; day_range = collect(50:1:100), nc_dir = ".", level = 1, fig_dir = ".")

    file_paths = filter(endswith(".nc"), readdir(nc_dir; join = true))

    var_array = Array([])
    ct = 0
    for day in day_range
        ncfile = filter(x -> contains(x, "day"*string(day)*"."), file_paths)
        if !isempty(ncfile)
            ct += 1
            nc = NCDataset(ncfile[1], "r")
            global lon = nc["lon"][:]
            global lat = nc["lat"][:]
            global z = nc["z"][:]
            var_array = isempty(var_array) ? nc[nc_name(var_name)][:] : cat(var_array, nc[nc_name(var_name)][:], dims=4)
        end
    end
    # lon-lat
    subplots = []
    push!(subplots, Plots.contourf(lon, lat, transpose(mean(var_array , dims = 4)[:,:,1,1]), xlabel = "Longitude (deg E)", ylabel = "Latitude (deg N)"))

    # lat-z
    push!(subplots, Plots.contourf(lat, z, transpose(mean(var_array , dims = (1,4))[1,:,:,1]), xlabel = "Latitude (deg N)", ylabel = "Height (m)"))
    plot(subplots..., title = val2string(var_name) * "["* get_units(var_name) *"]", size = (800,300), left_margin = 10Plots.mm, bottom_margin = 10Plots.mm)
    png(fig_dir*"/time_mean_"*val2string(var_name)*".png")

end

nc_name(::Val{:temperature}) = "temperature"
get_units(::Val{:temperature}) = "K"
nc_name(::Val{:zonal_wind}) = "u"
get_units(::Val{:zonal_wind}) = "m/s"
nc_name(::Val{:meridional_wind}) = "v"
get_units(::Val{:meridional_wind}) = "m/s"
nc_name(::Val{:density}) = "rho"
get_units(::Val{:density}) = "kg/m^3"
get_units(::Val{:mass_streamfunction}) = "kg/s"

function plot_time_mean(n::Val{:mass_streamfunction}; day_range = collect(50:1:100), nc_dir = ".", level = 1, fig_dir = ".")

    file_paths = filter(endswith(".nc"), readdir(nc_dir; join = true))

    v_array = Array([])
    rho_array = Array([])
    ct = 0
    for day in day_range
        ncfile = filter(x -> contains(x, "day"*string(day)*"."), file_paths)
        if !isempty(ncfile)
            ct += 1
            nc = NCDataset(ncfile[1], "r")
            global lon = nc["lon"][:]
            global lat = nc["lat"][:]
            global z = nc["z"][:]
            v_array = isempty(v_array) ? nc[nc_name(Val(:meridional_wind))][:] : cat(v_array, nc[nc_name(Val(:meridional_wind))][:], dims=4)
            rho_array = isempty(rho_array) ? nc[nc_name(Val(:density))][:] : cat(rho_array, nc[nc_name(Val(:density))][:], dims=4)
        end
    end

    dz = zeros(size(z));
    streamfunction = zeros(size(v_array));
    for i in collect(1:1:length(z))
        if i == 1
            dz = 0.5 * (z[2] + z[1]) - z[1]
        elseif i == length(z)
            dz = z[end] - 0.5 * (z[end] + z[end-1])
        else
            dz =  0.5 * (z[i+1] - z[i-1])
        end
        if i == 1
            @. streamfunction[:,:,i,:] = 2π * FT(6371e3) .* cos.(lat)[1,:,1,1] .* v_array[:,:,i,:] .* rho_array[:,:,i,:]*dz
        else
            @. streamfunction[:,:,i,:] = streamfunction[:,:,i-1,:] + 2π * FT(6371e3) .* cos.(lat)[1,:,1,1] .* v_array[:,:,i,:] .* rho_array[:,:,i,:]*dz
        end
    end

    # lon-lat
    subplots = []
    push!(subplots, Plots.contourf(lon, lat, transpose(mean(streamfunction , dims = 4)[:,:,1,1]), xlabel = "Longitude (deg E)", ylabel = "Latitude (deg N)"))

    # lat-z
    push!(subplots, Plots.contourf(lat, z, transpose(mean(streamfunction , dims = (1,4))[1,:,:,1]), xlabel = "Latitude (deg N)", ylabel = "Height (m)", ylims = (0,25e3)))
    plot(subplots..., title = val2string(n) * "["* get_units(n) *"]", size = (800,300), left_margin = 10Plots.mm, bottom_margin = 10Plots.mm)
    png(fig_dir*"/time_mean_"*val2string(n)*".png")

end