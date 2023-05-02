using Plots
using Statistics

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
            var_array = isempty(var_array) ? nc[nc_name(Val(var_name))][:] : cat(var_array, nc[nc_name(Val(var_name))][:], dims=4)
        end
    end
    # lon-lat
    subplots = []
    push!(subplots, Plots.contourf(lon, lat, transpose(mean(var_array , dims = 4)[:,:,1,1]), xlabel = "Longitude (deg E)", ylabel = "Latitude (deg N)"))

    # lat-z
    push!(subplots, Plots.contourf(lat, z, transpose(mean(var_array , dims = (1,4))[1,:,:,1]), xlabel = "Latitude (deg N)", ylabel = "Height (m)"))
    plot(subplots..., title = string(var_name) * "["* get_units(Val(var_name)) *"]", size = (800,300), left_margin = 10Plots.mm, bottom_margin = 10Plots.mm)
    png(fig_dir*"/time_mean_"*string(var_name)*".png")

end

nc_name(::Val{:temperature}) = "temperature"
get_units(::Val{:temperature}) = "K"
nc_name(::Val{:zonal_wind}) = "u"
get_units(::Val{:zonal_wind}) = "m/s"
nc_name(::Val{:meridional_wind}) = "v"
get_units(::Val{:zonal_wind}) = "m/s"
