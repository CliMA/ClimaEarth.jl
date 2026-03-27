using Oceananigans
using GLMakie

ocean_filename = "aquaplanet_ocean.jld2"
atmos_filename = "aquaplanet_atmos.jld2"

ζot = FieldTimeSeries(ocean_filename, "ζo")
Tot = FieldTimeSeries(ocean_filename, "To")

uat = FieldTimeSeries(atmos_filename, "ua")
Tat = FieldTimeSeries(atmos_filename, "Ta")
Nt = length(Tot)

fig = Figure(size=(1300, 800))
axzo = Axis(fig[1, 2], xlabel="Longitude (deg)", ylabel="Latitude (deg)")
axTo = Axis(fig[2, 2], xlabel="Longitude (deg)", ylabel="Latitude (deg)")

axua = Axis(fig[1, 3], xlabel="Longitude (deg)", ylabel="Latitude (deg)")
axTa = Axis(fig[2, 3], xlabel="Longitude (deg)", ylabel="Latitude (deg)")

slider = Slider(fig[3, 2:3], range=1:Nt, startvalue=1)
n = slider.value

ζon = @lift ζot[$n]
Ton = @lift Tot[$n]
Tan = @lift Tat[$n]
uan = @lift uat[$n]

ualim = 40
ζolim = 5e-5

hm = heatmap!(axzo, ζon, colormap=:balance, colorrange=(-ζolim, ζolim))
Colorbar(fig[1, 1], hm, flipaxis=false, label="Ocean vorticity (s⁻¹)")

hm = heatmap!(axTo, Ton, colormap=:magma, colorrange=(0, 32))
Colorbar(fig[2, 1], hm, flipaxis=true, label="Ocean temperature (ᵒC)")

hm = heatmap!(axTa, Tan, colormap=:magma, colorrange=(270, 310))
Colorbar(fig[2, 4], hm, flipaxis=false, label="Atmosphere temperature (K)")

hm = heatmap!(axua, uan, colormap=:balance, colorrange=(-ualim, ualim))
Colorbar(fig[1, 4], hm, flipaxis=true, label="Atmosphere zonal velocity (m s⁻¹)")

title = @lift string("Coupled aquapanet after ", prettytime(Tot.times[$n]))
Label(fig[0, 2:3], title)

record(fig, "aquaplanet.mp4", 1:Nt, framerate=12) do nn
    @info "Drawing frame $nn of $Nt"
    n[] = nn
end

display(fig)
