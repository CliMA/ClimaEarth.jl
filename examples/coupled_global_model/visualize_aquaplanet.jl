using Oceananigans
using GLMakie

filename = "aquaplanet_ocean.jld2"

ζt = FieldTimeSeries(filename, "ζ")
Tt = FieldTimeSeries(filename, "T")
Nt = length(Tt)

fig = Figure()
axz = Axis(fig[1, 1])
axT = Axis(fig[2, 1])
slider = Slider(fig[3, 1], range=1:Nt, startvalue=1)
n = slider.value

ζn = @lift ζt[$n]
Tn = @lift Tt[$n]

heatmap!(axz, ζn)
heatmap!(axT, Tn)

display(fig)
