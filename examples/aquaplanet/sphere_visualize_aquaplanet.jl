using Oceananigans
using GLMakie
using Printf

ocean_filename = "aquaplanet_ocean.jld2"
atmos_filename = "aquaplanet_atmos.jld2"

function geographic2cartesian(λ, φ, r=1)
    Nλ = length(λ)
    Nφ = length(φ)
    λ = repeat(reshape(λ, Nλ, 1), 1, Nφ)
    φ = repeat(reshape(φ, 1, Nφ), Nλ, 1)

    λ_azimuthal = λ .+ 180  # Convert to λ ∈ [0°, 360°]
    φ_azimuthal = 90 .- φ   # Convert to φ ∈ [0°, 180°] (0° at north pole)

    x = @. r * cosd(λ_azimuthal) * sind(φ_azimuthal)
    y = @. r * sind(λ_azimuthal) * sind(φ_azimuthal)
    z = @. r * cosd(φ_azimuthal)

    return x, y, z
end

ζot = FieldTimeSeries(ocean_filename, "ζo")
Tot = FieldTimeSeries(ocean_filename, "To")
Tot .+= 273.15 # convert ocean SST from ᵒC to ᵒK
uat = FieldTimeSeries(atmos_filename, "ua")
Tat = FieldTimeSeries(atmos_filename, "Ta")

Nt = length(Tat)

grid = uat.grid
λ = λnodes(grid, Center(), Center(), Center())
φ = φnodes(grid, Center(), Center(), Center())
x, y, z = geographic2cartesian(λ, φ)

n = Observable(1)
ζon = @lift interior(ζot[$n], :, :, 1)
uan = @lift interior(uat[$n], :, :, 1)
Tan = @lift interior(Tat[$n], :, :, 1)
Ton = @lift interior(Tot[$n], :, :, 1)

fig = Figure(size=(1400, 1400))
kw = (elevation=0.6, azimuth=5, aspect=:equal)

axua = Axis3(fig[1, 1]; kw...)
axTa = Axis3(fig[1, 2]; kw...)
axζo = Axis3(fig[3, 1]; kw...)
axTo = Axis3(fig[3, 2]; kw...)

sf = surface!(axua, x, y, z, color=uan, colorrange=(-40, 40), colormap=:balance, nan_color=:lightgray)
Colorbar(fig[2, 1], sf, width=Relative(0.6), vertical=false, label = "Atmosphere zonal velocity (m s⁻¹)", labelsize=20)

sf = surface!(axTa, x, y, z, color=Tan, colorrange=(270, 310), colormap=:magma, nan_color=:lightgray)
Colorbar(fig[2, 2], sf, width=Relative(0.6), vertical=false, label = "Atmosphere surface temperature (ᵒK)", labelsize=20)

sf = surface!(axζo, x, y, z, color=ζon, colorrange=(-5e-5, 5e-5), colormap=:balance, nan_color=:lightgray)
Colorbar(fig[4, 1], sf, width=Relative(0.6), vertical=false, label = "Ocean vorticity (s⁻¹)", labelsize=20)

sf = surface!(axTo, x, y, z, color=Ton, colorrange=(270, 305), colormap=:magma, nan_color=:lightgray)
Colorbar(fig[4, 2], sf, width=Relative(0.6), vertical=false, label = "Ocean SST (ᵒK)", labelsize=20)

t = uat.times
title = @lift @sprintf("%1.2f days", t[$n] / Oceananigans.Units.day)
Label(fig[0, :], title, fontsize=24)

for ax in (axua, axTa, axζo, axTo)
    hidedecorations!(ax)
    hidespines!(ax)
end

rowgap!(fig.layout, 1, Relative(-0.1))
rowgap!(fig.layout, 3, Relative(-0.1))

colgap!(fig.layout, 1, Relative(-0.07))
# colgap!(fig.layout, 2, Relative(-0.07))

fig

record(fig, "aquaplanet_sphere.mp4", 1:Nt, framerate = 16) do nn
    @info "Drawing frame $nn of $Nt..."
    n[] = nn
end
