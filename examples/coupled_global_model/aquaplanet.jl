import ClimaOcean as Ocean
import ClimaAtmos as Atmos
import Oceananigans

using ClimaEarth
using Printf
using GLMakie

####
#### Atmosphere simulation
####

atmosphere_configuration = Atmos.AtmosConfig("simple_atmos_simulation.yml")
atmosphere = Atmos.get_simulation(atmosphere_configuration)
Δt = atmosphere.integrator.dt # set in yml file

####
#### A near-global ocean
####

arch = Oceananigans.CPU()
Nx = 360
Ny = 160
Nz = 30
z_faces = Ocean.exponential_z_faces(; Nz, h=30, depth=6000)

grid = Oceananigans.LatitudeLongitudeGrid(arch;
                                          size = (Nx, Ny, Nz),
                                          longitude = (0, 360),
                                          latitude = (-80, 80),
                                          z = z_faces,
                                          halo = (7, 7, 7))
                  
momentum_advection   = Oceananigans.WENOVectorInvariant()
tracer_advection     = Oceananigans.WENO()
free_surface         = Oceananigans.SplitExplicitFreeSurface(grid; substeps=30)
vertical_mixing      = Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivity()
horizontal_viscosity = Oceananigans.HorizontalScalarDiffusivity(ν=2000)
closure = vertical_mixing #(vertical_mixing, horizontal_viscosity)

ocean = Ocean.ocean_simulation(grid; momentum_advection, tracer_advection,
                               closure, free_surface, warn=false)

# Set up initial conditions for temperature and salinity
Tatm(λ, φ, z=0) = 30 * cosd(φ)
Tᵢ(λ, φ, z) = 30 * (1 - tanh((abs(φ) - 40) / 5)) / 2 + rand()
Sᵢ(λ, φ, z) = 30 - 5e-3 * z + rand()
Oceananigans.set!(ocean.model, T=Tᵢ, S=Sᵢ)

# Ocean.set!(ocean.model, T=Ocean.ECCOMetadata(:temperature),
#                         S=Ocean.ECCOMetadata(:salinity))

#####
##### The coupled model
#####

radiation  = Ocean.Radiation(ocean_albedo=0.03)
sea_ice    = Ocean.FreezingLimitedOceanTemperature()
model      = Ocean.OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)
simulation = Ocean.Simulation(model; Δt, stop_time=360 * Oceananigans.Units.days)

# OceanSeaIceModel interpolates the atmosphere to the ocean grid:
heatmap(model.interfaces.near_surface_atmosphere_state.T)
display(current_figure())

#####
##### Set up some callbaks + diagnostics and run the simulation
#####

wallclock = Ref(time_ns())

function progress(sim)
    uo, vo, wo = sim.model.ocean.model.velocities

    max_uo = maximum(abs, uo)
    max_vo = maximum(abs, vo)
    max_wo = maximum(abs, wo)

    ρτx = sim.model.interfaces.atmosphere_ocean_interface.fluxes.x_momentum
    ρτy = sim.model.interfaces.atmosphere_ocean_interface.fluxes.y_momentum
    ΣQ = sim.model.interfaces.net_fluxes.ocean_surface.Q

    max_ρτx = maximum(abs, ρτx)
    max_ρτy = maximum(abs, ρτy)
    max_ΣQ = maximum(abs, ΣQ)

    elapsed = 1e-9 * (time_ns() - wallclock[])
    sdpd = sim.Δt / elapsed

    msg = @sprintf("Iter: %d, time: %s, SDPD: %.2f, max|uo|: (%.2e, %.2e, %.2e) (m s⁻¹)",
                   iteration(sim), prettytime(sim), sdpd, max_uo, max_vo, max_wo)

    msg *= @sprintf(", max(ρτ): (%.2e, %.2e) N m⁻², max(ΣQ): %.2e W m⁻²",
                    max_ρτx, max_ρτy, max_ΣQ)

    @info msg

    wallclock[] = time_ns()

    return nothing
end

Oceananigans.add_callback!(simulation, progress, Oceananigans.IterationInterval(10))

using Oceananigans: ∂x, ∂y
uo, vo, wo = ocean.model.velocities
T = ocean.model.tracers.T
S = ocean.model.tracers.S
ζ = ∂x(vo) - ∂y(uo)
outputs = (; uo, vo, wo, ζ, T, S)
surface_writer = Oceananigans.JLD2OutputWriter(ocean.model, outputs,
                                               schedule = IterationInterval(10),
                                               filename = "aquaplanet_ocean.jld2",
                                               indices = (:, :, 30),
                                               overwrite_existing = true)  
simulation.output_writers[:ocean] = surface_writer

run!(simulation)

