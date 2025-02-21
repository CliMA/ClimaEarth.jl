import ClimaOcean
import OrthogonalSphericalShellGrids

using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using Printf

####
#### A near-global ocean
####

arch = Oceananigans.CPU()
Nx = 180 # longitudinal direction -> 250 points is about 1.5ᵒ resolution
Ny = 85  # meridional direction -> same thing, 48 points is about 1.5ᵒ resolution
Nz = 40
z_faces = ClimaOcean.exponential_z_faces(; Nz, h=30, depth=4000)
# z_faces = ZStarVerticalCoordinate(z_faces)

grid = OrthogonalSphericalShellGrids.TripolarGrid(arch; size=(Nx, Ny, Nz), z=z_faces, halo=(4, 4, 4))

url = "https://www.dropbox.com/scl/fi/zy1cu64ybn93l67rjgiq0/\
       Downsampled_ETOPO_2022.nc?\
       rlkey=5upqvoxrnljj205amqf663vcw&st=ou8b32tt&dl=0"
filename = isfile("Downsampled_ETOPO_2022.nc") ?
                  "Downsampled_ETOPO_2022.nc" : download(url, "Downsampled_ETOPO_2022.nc")
                  
bottom_height = ClimaOcean.regrid_bathymetry(grid; filename, dir="./")
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height); active_cells_map=true)

momentum_advection = VectorInvariant()
tracer_advection = WENO()
free_surface = SplitExplicitFreeSurface(grid; substeps=30)
vertical_mixing = CATKEVerticalDiffusivity()
horizontal_viscosity = HorizontalScalarDiffusivity(ν=2000)
closure = (vertical_mixing, horizontal_viscosity)

ocean = ClimaOcean.ocean_simulation(grid;
                                    momentum_advection,
                                    tracer_advection,
                                    closure,
                                    free_surface)

# Set up initial conditions for temperature and salinity
ClimaOcean.set!(ocean.model, T=ClimaOcean.ECCOMetadata(:temperature),
                             S=ClimaOcean.ECCOMetadata(:salinity))

# Add atmopshere here
backend = ClimaOcean.JRA55NetCDFBackend(41)
atmosphere = ClimaOcean.JRA55PrescribedAtmosphere(arch; backend)

# Coupled simulation
Δt = 30minutes
radiation = ClimaOcean.Radiation(ocean_albedo=0.03)
similarity_theory = ClimaOcean.SimilarityTheoryTurbulentFluxes(grid; maxiter=10)
sea_ice = ClimaOcean.FreezingLimitedOceanTemperature()

model = ClimaOcean.OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation, similarity_theory)
simulation = ClimaOcean.Simulation(model; Δt, stop_iteration=1000) #stop_time=400days)

# Diagnostics
save_interval = Int(6hours/Δt)
u, v, w = ocean.model.velocities
T = ocean.model.tracers.T
S = ocean.model.tracers.S
s = sqrt(u^2 + v^2)
η = ocean.model.free_surface.η

simulation.output_writers[:surface] = JLD2OutputWriter(ocean.model, (; T, S, s),
                                                       schedule = IterationInterval(save_interval),
                                                       indices = (:, :, grid.Nz),
                                                       overwrite_existing = true,
                                                       filename = "surface_fields.jld2")

#=
Q  = simulation.model.interfaces.total.ocean.heat
τx = simulation.model.interfaces.total.ocean.momentum.u
τy = simulation.model.interfaces.total.ocean.momentum.v
Fv = simulation.model.interfaces.total.ocean.tracers.S

simulation.output_writers[:fluxes] = JLD2OutputWriter(ocean.model, (; Q, τx, τy, Fv),
                                                      schedule = IterationInterval(save_interval),
                                                      overwrite_existing = true,
                                                      filename = "surface_fluxes.jld2")
=#

# Also, we add a callback to print a message about how the simulation is going
wall_time = Ref(time_ns())

function progress(sim)
    clock = sim.model.clock
    atmos = sim.model.atmosphere
    elapsed = 1e-9 * (time_ns() - wall_time[])

    ocean_model = sim.model.ocean.model
    max_u = maximum(abs, ocean_model.velocities.u)
    max_v = maximum(abs, ocean_model.velocities.v)
    max_T = maximum(ocean_model.tracers.T)
    min_S = minimum(ocean_model.tracers.S)

    @info @sprintf("Iteration: %d, time: %s, wall_time: %s, \
                    max|uh|: (%.2e %.2e) max(T): %.2e, min(S): %.2e\n",
                   clock.iteration, prettytime(clock.time),
                   prettytime(time), prettytime(elapsed),
                   max_u, max_v, max_T, min_S)
    wall_time[] = time_ns()
    return nothing
end

add_callback!(simulation, progress, IterationInterval(10))

ClimaOcean.run!(simulation)

