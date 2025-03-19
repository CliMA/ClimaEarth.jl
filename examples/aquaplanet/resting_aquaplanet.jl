import ClimaOcean as Ocean
import ClimaAtmos as Atmos
import Oceananigans
import Thermodynamics

using ClimaEarth
using Printf
# using GLMakie

output_prefix = "resting_aquaplanet"

# Config docs:
# https://clima.github.io/ClimaAtmos.jl/dev/config/
config_dict = Dict(
    #"device" => "auto", #"CUDADevice",
    "device" => "CUDADevice",
    #"z_max" => 60000.0,
    "z_max" => 30000.0,
    "z_elem" => 31,
    "h_elem" => 18, # h_elem = 30 => ~ 1 degree
    "dz_bottom" => 50.0,
    #"dt" => "100secs",
    "dt" => "200secs",
    "topography" => "NoWarp", # "Earth"
    "rayleigh_sponge" => true,
    "implicit_diffusion" => true,
    "approximate_linear_solve_iters" => 2,
    "moist" =>  "equil",
    "surface_setup" => "PrescribedSurface", #DefaultMoninObukhov"
    #"initial_condition" => "MoistBaroclinicWave",
    "initial_condition" => "DecayingProfile",
    #"initial_condition" => "DryBaroclinicWave",
    "vert_diff" => "DecayWithHeightDiffusion",
    "precip_model" =>  "0M",
    "cloud_model" =>  "grid_scale",
    "log_progress" =>  false,
    "output_dir" => output_prefix, 
    "dt_save_to_sol" => "Inf",
    "toml" => ["sphere_aquaplanet.toml"]
)

atmosphere = Atmos.AtmosSimulation(config_dict)
Δt = atmosphere.integrator.dt # set in yml file

####
#### A near-global ocean
####

arch = if Atmos.ClimaComms.device(atmosphere) isa Atmos.ClimaComms.CUDADevice
    Oceananigans.GPU()
else
    Oceananigans.CPU()
end
    
Nx = 360
Ny = 170
Nz = 30
z_faces = Ocean.exponential_z_faces(; Nz, h=30, depth=6000)

grid = Oceananigans.LatitudeLongitudeGrid(arch;
                                          size = (Nx, Ny, Nz),
                                          longitude = (0, 360),
                                          latitude = (-85, 85),
                                          z = z_faces,
                                          halo = (7, 7, 7))
                  
momentum_advection   = Oceananigans.WENOVectorInvariant()
tracer_advection     = Oceananigans.WENO()
free_surface         = Oceananigans.SplitExplicitFreeSurface(grid; substeps=30)
vertical_mixing      = Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivity()
horizontal_viscosity = Oceananigans.HorizontalScalarDiffusivity(ν=2000)
closure = vertical_mixing #(vertical_mixing, horizontal_viscosity)

ocean = Ocean.ocean_simulation(grid; momentum_advection, tracer_advection,
                               closure, free_surface)

# Set up initial conditions for temperature and salinity
Tᵢ(λ, φ, z) = 30 #* (1 - tanh((abs(φ) - 40) / 5)) / 2 + rand()
Sᵢ(λ, φ, z) = 30 - 5e-3 * z #+ 1e-2 * rand()
Oceananigans.set!(ocean.model, T=Tᵢ, S=Sᵢ)

# Ocean.set!(ocean.model, T=Ocean.ECCOMetadata(:temperature),
#                         S=Ocean.ECCOMetadata(:salinity))

#####
##### The coupled model
#####

radiation  = Ocean.Radiation(ocean_albedo=0.03)
sea_ice    = Ocean.FreezingLimitedOceanTemperature()
model      = Ocean.OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)
simulation = Ocean.Simulation(model; Δt, stop_time=60 * Oceananigans.Units.days)

#####
##### Set up some callbaks + diagnostics and run the simulation
#####

wallclock = Ref(time_ns())

function progress(sim)
    atmosphere = sim.model.atmosphere
    ocean = sim.model.ocean
    tsa = atmosphere.integrator.p.precomputed.ᶜts
    elapsed = 1e-9 * (time_ns() - wallclock[])
    sdpd = sim.Δt / elapsed
    uo, vo, wo = ocean.model.velocities

    max_uo = maximum(abs, uo)
    max_vo = maximum(abs, vo)
    max_wo = maximum(abs, wo)

    msg = @sprintf("Iter: %d, time: %s, SDPD: %.1f, max|uo|: (%.2e, %.2e, %.2e) (m s⁻¹)",
                   Oceananigans.iteration(sim),
                   Oceananigans.prettytime(sim),
                   sdpd, max_uo, max_vo, max_wo)

    ua = sim.model.interfaces.exchanger.exchange_atmosphere_state.u
    va = sim.model.interfaces.exchanger.exchange_atmosphere_state.v
    max_ua = maximum(abs, ua)
    max_va = maximum(abs, va)

    # Ta = sim.model.interfaces.exchanger.exchange_atmosphere_state.T
    # min_Ta = minimum(Ta)
    # max_Ta = maximum(Ta)

    ℂa = atmosphere.integrator.p.params.thermodynamics_params
    
    Ta = Thermodynamics.air_temperature.(ℂa, tsa)
    min_Ta, max_Ta = minimum(Ta), maximum(Ta)

    qa = sim.model.interfaces.exchanger.exchange_atmosphere_state.q
    max_qa = maximum(qa)

    msg *= @sprintf("\n   max|ua₀|: (%.1e, %.1e) (m s⁻¹), extrema(Ta): (%d, %d) K, max(qa): %.1e",
                    max_ua, max_va, min_Ta, max_Ta, max_qa) 

    pa = Thermodynamics.air_pressure.(ℂa, tsa)
    ρa = Thermodynamics.air_density.(ℂa, tsa)
    min_ρa, max_ρa = minimum(ρa), maximum(ρa)
    min_pa, max_pa = minimum(pa), maximum(pa)

    msg *= @sprintf("\n   extrema(pa): (%d, %d) Pa, extrema(ρa): (%.6f, %.6f) kg m⁻³",
                    min_pa, max_pa, min_ρa, max_ρa)

    ρτx = sim.model.interfaces.atmosphere_ocean_interface.fluxes.x_momentum
    ρτy = sim.model.interfaces.atmosphere_ocean_interface.fluxes.y_momentum
    Qv = sim.model.interfaces.atmosphere_ocean_interface.fluxes.latent_heat
    Qc = sim.model.interfaces.atmosphere_ocean_interface.fluxes.sensible_heat

    max_ρτxa = maximum(abs, atmosphere.integrator.p.precomputed.sfc_conditions.ρ_flux_uₕ.components.data.:1)
    max_ρτya = maximum(abs, atmosphere.integrator.p.precomputed.sfc_conditions.ρ_flux_uₕ.components.data.:2)

    ΣQ = sim.model.interfaces.net_fluxes.ocean_surface.Q

    max_ρτxo = maximum(abs, ρτx)
    max_ρτyo = maximum(abs, ρτy)
    max_ΣQ = maximum(abs, ΣQ)
    max_Qc = maximum(abs, Qc)
    max_Qv = maximum(abs, Qv)

    msg *= @sprintf("\n   max(ρτo): (%.2f, %.2f) N m⁻², max(ρτa): (%.2f, %.2f) N m⁻², max(ΣQ): %d W m⁻²",
                    max_ρτxo, max_ρτyo,
                    max_ρτxa, max_ρτya, max_ΣQ)

    msg *= @sprintf(", max(Qc): %d W m⁻², max(Qv): %d W m⁻²",
                    max_Qc, max_Qv)

    @info msg

    wallclock[] = time_ns()

    return nothing
end

Oceananigans.add_callback!(simulation, progress, Oceananigans.IterationInterval(10))

using Oceananigans: ∂x, ∂y

ua = model.interfaces.exchanger.exchange_atmosphere_state.u
va = model.interfaces.exchanger.exchange_atmosphere_state.v
Ta = model.interfaces.exchanger.exchange_atmosphere_state.T
qa = model.interfaces.exchanger.exchange_atmosphere_state.q

uo, vo, wo = ocean.model.velocities
To = ocean.model.tracers.T
So = ocean.model.tracers.S
ζo = ∂x(vo) - ∂y(uo)
#ζa = ∂x(va) - ∂y(ua)
#outputs = (; uo, vo, wo, ζo, To, So, ua, va, ζa, Ta, qa)

outputs = (; uo, vo, wo, ζo, To, So) #, ua, va, Ta, qa)
ocean_writer = Oceananigans.JLD2Writer(ocean.model, outputs,
                                       schedule = Oceananigans.IterationInterval(108),
                                       filename = output_prefix * "_ocean.jld2",
                                       indices = (:, :, 30),
                                       overwrite_existing = true)  

outputs = (; ua, va, Ta, qa)
atmos_writer = Oceananigans.JLD2Writer(ocean.model, outputs,
                                       schedule = Oceananigans.IterationInterval(108),
                                       filename = output_prefix * "_atmos.jld2",
                                       overwrite_existing = true)  

simulation.output_writers[:ocean] = ocean_writer
simulation.output_writers[:atmos] = atmos_writer

Oceananigans.run!(simulation)

