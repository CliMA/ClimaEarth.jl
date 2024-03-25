using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using GLMakie

# Configuration
Nz                         = 100
Lz                         = 400
arch                       = CPU()
float_type                 = Float64
gravitational_acceleration = g = 9.81
ocean_reference_density    = ρᵣ = 1024
ocean_heat_capactity       = cₚ = 3991
coriolis_parameter         = 1e-4
temperature_gradient       = 1e-5 / (2e-4 * g)
salinity_gradient          = 1e-8 / (4e-8 * g)
surface_temperature        = 20
surface_salinity           = 35
heat_flux                  = 200 # W m⁻²
temperature_flux           = heat_flux / (ρᵣ * cₚ)

# Grid
grid = RectilinearGrid(arch, float_type,
                       size = Nz,
                       z = (-Lz, 0),
                       topology = (Flat, Flat, Bounded))

# Turbulence closure
closure = CATKEVerticalDiffusivity()

# Boundary conditions
JT = Field{Center, Center, Nothing}(grid)
JS = Field{Center, Center, Nothing}(grid)
Ju = Field{Face, Center, Nothing}(grid)
Jv = Field{Center, Face, Nothing}(grid)

top_T_bc = FluxBoundaryCondition(JT)
top_S_bc = FluxBoundaryCondition(JS)
top_u_bc = FluxBoundaryCondition(Ju)
top_v_bc = FluxBoundaryCondition(Jv)

T_bcs = FieldBoundaryConditions(top=top_T_bc)
S_bcs = FieldBoundaryConditions(top=top_S_bc)
u_bcs = FieldBoundaryConditions(top=top_u_bc)
v_bcs = FieldBoundaryConditions(top=top_v_bc)

# Buoyancy
#equation_of_state = TEOS10EquationOfState(float_type, reference_density=ocean_reference_density)
equation_of_state = LinearEquationOfState(float_type; thermal_expansion=2e-4, haline_contraction=8e-5)
buoyancy = SeawaterBuoyancy(; equation_of_state)

# Coriolis force
coriolis = FPlane(float_type, f=coriolis_parameter)

model = HydrostaticFreeSurfaceModel(; grid, closure, buoyancy,
                                    momentum_advection = nothing,
                                    tracer_advection = nothing,
                                    tracers = (:T, :S, :e),
                                    boundary_conditions = (T=T_bcs, S=S_bcs, u=u_bcs, v=v_bcs))

T₀ = surface_temperature
S₀ = surface_salinity
dTdz = temperature_gradient
dSdz = salinity_gradient
Tᵢ(z) = T₀ + dTdz * z
Sᵢ(z) = S₀ + dSdz * z
set!(model, T=Tᵢ, S=Sᵢ, e=1e-6)

simulation = Simulation(model, Δt=1minutes, stop_iteration=10) #stop_time=1day)

function progress(sim)
    msg = string("Iter: ", iteration(sim), ", time: ", prettytime(sim))
    @info msg
end

add_callback!(simulation, progress, IterationInterval(10))

set!(JT, temperature_flux)

run!(simulation)

u, v, w = model.velocities
T = model.tracers.T
S = model.tracers.S
e = model.tracers.e
z = znodes(T)

fig = Figure(size=(1200, 300))
axT = Axis(fig[1, 1], xlabel="Temperature (ᵒC)", ylabel="z (m)")
axS = Axis(fig[1, 2], xlabel="Salinity (g kg⁻¹)", ylabel="z (m)")
axe = Axis(fig[1, 3], xlabel="Turbulent kinetic energy (m² s⁻²)", ylabel="z (m)")
axu = Axis(fig[1, 4], xlabel="Velocities (m s⁻¹)", ylabel="z (m)")

Tn = interior(T, 1, 1, :)
Sn = interior(S, 1, 1, :)
en = interior(e, 1, 1, :)
un = interior(u, 1, 1, :)
vn = interior(v, 1, 1, :)

lines!(axT, Tn, z, color=:black)
lines!(axS, Sn, z, color=:black)
lines!(axe, en, z, color=:black)
lines!(axu, un, z, color=:black)
lines!(axu, vn, z, color=:black, linestyle=:dash)

display(fig)

