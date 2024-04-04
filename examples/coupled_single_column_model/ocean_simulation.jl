
using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using GLMakie

function default_ocean_grid(arch, float_type)
    return RectilinearGrid(arch, float_type,
                           size = 100,
                           z = (-400, 0),
                           topology = (Flat, Flat, Bounded))
end

function ocean_single_column_simulation(float_type=Float32;
    architecture               = CPU(),
    grid                       = default_ocean_grid(architecture, float_type),
    gravitational_acceleration = 9.81,
    closure                    = CATKEVerticalDiffusivity(float_type),
    reference_density          = 1024,
    coriolis_parameter         = 1e-4)

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
    equation_of_state = TEOS10EquationOfState(float_type; reference_density)
    buoyancy = SeawaterBuoyancy(; equation_of_state, gravitational_acceleration)

    # Coriolis force
    coriolis = FPlane(float_type, f=coriolis_parameter)

    model = HydrostaticFreeSurfaceModel(; grid, closure, buoyancy, coriolis,
                                        momentum_advection = nothing,
                                        tracer_advection = nothing,
                                        tracers = (:T, :S, :e),
                                        boundary_conditions = (T=T_bcs, S=S_bcs, u=u_bcs, v=v_bcs))

    simulation = Simulation(model, Δt=1minutes)
    
    return simulation
end

#=
# T₀ = surface_temperature
    S₀ = surface_salinity
    dTdz = temperature_gradient
    dSdz = salinity_gradient
    Tᵢ(z) = T₀ + dTdz * z
    Sᵢ(z) = S₀ + dSdz * z
    set!(model, T=Tᵢ, S=Sᵢ, e=1e-6)


    display(fig)
temperature_gradient       = 1e-5 / (2e-4 * g),
    salinity_gradient          = 1e-8 / (4e-8 * g),
    surface_temperature        = 20
    surface_salinity           = 35
    heat_flux                  = 200 # W m⁻²
    temperature_flux           = heat_flux / (ρᵣ * cₚ)
=#
