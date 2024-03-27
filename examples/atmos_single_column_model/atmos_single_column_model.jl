using ClimaAtmos
using ClimaCore
using Thermodynamics
using ClimaCore.Fields: level
using ClimaCore.Geometry: UVVector
using ClimaCore.Geometry: Covariant12Vector
using SurfaceFluxes
using StaticArrays: SVector
using Thermodynamics: Liquid
using Thermodynamics: PhaseNonEquil_ρTq, PhasePartition

saturation_specific_humidity = Thermodynamics.q_vap_saturation_generic

function surface_value(f)
    f1 = level(f, 1) # view into `f` at surface (level 1)
    return first(parent(f1)) # surface value
end

config_filename = "config.yml"
toml_filename = "toml/single_column_precipitation_test.toml"

yml_config = """
config: "column"
initial_condition: "SimplePlume"
surface_setup: "DefaultMoninObukhov"
toml: [$(toml_filename)]
z_elem: 200
z_max: 10000.0
z_stretch: false
moist: "nonequil"
dt_save_to_sol: "Inf"
precip_model: "1M"
vert_diff: "FriersonDiffusion"
implicit_diffusion: true
t_end: "10days"
z_elem: 100
t_end: "10days"
dt_save_to_sol: "Inf"
dt: "15secs"
job_id: "test"
hyperdiff": nothing
"""

#toml: [$(toml_filename)]

open(config_filename, "w") do io
    write(io, yml_config)
end

config = ClimaAtmos.AtmosConfig(config_filename)
simulation = ClimaAtmos.get_simulation(config)

# TODO: figure out how to set initial conditions
# center_space = axes(simulation.integrator.u.c)
# set!(simulation, rho=1.2, T=DryAdiabaticProfile())
# u = 1
Ûʰ = simulation.integrator.u.c.uₕ # covariant velocity
Uʰ = UVVector.(Ûʰ) # contravariant velocity
Uʰ.components.data.:1 .= 1 # set u = 1
Uʰ.components.data.:2 .= 0 # set v = 0
Ûʰ .= Covariant12Vector.(Uʰ)

# Extract atmospheric thermodynamic state
Ûʰ₁ = level(Ûʰ, 1)
Uʰ₁ = ClimaCore.Geometry.UVVector.(Ûʰ₁)
u₁ = parent(Uʰ₁)[1]
v₁ = parent(Uʰ₁)[2]

ℂ = ClimaAtmos.Parameters.thermodynamics_params(simulation.integrator.p.params) 
𝒯 = simulation.integrator.p.precomputed.ᶜts
𝒯₁ = level(𝒯, 1)
ρ₁ = first(parent(𝒯₁.ρ))
e₁ = first(parent(𝒯₁.e_int))
q₁ = first(parent(𝒯₁.q))
𝒯₁ = PhaseEquil_ρeq(ℂ, ρ₁, e₁, q₁)

# Compute ocean thermodynamic state, including saturation specific humidity
FT = typeof(u₁)
u₀ = convert(FT, 0)
v₀ = convert(FT, 0)
z₀ = convert(FT, 0)
T₀ = convert(FT, 293.0)
x_H₂O = convert(FT, 0.98) # mole fraction of H₂O in seawater

cvₘ = Thermodynamics.cv_m(ℂ, 𝒯₁)
R₁ = Thermodynamics.gas_constant_air(ℂ, 𝒯₁)
κ₁ = cvₘ / R₁
ρ₁ = Thermodynamics.air_density(ℂ, 𝒯₁)
T₁ = Thermodynamics.air_temperature(ℂ, 𝒯₁)
ρ₀ = ρ₁ * (T₀ / T₁)^κ₁

q★ = x_H₂O * saturation_specific_humidity(ℂ, T₀, ρ₀, Liquid())
q₀ = PhasePartition(q★)
𝒯₀ = PhaseNonEquil_ρTq(ℂ, ρ₀, T₀, q₀)

values = SurfaceFluxes.StateValues(z₀, SVector(u₀, v₀), 𝒯₀)

#=
ρe1 = surface_value(simulation.integrator.u.c.T)

# or
Nt = 100
Δt = 60.0

for n = 1:Nt
    step!(simulation, Δt, true)
end
=#



# sol_res = ClimaAtmos.solve_atmos!(simulation)

#=
struct ConstantWinds <: InitialCondition
    u :: Float64
    v :: Float64
end

function (ic::ConstantWinds)(params)
    function local_state(local_geometry)
        FT = eltype(params)
        thermo_params = CAP.thermodynamics_params(params)
        temp_profile = DryAdiabaticProfile{FT}(thermo_params, FT(310), FT(290))

        (; z) = local_geometry.coordinates
        T, p = temp_profile(thermo_params, z)
        q_tot = FT(0)
        u = convert(FT, ic.u)
        v = convert(FT, ic.v)

        return LocalState(;
            params,
            geometry = local_geometry,
            thermo_state = TD.PhaseEquil_pTq(thermo_params, p, T, q_tot),
            velocity = Geometry.UVVector(u, v)
        )
    end
    return local_state
end
=#
