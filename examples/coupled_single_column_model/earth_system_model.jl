using StaticArrays: SVector

using SurfaceFluxes
using SurfaceFluxes: ValuesOnly, surface_conditions

import Thermodynamics
using Thermodynamics: Liquid
using Thermodynamics: PhaseNonEquil_ρTq, PhasePartition 
using Thermodynamics: saturation_vapor_pressure, q_vap_saturation_from_density

using ClimaAtmos.SurfaceConditions: scalar_flux, vector_flux
using ClimaCore.Utilities: half
using Oceananigans.Models: AbstractModel

const saturation_vapor_specific_humidity = Thermodynamics.q_vap_saturation_generic

struct SingleColumnAtmosOceanModel{C, G, A, O, F, P} <: AbstractModel{Nothing}
    clock :: C
    grid :: G
    atmos :: A
    ocean :: O
    surface :: F
    parameters :: P
end

function SingleColumnAtmosOceanModel(atmos, ocean)
    grid = ocean.model.grid
    clock = ocean.model.clock
    surface = nothing

    ℂᵀ = ClimaAtmos.Parameters.thermodynamics_params(atmos.integrator.p.params) 
    ℂᴶ = ClimaAtmos.Parameters.surface_fluxes_params(atmos.integrator.p.params) 

    parameters = (
        thermodynamics = ℂᵀ,
        surface_fluxes = ℂᴶ,
    )

    return SingleColumnAtmosOceanModel(clock, grid, atmos, ocean, surface, parameters)
end

"""
    air_sea_turbulent_fluxes(ℋₐ, Tₒ, uₐ, vₐ, uₒ, vₒ, zₐ, zₒ, ℂᵀ, ℂᴶ)

Compute turbluent air-sea fluxes of sensible heat, latent heat, water vapor, and momentum.
"""
@inline function air_sea_turbulent_fluxes(ℂᵀ, ℂᴶ, ℋₐ, uₐ, vₐ, zₐ, Tₒ, uₒ, vₒ, zₒ)

    # Assemble atmosphere dynamical state
    𝕌ₐ = SurfaceFluxes.StateValues(zₐ, SVector(uₐ, vₐ), ℋₐ)

    # Extrapolate to get surface density
    cvₘ = Thermodynamics.cv_m(ℂᵀ, ℋₐ)
    Rₐ = Thermodynamics.gas_constant_air(ℂᵀ, ℋₐ)
    κₐ = cvₘ / Rₐ
    ρₐ = Thermodynamics.air_density(ℂᵀ, ℋₐ)
    Tₐ = Thermodynamics.air_temperature(ℂᵀ, ℋₐ)
    ρₛ = ρₐ * (Tₒ / Tₐ)^κₐ

    # Compute saturation vapor specific humidity over seawater
    FT = typeof(ρₛ)

    # Reference value for the mole fraction of H₂O in seawater.
    # The actual value depends on salinity.
    x_H₂O = convert(FT, 0.98) 
    q★ = x_H₂O * saturation_vapor_specific_humidity(ℂᵀ, Tₒ, ρₛ, Liquid())
    
    qₛ = PhasePartition(q★)
    ℋₛ = PhaseNonEquil_ρTq(ℂᵀ, ρₛ, Tₒ, qₛ)
    𝕌ₛ = SurfaceFluxes.StateValues(zₒ, SVector(uₒ, vₒ), ℋₛ)

    # Roughness and gustiness
    ℓₘ = convert(FT, 1e-4)
    ℓₕ = convert(FT, 1e-4)
    Ug = convert(FT, 1e-1)

    input = ValuesOnly(𝕌ₐ, 𝕌ₛ, ℓₘ, ℓₕ, gustiness=Ug)
    fluxes = surface_conditions(ℂᴶ, input)

    return (sensible_heat = fluxes.shf,
            latent_heat = fluxes.lhf,
            vapor = fluxes.evaporation,
            x_momentum = fluxes.ρτxz,
            y_momentum = fluxes.ρτyz)
end

function time_step!(scm::SingleColumnAtmosOceanModel, Δt)

    atmos = scm.atmos
    ocean = scm.ocean

    # Compute fluxes
    Ûʰ = atmos.integrator.u.c.uₕ # covariant velocity
    Ûʰ₁ = level(Ûʰ, 1)
    ℋ = atmos.integrator.p.precomputed.ᶜts
    ℋ₁ = level(ℋ, 1)

    Ûʰ = atmos.integrator.u.c.uₕ # covariant velocity
    Ûʰ₁ = level(Ûʰ, 1)
    Uʰ₁ = ClimaCore.Geometry.UVVector.(Ûʰ₁)
    u₁ = Uʰ₁.components.data.:1
    v₁ = Uʰ₁.components.data.:2

    X = coordinate_field(atmos.integrator.u.c)
    z₁ = level(X.z, 1)

    # Build ocean state
    surface_layer_space = axes(u₁)
    z₀ = zeros(surface_layer_space)
    u₀ = zeros(surface_layer_space)
    v₀ = zeros(surface_layer_space)
    T₀ = zeros(surface_layer_space)

    Nz = size(ocean.model.grid, 3)
    uₒ = ocean.model.velocities.u[1, 1, Nz]
    vₒ = ocean.model.velocities.v[1, 1, Nz]
    Tₒ = ocean.model.tracers.T[1, 1, Nz]
    Sₒ = ocean.model.tracers.S[1, 1, Nz]

    z₀ .= 0
    u₀ .= uₒ
    v₀ .= vₒ
    T₀ .= Tₒ .+ 273.15

    ℂᵀ = scm.parameters.thermodynamics
    ℂᴶ = scm.parameters.surface_fluxes

    turbulent_fluxes = air_sea_turbulent_fluxes.(Ref(ℂᵀ), Ref(ℂᴶ),
                                                 ℋ₁, u₁, v₁, z₁,
                                                 T₀, u₀, v₀, z₀)

    @show turbulent_fluxes 

    # 1. Set turbulent fluxes in ocean
    # 2. Set turbulent fluxes in atmos
    
    c = atmos.integrator.p.scratch.ᶠtemp_scalar
    𝒢 = Fields.level(Fields.local_geometry_field(c), half)
    ρwh = atmos.integrator.p.precomputed.sfc_conditions.ρ_flux_h_tot 
    ρwq = atmos.integrator.p.precomputed.sfc_conditions.ρ_flux_q_tot 
    ρτ = atmos.integrator.p.precomputed.sfc_conditions.ρ_flux_uₕ

    ρτxz = turbulent_fluxes.x_momentum
    ρτyz = turbulent_fluxes.y_momentum
    F = turbulent_fluxes.vapor
    Qv = turbulent_fluxes.latent_heat
    Qc = turbulent_fluxes.sensible_heat
    ΣQ = Qv + Qc

    @. ρτ = vector_flux(ρτxz, ρτyz, 𝒢)
    @. ρwq = scalar_flux(F, 𝒢)
    @. ρwh = scalar_flux(ΣQ, 𝒢)
    
    # 3. Set radiative fluxes in ocean
    # 4. Update surface state in atmos
    radiation = atmos.integrator.p.radiation.radiation_model
    radiation.surface_temperature .= T₀
    radiation.direct_sw_surface_albedo .= 0.03
    radiation.diffuse_sw_surface_albedo .= 0.03

    # Step forward atmos
    # step!(scm.atmos, Δt, true)
    #
    # Step forward ocean
    # ocean.Δt = Δt
    # time_step!(scm.ocean)

    return nothing
end

