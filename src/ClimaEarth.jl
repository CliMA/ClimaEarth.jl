module ClimaEarth

using OffsetArrays
using KernelAbstractions

import ClimaAtmos as CA
import ClimaOcean as CO
import Oceananigans as OC

# This document lays out the functions that must be extended to
# use an atmospheric simulation in ClimaOcean.

import Oceananigans.TimeSteppers: time_step!
import Oceananigans.Models: update_model_field_time_series!

# Make sure the atmospheric parameters from SpeedyWeather can be used in the compute fluxes function
import ClimaOcean.OceanSeaIceModels.PrescribedAtmospheres: 
    thermodynamics_parameters, 
    boundary_layer_height, 
    surface_layer_height

import ClimaOcean.OceanSeaIceModels: compute_net_atmosphere_fluxes!

import ClimaOcean.OceanSeaIceModels.InterfaceComputations:
    atmosphere_exchanger,
    initialize!,
    StateExchanger,
    interpolate_atmosphere_state!

using ClimaOcean.OceanSeaIceModels: OceanSeaIceModel

const ClimaCoupledModel = OceanSeaIceModel{<:Any, <:CA.AtmosSimulation}

# This can be left blank:
update_model_field_time_series!(::CA.AtmosSimulation, time) = nothing

# Take one time-step
function time_step!(atmos::CA.AtmosSimulation, Δt)
    # TODO: check if the time-step can be changed.
    @assert Δt == atmos.integrator.dt
    CA.SciMLBase.step!(atmos.integrator)
    return nothing
end

# The height of near-surface variables used in the turbulent flux solver
surface_layer_height(s::CA.AtmosSimulation) = 10 # meters, for example

# This is a parameter that is used in the computation of the fluxes,
# It probably should not be here but in the similarity theory type.
boundary_layer_height(atmos::CA.AtmosSimulation) = 600

# Note: possibly, can use the atmos thermodynamic parameters directly here.
thermodynamics_parameters(atmos::CA.AtmosSimulation) = 
    atmos.integrator.p.params.thermodynamics_params

Base.summary(::CA.AtmosSimulation) = "ClimaAtmos.AtmosSimulation"

import ClimaCore as CC
using Oceananigans.Grids: λnodes, φnodes, LatitudeLongitudeGrid
using Oceananigans.Fields: Center
using Thermodynamics
Thermodynamics.print_warning() = false
#Thermodynamics.print_warning() = true

"""
    interpolate_atmospheric_state!(surface_atmosphere_state, 
                                        interpolated_prescribed_freshwater_flux, 
                                        atmos::ClimaAtmosSimulation, 
                                        grid, clock)

Interpolate the atmospheric state in `atmos` to `surface_atmospheric_state`, a
the collection of `Field`s needed to compute turbulent fluxes.
"""
function interpolate_atmosphere_state!(interfaces, 
                                       atmosphere::CA.AtmosSimulation, 
                                       coupled_model)

    remapper = interfaces.exchanger.atmosphere_exchanger
    exchange_atmosphere_state = interfaces.exchanger.exchange_atmosphere_state

    ue  = parent(exchange_atmosphere_state.u)
    ve  = parent(exchange_atmosphere_state.v)
    Te  = parent(exchange_atmosphere_state.T)
    qe  = parent(exchange_atmosphere_state.q)
    pe  = parent(exchange_atmosphere_state.p)

    ue = dropdims(ue, dims=3)
    ve = dropdims(ve, dims=3)
    Te = dropdims(Te, dims=3)
    qe = dropdims(qe, dims=3)
    pe = dropdims(pe, dims=3)

    Uah = CC.Geometry.UVVector.(CC.Spaces.level(atmosphere.integrator.u.c.uₕ, 1))
    ua = Uah.components.data.:1
    va = Uah.components.data.:2

    tsa = CC.Spaces.level(atmosphere.integrator.p.precomputed.ᶜts, 1)
    ℂa = atmosphere.integrator.p.params.thermodynamics_params
    Ta = Thermodynamics.air_temperature.(ℂa, tsa)
    pa = Thermodynamics.air_pressure.(ℂa, tsa)
    qa = Thermodynamics.total_specific_humidity.(ℂa, tsa)

    CC.Remapping.interpolate!(ue, remapper, ua)
    CC.Remapping.interpolate!(ve, remapper, va)
    CC.Remapping.interpolate!(Te, remapper, Ta)
    CC.Remapping.interpolate!(pe, remapper, pa)
    CC.Remapping.interpolate!(qe, remapper, qa)

    #=
    # This is needed, unless the above computations include the halos.
    # OC.fill_halo_regions!(exchange_atmosphere_state.u)
    # OC.fill_halo_regions!(exchange_atmosphere_state.v)
    # OC.fill_halo_regions!(exchange_atmosphere_state.T)
    # OC.fill_halo_regions!(exchange_atmosphere_state.q)
    # OC.fill_halo_regions!(exchange_atmosphere_state.p)
    =#

    return nothing
end

# Note: this just copies, for now.
@kernel function _interpolate_atmosphere_state!(exchange_state, atmos_state)
    i, j = @index(Global, NTuple)
    @inbounds begin
        exchange_state.u[i, j, 1] = atmos_state.u[i, j]
        exchange_state.v[i, j, 1] = atmos_state.v[i, j]
        exchange_state.T[i, j, 1] = atmos_state.T[i, j]
        exchange_state.q[i, j, 1] = atmos_state.q[i, j]
        exchange_state.p[i, j, 1] = atmos_state.p[i, j]
    end
end

function atmosphere_exchanger(atmosphere::CA.AtmosSimulation, exchange_grid)
    λ = λnodes(exchange_grid, Center(), Center(), Center(), with_halos=true)
    φ = φnodes(exchange_grid, Center(), Center(), Center(), with_halos=true)

    if exchange_grid isa LatitudeLongitudeGrid
        λ = reshape(λ, length(λ), 1)
        φ = reshape(φ, 1, length(φ))
    end

    Xh = @. CC.Geometry.LatLongPoint(φ, λ)
    space = axes(atmosphere.integrator.u.c)
    surface = CC.Spaces.level(space, 1)

    # Note: buffer_length gives the maximum number of variables that can be remapped
    # within a single kernel.
    #atmos_to_exchange_remapper = CC.Remapping.Remapper(space, Xh, nothing, buffer_length=5)
    atmos_to_exchange_remapper = CC.Remapping.Remapper(surface, Xh, nothing, buffer_length=5)

    return atmos_to_exchange_remapper
end

initialize!(::StateExchanger, ::CA.AtmosSimulation) = nothing

#=
struct AtmosOceanExchanger{A2E, E2A}
    atmos_to_exchange_regridder :: A2E
    exchange_to_atmos_regridder :: E2A
end

function atmosphere_exchanger(atmosphere::CA.AtmosSimulation, exchange_grid)
    space3 = axes(atmosphere.integrator.p.precomputed.sfc_conditions.ts)
    space2 = CC.Spaces.SpectralElementSpace2D(space3.grid.full_grid.horizontal_grid)
    regridder = ClimaUtilities.Regridders.InterpolationsRegridder(space2)

end
=#

using ClimaUtilities
using ClimaCore.Utilities: half

function compute_net_atmopshere_fluxes!(coupled_model::ClimaCoupledModel)
    atmosphere = coupled_model.atmosphere
    ocean = coupled_model.ocean
    ocean_grid = ocean.model.grid
    interfaces = coupled_model.interfaces

    # Regridding back to atmos
    space3 = axes(atmosphere.integrator.p.precomputed.sfc_conditions.ts)
    space2 = CC.Spaces.SpectralElementSpace2D(space3.grid.full_grid.horizontal_grid)
    regridder = ClimaUtilities.Regridders.InterpolationsRegridder(space2)

    if ocean_grid isa LatitudeLongitudeGrid
        λo = λnodes(ocean_grid, Center(), Center(), Center())
        φo = φnodes(ocean_grid, Center(), Center(), Center())

        Qv_ao = interfaces.atmosphere_ocean_interface.fluxes.latent_heat
        Qc_ao = interfaces.atmosphere_ocean_interface.fluxes.sensible_heat
        Fv_ao = interfaces.atmosphere_ocean_interface.fluxes.water_vapor
        ρτx_ao = interfaces.atmosphere_ocean_interface.fluxes.x_momentum
        ρτy_ao = interfaces.atmosphere_ocean_interface.fluxes.y_momentum
    elseif ocean_grid isa OrthogonalSphericalShellGrid
        # One quick and dirty option: https://github.com/CliMA/OrthogonalSphericalShellGrids.jl/pull/29
        error("Not supported yet!")
    end

    # Regrid to ClimaCore grid
    Qv  = ClimaUtilities.Regridders.regrid(regridder, interior(Qv_ao, :, :, 1),  (λi, φi))
    Qc  = ClimaUtilities.Regridders.regrid(regridder, interior(Qc_ao, :, :, 1),  (λi, φi))
    Fv  = ClimaUtilities.Regridders.regrid(regridder, interior(Fv_ao, :, :, 1),  (λi, φi))
    ρτx = ClimaUtilities.Regridders.regrid(regridder, interior(ρτx_ao, :, :, 1), (λi, φi))
    ρτy = ClimaUtilities.Regridders.regrid(regridder, interior(ρτy_ao, :, :, 1), (λi, φi))

    # Project onto a vector!
    # :eyes https://github.com/CliMA/ClimaEarth.jl/pull/5/files
    c = atmos.integrator.p.scratch.ᶠtemp_scalar
    𝒢 = ClimaCore.Fields.level(ClimaCore.Fields.local_geometry_field(c), half)
    ρwh = atmos.integrator.p.precomputed.sfc_conditions.ρ_flux_h_tot
    @. ρwh = CA.SurfaceConditions.vector_from_component(Qv, 𝒢) +
             CA.SurfaceConditions.vector_from_component(Qc, 𝒢)

    # Mass or volume flux: check units
    ρwq = atmos.integrator.p.precomputed.sfc_conditions.ρ_flux_q_tot
    @. ρwq = CA.SurfaceConditions.vector_from_component(Fv, 𝒢)
    
    # TODO: validate this?
    ρτ = atmosphere.integrator.p.precomputed.sfc_conditions.ρ_flux_uₕ  
    @. ρτ = CA.SurfaceConditions.tensor_from_components(ρτx, ρτy, 𝒢)

    return nothing
end

end # module ClimaEarth

