module ClimaEarth

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

    exchange_atmosphere_state = interfaces.exchanger.exchange_atmosphere_state
    grid = exchange_atmosphere_state.u.grid

    λ = λnodes(grid, Center(), Center(), Center())
    φ = φnodes(grid, Center(), Center(), Center())

    Nx, Ny, Nz = size(grid)

    if grid isa LatitudeLongitudeGrid
        λ = reshape(λ, Nx, 1)
        φ = reshape(φ, 1, Ny)
    end

    hcoords = @. CC.Geometry.LatLongPoint(φ, λ)

    uah = CC.Geometry.UVVector.(CC.Spaces.level(atmosphere.integrator.u.c.uₕ, 1))
    ui = CC.Remapping.interpolate(uah.components.data.:1, hcoords, nothing)
    vi = CC.Remapping.interpolate(uah.components.data.:2, hcoords, nothing)

    # T, q, p from thermodynamic state ts?
    tsa = CC.Spaces.level(atmosphere.integrator.p.precomputed.ᶜts, 1)
    ℂa = atmosphere.integrator.p.params.thermodynamics_params
    Ta = Thermodynamics.air_temperature.(ℂa, tsa)
    pa = Thermodynamics.air_pressure.(ℂa, tsa)
    qa = Thermodynamics.total_specific_humidity.(ℂa, tsa)

    Ti = CC.Remapping.interpolate(Ta, hcoords, nothing)
    pi = CC.Remapping.interpolate(pa, hcoords, nothing)
    qi = CC.Remapping.interpolate(qa, hcoords, nothing)

    # interior(exchange_atmosphere_state.u, :, :, 1) ?
    # or write a kernel and set them all.
    exchange_atmosphere_state.u  .= ui
    exchange_atmosphere_state.v  .= vi
    exchange_atmosphere_state.T  .= Ti
    exchange_atmosphere_state.q  .= qi
    exchange_atmosphere_state.p  .= pi
    exchange_atmosphere_state.Qs .= 0
    exchange_atmosphere_state.Qℓ .= 0

    return nothing
end

atmosphere_exchanger(atmosphere::CA.AtmosSimulation, exchange_grid) = nothing

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
    @. ρwh = CA.SurfaceConditions.vector_from_component(Qv, 𝒢) + CA.SurfaceConditions.vector_from_component(Qc, 𝒢)

    # Mass or volume flux: check units
    # ρwq = atmos.integrator.p.precomputed.sfc_conditions.ρ_flux_q_tot
    # @. ρwq = CA.SurfaceConditions.vector_from_component(Fv, 𝒢)
    
    # TODO: validate this?
    ρτ = atmosphere.integrator.p.precomputed.sfc_conditions.ρ_flux_uₕ  
    @. ρτ = CA.SurfaceConditions.tensor_from_components(ρτx, ρτy, 𝒢)

    return nothing
end

end # module ClimaEarth

