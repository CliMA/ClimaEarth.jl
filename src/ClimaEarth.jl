module ClimaEarth

using OffsetArrays
using KernelAbstractions
using Statistics

import ClimaAtmos as CA
import ClimaOcean as CO
import Oceananigans as OC

import CUDA

# This document lays out the functions that must be extended to
# use an atmospheric simulation in ClimaOcean.

import Oceananigans.TimeSteppers: time_step!
import Oceananigans.Models: update_model_field_time_series!

# Make sure the atmospheric parameters from SpeedyWeather can be used in the compute fluxes function
import ClimaOcean.OceanSeaIceModels.PrescribedAtmospheres: 
    thermodynamics_parameters, 
    boundary_layer_height, 
    surface_layer_height

import ClimaOcean.OceanSeaIceModels:
    compute_net_atmosphere_fluxes!

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

    interpolator = interfaces.exchanger.atmosphere_exchanger.to_exchange_interp
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

    # TODO: can we avoid allocating for Ta, pa, qa?
    tsa = CC.Spaces.level(atmosphere.integrator.p.precomputed.ᶜts, 1)
    ℂa = atmosphere.integrator.p.params.thermodynamics_params
    Ta = Thermodynamics.air_temperature.(ℂa, tsa)
    pa = Thermodynamics.air_pressure.(ℂa, tsa)
    qa = Thermodynamics.total_specific_humidity.(ℂa, tsa)

    #=
    # TODO: make this work without allocation
    #       make sure that Remapper(args...; buffer_length=5)
    #       or whatever it needs to be
    exchange_fields = cat(ue, ve, Te, pe, qe, dims=3)
    atmos_fields    = [ua, va, Ta, pa, qa]
    CC.Remapping.interpolate!(exchange_fields, remapper, atmos_fields)
    =#

    CC.Remapping.interpolate!(ue, interpolator, ua)
    CC.Remapping.interpolate!(ve, interpolator, va)
    CC.Remapping.interpolate!(Te, interpolator, Ta)
    CC.Remapping.interpolate!(pe, interpolator, pa)
    CC.Remapping.interpolate!(qe, interpolator, qa)

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

#=
mutable struct AtmosphereExchanger
    atmosphere_to_exchange
    exchange_to_atmos
end
=#

function atmosphere_exchanger(atmosphere::CA.AtmosSimulation, exchange_grid)
    λ = λnodes(exchange_grid, Center(), Center(), Center(), with_halos=true)
    φ = φnodes(exchange_grid, Center(), Center(), Center(), with_halos=true)

    if exchange_grid isa LatitudeLongitudeGrid
        λ = reshape(λ, length(λ), 1)
        φ = reshape(φ, 1, length(φ))
    end

    Xh = @. CC.Geometry.LatLongPoint(φ, λ)
    space = axes(atmosphere.integrator.u.c)
    first_level = CC.Spaces.level(space, 1)

    # Note: buffer_length gives the maximum number of variables that can be remapped
    # within a single kernel.
    to_exchange_interp = CC.Remapping.Remapper(first_level, Xh, nothing, buffer_length=1)

    # Make a remapper for exchange_to_atmos regridding
    space3 = axes(atmosphere.integrator.p.precomputed.sfc_conditions.ts)
    space2 = CC.Spaces.SpectralElementSpace2D(space3.grid.full_grid.horizontal_grid)
    regridder = ClimaUtilities.Regridders.InterpolationsRegridder(space2)
    atmos_surface_points = regridder.coordinates

    if exchange_grid isa OC.Grids.OrthogonalSphericalShellGrid
        # One quick and dirty option: https://github.com/CliMA/OrthogonalSphericalShellGrids.jl/pull/29
        error("Not supported yet!")
    end

    dummy_flux = OC.Field{OC.Center, OC.Center, Nothing}(exchange_grid)
    Qc_a = map_interpolate(atmos_surface_points, dummy_flux)
    Qv_a = map_interpolate(atmos_surface_points, dummy_flux)
    Fv_a = map_interpolate(atmos_surface_points, dummy_flux)
    ρτx_a = map_interpolate(atmos_surface_points, dummy_flux)
    ρτy_a = map_interpolate(atmos_surface_points, dummy_flux)
    turbulent_atmosphere_surface_fluxes = (; Qc_a, Qv_a, Fv_a, ρτx_a, ρτy_a)

    return (; to_exchange_interp, turbulent_atmosphere_surface_fluxes, atmos_surface_points)
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

@inline to_node(pt::CA.ClimaCore.Geometry.LatLongPoint) = pt.long, pt.lat

instantiate(L) = L()

function map_interpolate(points, oc_field::OC.Field) #, loc, grid)
    loc = map(instantiate, OC.Fields.location(oc_field))
    grid = oc_field.grid
    data = oc_field.data

    map(points) do pt
        FT = eltype(pt)
        fᵢ = OC.Fields.interpolate(to_node(pt), data, loc, grid)
        convert(FT, fᵢ)
    end
end

function map_interpolate!(cc_field, points, oc_field::OC.Field)
    loc = map(instantiate, OC.Fields.location(oc_field))
    grid = oc_field.grid
    data = oc_field.data

    map!(cc_field, points) do pt
        FT = eltype(pt)
        fᵢ = OC.Fields.interpolate(to_node(pt), data, loc, grid)
        convert(FT, fᵢ)
    end

    return nothing
end

function compute_net_atmosphere_fluxes!(coupled_model::ClimaCoupledModel)
    atmosphere = coupled_model.atmosphere
    ocean = coupled_model.ocean
    ocean_grid = ocean.model.grid
    interfaces = coupled_model.interfaces
    exchanger = interfaces.exchanger

    atmos_surface_points = exchanger.atmosphere_exchanger.atmos_surface_points
    (; Qc_a, Qv_a, Fv_a, ρτx_a, ρτy_a) = exchanger.atmosphere_exchanger.turbulent_atmosphere_surface_fluxes 

    Qv_e = interfaces.atmosphere_ocean_interface.fluxes.latent_heat
    Qc_e = interfaces.atmosphere_ocean_interface.fluxes.sensible_heat
    Fv_e = interfaces.atmosphere_ocean_interface.fluxes.water_vapor
    ρτx_e = interfaces.atmosphere_ocean_interface.fluxes.x_momentum
    ρτy_e = interfaces.atmosphere_ocean_interface.fluxes.y_momentum

    map_interpolate!(Qc_a,  atmos_surface_points, Qc_e)
    map_interpolate!(Qv_a,  atmos_surface_points, Qv_e)
    map_interpolate!(Fv_a,  atmos_surface_points, Fv_e)
    map_interpolate!(ρτx_a, atmos_surface_points, ρτx_e)
    map_interpolate!(ρτy_a, atmos_surface_points, ρτy_e)

    # Project onto a vector...
    # :eyes https://github.com/CliMA/ClimaEarth.jl/pull/5/files
    c = atmosphere.integrator.p.scratch.ᶠtemp_scalar
    𝒢 = CC.Fields.level(CC.Fields.local_geometry_field(c), half)
    ρwh = atmosphere.integrator.p.precomputed.sfc_conditions.ρ_flux_h_tot
    @. ρwh = CA.SurfaceConditions.vector_from_component(Qv_a, 𝒢) +
             CA.SurfaceConditions.vector_from_component(Qc_a, 𝒢)

    # Mass or volume flux: check units
    ρwq = atmosphere.integrator.p.precomputed.sfc_conditions.ρ_flux_q_tot
    @. ρwq = CA.SurfaceConditions.vector_from_component(Fv_a, 𝒢)
    
    # TODO: validate this?
    ρτ = atmosphere.integrator.p.precomputed.sfc_conditions.ρ_flux_uₕ  
    @. ρτ = CA.SurfaceConditions.tensor_from_components(ρτx_a, ρτy_a, 𝒢)

    return nothing
end

end # module ClimaEarth

