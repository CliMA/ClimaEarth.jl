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

    λ = λnodes(grid, Center(), Center(), Center(), with_halos=true)
    φ = φnodes(grid, Center(), Center(), Center(), with_halos=true)

    if grid isa LatitudeLongitudeGrid
        # Note that we have to include halos here, or fill_halo_regions! later.
        Nx, Ny, Nz = size(grid)
        λ = λ[0:Nx+1]
        φ = φ[0:Ny+1]
        λ = reshape(λ, Nx+2, 1)
        φ = reshape(φ, 1, Ny+2)
    end

    # Relevant code for moving interpolation into a kernel:
    #
    # https://github.com/CliMA/ClimaCore.jl/ \
    #   blob/7624132c30836365a1be76865d19c7722070b98d/ext/cuda/remapping_interpolate_array.jl#L121
    #
    Xh = @. CC.Geometry.LatLongPoint(φ, λ)

    uah = CC.Geometry.UVVector.(CC.Spaces.level(atmosphere.integrator.u.c.uₕ, 1))
    ui = CC.Remapping.interpolate(uah.components.data.:1, Xh, nothing)
    vi = CC.Remapping.interpolate(uah.components.data.:2, Xh, nothing)

    # T, q, p from thermodynamic state ts?
    tsa = CC.Spaces.level(atmosphere.integrator.p.precomputed.ᶜts, 1)
    ℂa = atmosphere.integrator.p.params.thermodynamics_params
    Ta = Thermodynamics.air_temperature.(ℂa, tsa)
    pa = Thermodynamics.air_pressure.(ℂa, tsa)
    qa = Thermodynamics.total_specific_humidity.(ℂa, tsa)

    Ti = CC.Remapping.interpolate(Ta, Xh, nothing)
    pi = CC.Remapping.interpolate(pa, Xh, nothing)
    qi = CC.Remapping.interpolate(qa, Xh, nothing)

    ui = OffsetArray(ui, -1, -1)
    vi = OffsetArray(vi, -1, -1)
    Ti = OffsetArray(Ti, -1, -1)
    qi = OffsetArray(qi, -1, -1)
    pi = OffsetArray(pi, -1, -1)

    #=
    ui = on_architecture(arch, ui)
    vi = on_architecture(arch, vi)
    Ti = on_architecture(arch, Ti)
    qi = on_architecture(arch, qi)
    pi = on_architecture(arch, pi)
    
    kp = CO.OceanSeaIceModels.InterfaceComputations.interface_kernel_parameters(grid)
    arch = OC.architecture(grid)
    atmos_state = (u=ui, v=vi, T=Ti, q=qi, p=pi)
    OC.Utils.launch!(arch, grid, kp,
                     _interpolate_atmosphere_state!,
                     exchange_atmosphere_state,
                     atmos_state)
    =#

    ue  = view(exchange_atmosphere_state.u.data,  0:Nx+1, 0:Ny+1, 1)
    ve  = view(exchange_atmosphere_state.v.data,  0:Nx+1, 0:Ny+1, 1)
    Te  = view(exchange_atmosphere_state.T.data,  0:Nx+1, 0:Ny+1, 1)
    qe  = view(exchange_atmosphere_state.q.data,  0:Nx+1, 0:Ny+1, 1)
    pe  = view(exchange_atmosphere_state.p.data,  0:Nx+1, 0:Ny+1, 1)
    Qse = view(exchange_atmosphere_state.Qs.data, 0:Nx+1, 0:Ny+1, 1)
    Qℓe = view(exchange_atmosphere_state.Qℓ.data, 0:Nx+1, 0:Ny+1, 1)

    copyto!(ue, ui)
    copyto!(ve, vi)
    copyto!(Te, Ti)
    copyto!(qe, qi)
    copyto!(pe, pi)

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

