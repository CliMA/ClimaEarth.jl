# This document lays out the functions that must be extended to
# use an atmospheric simulation in ClimaOcean.

using ClimaOcean

import Oceananigans: time_step!
import Oceananigans.Models: update_model_field_time_series!

# Make sure the atmospheric parameters from SpeedyWeather can be used in the compute fluxes function
import ClimaOcean.OceanSeaIceModels.Atmospheres: thermodynamics_parameters, 
                                                 boundary_layer_height, 
                                                 surface_layer_height,
                                                 regrid_fluxes_to_atmospheric_model!, 
                                                 interpolate_atmospheric_state!

# Assuming there is a type called ClimaAtmosSimulation encapsulating
# the atmosphere simulation:

# This can be left blank
update_model_field_time_series!(::ClimaAtmosSimulation, time) = nothing

# Take one time-step
time_step!(atmos::ClimaAtmosSimulation) = nothing

# The height of surface variables, used by the turbulent flux solver
surface_layer_height(s::ClimaAtmosSimulation) = 10 # meters, for example

# This is a parameter that is used in the computation of the fluxes,
# It probably should not be here but in the similarity theory type.
boundary_layer_height(atmos::ClimaAtmosSimulation) = 600

# Note: possibly, can use the atmos thermodynamic parameters directly here.
thermodynamics_parameters(atmos::ClimaAtmosSimulation) = 
    PrescribedAtmosphereThermodynamicsParameters(Float32)

"""
    interpolate_atmospheric_state!(surface_atmosphere_state, 
                                        interpolated_prescribed_freshwater_flux, 
                                        atmos::ClimaAtmosSimulation, 
                                        grid, clock)

Interpolate the atmospheric state in `atmos` to `surface_atmospheric_state`, a
the collection of `Field`s needed to compute turbulent fluxes.
"""
function interpolate_atmospheric_state!(surface_atmosphere_state, 
                                        interpolated_prescribed_freshwater_flux, 
                                        atmos::ClimaAtmosSimulation, 
                                        grid, clock)

    #=
    # Get the atmospheric state on the ocean grid
    #
    # These are required if the interpolation must be done on the CPU
    # Otherwise, these are unnecessary.
    ua = on_architecture(CPU(), surface_atmosphere_state.u)
    va = on_architecture(CPU(), surface_atmosphere_state.v)
    Ta = on_architecture(CPU(), surface_atmosphere_state.T)
    qa = on_architecture(CPU(), surface_atmosphere_state.q)
    pa = on_architecture(CPU(), surface_atmosphere_state.p)
    Qs = on_architecture(CPU(), surface_atmosphere_state.Qs)
    Qℓ = on_architecture(CPU(), surface_atmosphere_state.Qℓ)

    λ,  φ,  _ = Oceananigans.Grids.nodes(grid, Center(), Center(), Center(), with_halos=true) 

    λ = Array(vec(on_architecture(CPU(), λ)))
    φ = Array(vec(on_architecture(CPU(), φ)))

    spectral_grid = atmos.model.spectral_grid
    interpolator = RingGrids.AnvilInterpolator(Float32, spectral_grid.Grid, spectral_grid.nlat_half, length(λ))
    RingGrids.update_locator!(interpolator, λ, φ)

    RingGrids.interpolate!(vec(view(ua, :, :, 1)), atmos.diagnostic_variables.grid.u_grid[:, end],            interpolator)
    RingGrids.interpolate!(vec(view(va, :, :, 1)), atmos.diagnostic_variables.grid.v_grid[:, end],            interpolator)
    RingGrids.interpolate!(vec(view(Ta, :, :, 1)), atmos.diagnostic_variables.grid.temp_grid[:, end],         interpolator)
    RingGrids.interpolate!(vec(view(qa, :, :, 1)), atmos.diagnostic_variables.grid.humid_grid[:, end],        interpolator)
    RingGrids.interpolate!(vec(view(pa, :, :, 1)), exp.(atmos.diagnostic_variables.grid.pres_grid[:, end]),   interpolator)
    RingGrids.interpolate!(vec(view(Qs, :, :, 1)), atmos.diagnostic_variables.physics.surface_shortwave_down, interpolator)
    RingGrids.interpolate!(vec(view(Qℓ, :, :, 1)), atmos.diagnostic_variables.physics.surface_longwave_down,  interpolator)

    surface_atmosphere_state.u  .= ua 
    surface_atmosphere_state.v  .= va 
    surface_atmosphere_state.T  .= Ta 
    surface_atmosphere_state.q  .= qa 
    surface_atmosphere_state.p  .= pa 
    surface_atmosphere_state.Qs .= Qs
    surface_atmosphere_state.Qℓ .= Qℓ
    =#

    return nothing
end

"""
    regrid_fluxes_to_atmospheric_model!(atmos::ClimaAtmosSimulation, turbulent_fluxes)

Regrid `turbulent_fluxes` computed by ClimaOcean and pass into the `atmos` simulation.
"""
function regrid_fluxes_to_atmospheric_model!(atmos::ClimaAtmosSimulation, turbulent_fluxes)
    
    Qs = turbulent_fluxes.sensible_heat
    Mv = turbulent_fluxes.water_vapor

    # Example computation for SpeedyWeather, using prebuilt interpolator object:
    # regrid_flux!(atmos.diagnostic_variables.physics.sensible_heat_flux, Qs, interpolator, weights)
    # regrid_flux!(atmos.diagnostic_variables.physics.evaporative_flux,   Mv, interpolator, weights)

    return nothing
end

#=
# Utility function for the above example:
function regrid_flux!(speedy_flux, climaocean_flux, interpolator, weights)
    interpolate!(tmp_field, climaocean_flux, weights) 
    tmp_field = on_architecture(CPU(), interior(tmp_field, :, :, 1))
    tmp_field[isnan.(tmp_field)] .= 0 # We do not have antarctic
    flux = FullClenshawGrid(vec(reverse(tmp_field, dims=2)))
    RingGrids.interpolate!(speedy_flux, flux, interpolator)
    return nothing
end
=#
