using ClimaEarth
using ClimaAtmos
using ClimaOcean

# Define parameters common to all components
# ...


# Set up ocean component
# ocean specific parameters...
ocean = global_ocean_simulation(ocean_parameters...)

# Set up atmos component
# atmos specific parameters...
atmos = global_atmos_simulation(atmos_parameters...)

# Set up land component
# land specific parameters...
land = global_ocean_simulation(land_parameters...)

earth_system_model = EarthSystemModel(atmos, ocean, land; coupled_parameters...)
# Set the initial condition
