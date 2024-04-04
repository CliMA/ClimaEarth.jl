using ClimaCore
using ClimaCore.Fields: level, coordinate_field
using ClimaCore.Geometry: UVVector
using ClimaCore.Geometry: Covariant12Vector

include("earth_system_model.jl")
include("atmos_simulation.jl")
include("ocean_simulation.jl")

atmos = AtmosSimulation()
ocean = ocean_single_column_simulation()

# TODO: figure out how to set initial conditions
# center_space = axes(atmos.integrator.u.c)
# set!(atmos, rho=1.2, T=DryAdiabaticProfile())
# u = 1
Ûʰ = atmos.integrator.u.c.uₕ # covariant velocity
Uʰ = UVVector.(Ûʰ) # contravariant velocity
Uʰ.components.data.:1 .= 10 # set u = 10
Uʰ.components.data.:2 .= 0 # set v = 0
Ûʰ .= Covariant12Vector.(Uʰ)

scm = SingleColumnAtmosOceanModel(atmos, ocean)

time_step!(scm, 60)

