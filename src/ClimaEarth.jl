module ClimaEarth

"""
    initialize!(obj)

Called at the beginning of a simulation to initialize `obj`,
which may be `Simulation`, `Callback`, `Callback.schedule`, or any
number of other `obj`ects.
"""
function initialize! end

include("Utils.jl")
include("Schedules.jl")
include("Simulations/Simulations.jl")

using .Simulations: Simulation

end # module ClimaEarth

