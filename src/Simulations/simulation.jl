using OrderedCollections: OrderedDict

mutable struct Simulation{ML, DT, CL}
    model :: ML
    clock :: CL
    time_step :: DT
    stop_time :: DT
    stop_iteration :: Float64
    wall_time_limit :: Float64
    callbacks :: OrderedDict{Symbol, Callback}
    run_wall_time :: Float64
    running :: Bool
    initialized :: Bool
    verbose :: Bool
end

"""
    Simulation(model; time_step,
               verbose = true,
               stop_iteration = Inf,
               stop_time = Inf,
               wall_time_limit = Inf)

Construct a `Simulation` for a `model` with `time_step`.

Keyword arguments
=================

- `time_step`: Required keyword argument specifying the simulation time step. Can be a `Number`
               for constant time steps or a `TimeStepWizard` for adaptive time-stepping.

- `stop_iteration`: Stop the simulation after this many iterations.

- `stop_time`: Stop the simulation once this much model clock time has passed.

- `wall_time_limit`: Stop the simulation if it's been running for longer than this many
                     seconds of wall clock time.
"""
function Simulation(model; time_step,
                    verbose = true,
                    stop_iteration = Inf,
                    stop_time = Inf,
                    wall_time_limit = Inf)

   if stop_iteration == Inf && stop_time == Inf && wall_time_limit == Inf
       @warn "This simulation will run forever as stop iteration = stop time " *
             "= wall time limit = Inf."
   end

   callbacks = OrderedDict{Symbol, Callback}()

   callbacks[:stop_time_exceeded] = Callback(stop_time_exceeded)
   callbacks[:stop_iteration_exceeded] = Callback(stop_iteration_exceeded)
   callbacks[:wall_time_limit_exceeded] = Callback(wall_time_limit_exceeded)

   return Simulation(model,
                     time_step,
                     stop_time,
                     Float64(stop_iteration),
                     Float64(wall_time_limit),
                     callbacks,
                     0.0,
                     false,
                     false,
                     verbose)
end

function Base.show(io::IO, s::Simulation)
    modelstr = summary(s.model)
    return print(io, "Simulation of ", modelstr, "\n",
                     "├── Next time step: $(prettytime(s.time_step))", "\n",
                     "├── Elapsed wall time: $(prettytime(s.run_wall_time))", "\n",
                     "├── Wall time per iteration: $(prettytime(s.run_wall_time / iteration(s)))", "\n",
                     "├── Stop time: $(prettytime(s.stop_time))", "\n",
                     "├── Stop iteration : $(s.stop_iteration)", "\n",
                     "├── Wall time limit: $(s.wall_time_limit)", "\n",
                     "└── Callbacks: $(ordered_dict_show(s.callbacks, "│"))")
end

#####
##### Utilities
#####

"""
    iteration(sim::Simulation)

Return the current simulation iteration.
"""
function iteration end

#####
##### Default stop criteria callback functions
#####

wall_time_msg(sim) = string("Simulation is stopping after running for ", run_wall_time(sim), ".")

function stop_iteration_exceeded(sim)
    if sim.model.clock.iteration >= sim.stop_iteration
        if sim.verbose
            msg = string("Model iteration ", iteration(sim), " equals or exceeds stop iteration ", Int(sim.stop_iteration), ".")
            @info wall_time_msg(sim) 
            @info msg
        end

        sim.running = false 
    end

    return nothing
end

function stop_time_exceeded(sim)
    if sim.model.clock.time >= sim.stop_time
        if sim.verbose
            msg = string("Simulation time ", prettytime(sim), " equals or exceeds stop time ", prettytime(sim.stop_time), ".")
            @info wall_time_msg(sim) 
            @info msg
        end

        sim.running = false 
    end

    return nothing
end

function wall_time_limit_exceeded(sim)
    if sim.run_wall_time >= sim.wall_time_limit
        if sim.verbose
            msg = string("Simulation run time ", run_wall_time(sim), " equals or exceeds wall time limit ", prettytime(sim.wall_time_limit), ".")
            @info wall_time_msg(sim)
            @info msg
        end

        sim.running = false 
    end

    return nothing
end

