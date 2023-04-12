using Glob

#=
using Oceananigans.Fields: set!
using Oceananigans.OutputWriters: WindowedTimeAverage, checkpoint_superprefix
using Oceananigans.TimeSteppers: unit_time
using Oceananigans: AbstractModel, run_diagnostic!, write_output!

import Oceananigans: initialize!
import Oceananigans.OutputWriters: checkpoint_path, set!
import Oceananigans.TimeSteppers: time_step!
=#

using ..Utils: prettytime

import ClimaEarth: initialize!

# Simulations are for running

#####
##### Time-step "alignment" with output and callbacks scheduled on TimeInterval
#####

function collect_scheduled_activities(sim)
    writers = values(sim.output_writers)
    callbacks = values(sim.callbacks)
    return tuple(writers..., callbacks...)
end

function compute_schedule_aligned_time_step(sim, Δt)
    clock = sim.clock
    activities = collect_scheduled_activities(sim)

    for activity in activities
        aligned_Δt = aligned_time_step(activity.schedule, clock, Δt)
    end

    return aligned_Δt
end

"""
    aligned_time_step(sim, Δt)

Return a time step 'aligned' with `sim.stop_time`, output writer schedules, 
and callback schedules. Alignment with `sim.stop_time` takes precedence.
"""
function aligned_time_step(sim::Simulation, Δt)
    clock = sim.clock
    aligned_Δt = Δt

    # Align time step with output writing and callback execution
    aligned_Δt = compute_schedule_aligned_time_step(sim, aligned_Δt)
    
    # Align time step with simulation stop time
    aligned_Δt = min(aligned_Δt, unit_time(sim.stop_time - clock.time))

    # Temporary fix for https://github.com/CliMA/Oceananigans.jl/issues/1280
    aligned_Δt = aligned_Δt <= 0 ? Δt : aligned_Δt

    return aligned_Δt
end

"""
    run!(simulation; pickup=false)

Run a `simulation` until one of `simulation.stop_criteria` evaluates `true`.
The simulation will then stop.

# Picking simulations up from a checkpoint

Simulations are "picked up" from a checkpoint if `pickup` is either `true`, a `String`, or an
`Integer` greater than 0.

Picking up a simulation sets field and tendency data to the specified checkpoint,
leaving all other model properties unchanged.

Possible values for `pickup` are:

  * `pickup=true` picks a simulation up from the latest checkpoint associated with
    the `Checkpointer` in `simulation.output_writers`.

  * `pickup=iteration::Int` picks a simulation up from the checkpointed file associated
     with `iteration` and the `Checkpointer` in `simulation.output_writers`.

  * `pickup=filepath::String` picks a simulation up from checkpointer data in `filepath`.

Note that `pickup=true` and `pickup=iteration` fails if `simulation.output_writers` contains
more than one checkpointer.
"""
function run!(sim; pickup=false)

    #=
    if we_want_to_pickup(pickup)
        checkpoint_file_path = checkpoint_path(pickup, sim.output_writers)
        set!(sim.model, checkpoint_file_path)
    end
    =#

    sim.initialized = false
    sim.running = true
    sim.run_wall_time = 0.0

    while sim.running
        time_step!(sim)
    end

    return nothing
end

struct TendencyCallsite end
struct UpdateStateCallsite end

const ModelCallsite = Union{TendencyCallsite, UpdateStateCallsite}

""" Step `sim`ulation forward by one time step. """
function time_step!(sim::Simulation)

    start_time_step = time_ns()
    model_callbacks = Tuple(cb for cb in values(sim.callbacks) if cb isa ModelCallsite)

    if !(sim.initialized) # execute initialization step
        initialize!(sim)
        initialize!(sim.model)

        if sim.running # check that initialization didn't stop time-stepping
            if sim.verbose 
                @info "Executing initial time step..."
                start_time = time_ns()
            end

            Δt = aligned_time_step(sim, sim.time_step)
            time_step!(sim.model, Δt, callbacks=model_callbacks)

            if sim.verbose 
                elapsed_initial_step_time = prettytime(1e-9 * (time_ns() - start_time))
                @info "    ... initial time step complete ($elapsed_initial_step_time)."
            end
        else
            @warn "Simulation stopped during initialization."
        end

    else # business as usual...
        Δt = aligned_time_step(sim, sim.time_step)
        time_step!(sim.model, Δt, callbacks=model_callbacks)
    end

    for callback in values(sim.callbacks)
        callback.callsite isa TimeStepCallsite && callback.schedule(sim.model) && callback(sim)
    end

    end_time_step = time_ns()

    # Increment the wall clock
    sim.run_wall_time += 1e-9 * (end_time_step - start_time_step)

    return nothing
end

#####
##### Simulation initialization
#####

we_want_to_pickup(pickup::Bool) = pickup
we_want_to_pickup(pickup::Integer) = true
we_want_to_pickup(pickup::String) = true
we_want_to_pickup(pickup) = throw(ArgumentError("Cannot run! with pickup=$pickup"))

""" 
    initialize!(sim::Simulation, pickup=false)

Initialize a simulation:

- Update the auxiliary state of the simulation (filling halo regions, computing auxiliary fields)
- Evaluate all diagnostics, callbacks, and output writers if sim.model.clock.iteration == 0
- Add diagnostics that "depend" on output writers
"""
function initialize!(sim::Simulation)
    if sim.verbose
        @info "Initializing simulation..."
        start_time = time_ns()
    end

    model = sim.model
    clock = sim.clock

    update_state!(model)

    #=
    # Reset! the model time-stepper, evaluate all diagnostics, and write all output at first iteration
    if clock.iteration == 0
        reset!(sim.model.timestepper)

        for callback in values(sim.callbacks) 
            callback.callsite isa TimeStepCallsite && initialize!(callback, sim)
        end
    end
    =#

    if clock.iteration == 0
        for callback in values(sim.callbacks) 
            initialize!(callback, sim)
        end
    end

    sim.initialized = true

    if sim.verbose
        initialization_time = prettytime(1e-9 * (time_ns() - start_time))
        @info "    ... simulation initialization complete ($initialization_time)"
    end

    return nothing
end

