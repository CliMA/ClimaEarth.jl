#=
# Making the spaces:

@warn "perturb_initstate flag is ignored for single column configuration"
FT = eltype(params)
Δx = FT(1) # Note: This value shouldn't matter, since we only have 1 column.
quad = Quadratures.GL{1}()
horizontal_mesh = periodic_rectangle_mesh(;
    x_max = Δx,
    y_max = Δx,
    x_elem = 1,
    y_elem = 1,
)
if bubble
    @warn "Bubble correction not compatible with single column configuration. It will be switched off."
    bubble = false
end
h_space =
    make_horizontal_space(horizontal_mesh, quad, comms_ctx, bubble)
z_stretch = if parsed_args["z_stretch"]
    Meshes.GeneralizedExponentialStretching(dz_bottom, dz_top)
else
    Meshes.Uniform()
end
make_hybrid_spaces(h_space, z_max, z_elem, z_stretch; parsed_args)
=#


function AtmosSimulation(spaces,
                         config::AtmosConfig;
                         stretching)

    params = create_parameter_set(config)
    atmos = get_atmos(config, params)

    sim_info = get_sim_info(config)
    job_id = sim_info.job_id
    output_dir = sim_info.output_dir

    CP.log_parameter_information(
        config.toml_dict,
        joinpath(output_dir, "$(job_id)_parameters.toml"),
        strict = true,
    )
    YAML.write_file(joinpath(output_dir, "$job_id.yml"), config.parsed_args)

    if sim_info.restart
        s = @timed_str begin
            (Y, t_start) = get_state_restart(config)
            spaces = get_spaces_restart(Y)
        end
        @info "Allocating Y: $s"
    else
        spaces = get_spaces(config.parsed_args,
                            params,             # clima params stuff
                            config.comms_ctx)
    end

    
    initial_condition = get_initial_condition(config.parsed_args)
    surface_setup = get_surface_setup(config.parsed_args)

    if !sim_info.restart
        s = @timed_str begin
            Y = ICs.atmos_state(
                initial_condition(params),
                atmos,
                spaces.center_space,
                spaces.face_space,
            )
            t_start = Spaces.undertype(axes(Y.c))(0)
        end
        @info "Allocating Y: $s"
    end

    s = @timed_str begin
        p = build_cache(Y, atmos, params, surface_setup, sim_info)
    end
    @info "Allocating cache (p): $s"

    if config.parsed_args["discrete_hydrostatic_balance"]
        set_discrete_hydrostatic_balanced_state!(Y, p)
    end

    FT = Spaces.undertype(axes(Y.c))
    s = @timed_str begin
        ode_algo = ode_configuration(FT, config.parsed_args)
    end
    @info "ode_configuration: $s"

    s = @timed_str begin
        callback = get_callbacks(config, sim_info, atmos, params, Y, p, t_start)
    end
    @info "get_callbacks: $s"

    # Initialize diagnostics
    s = @timed_str begin
        diagnostics, writers =
            get_diagnostics(config.parsed_args, atmos, Spaces.axes(Y.c))
    end
    @info "initializing diagnostics: $s"

    length(diagnostics) > 0 && @info "Computing diagnostics:"

    # First, we convert all the ScheduledDiagnosticTime into ScheduledDiagnosticIteration,
    # ensuring that there is consistency in the timestep and the periods and translating
    # those periods that depended on the timestep
    diagnostics_iterations =
        [CAD.ScheduledDiagnosticIterations(d, sim_info.dt) for d in diagnostics]

    # For diagnostics that perform reductions, the storage is used for the values computed
    # at each call. Reductions also save the accumulated value in diagnostic_accumulators.
    diagnostic_storage = Dict()
    diagnostic_accumulators = Dict()
    diagnostic_counters = Dict()

    s = @timed_str begin
        diagnostics_functions = CAD.get_callbacks_from_diagnostics(
            diagnostics_iterations,
            diagnostic_storage,
            diagnostic_accumulators,
            diagnostic_counters,
            output_dir,
        )
    end
    @info "Prepared diagnostic callbacks: $s"

    # It would be nice to just pass the callbacks to the integrator. However, this leads to
    # a significant increase in compile time for reasons that are not known. For this
    # reason, we only add one callback to the integrator, and this function takes care of
    # executing the other callbacks. This single function is orchestrate_diagnostics

    orchestrate_diagnostics(integrator) =
        CAD.orchestrate_diagnostics(integrator, diagnostics_functions)

    diagnostic_callbacks =
        call_every_n_steps(orchestrate_diagnostics, skip_first = true)

    # The generic constructor for SciMLBase.CallbackSet has to split callbacks into discrete
    # and continuous. This is not hard, but can introduce significant latency. However, all
    # the callbacks in ClimaAtmos are discrete_callbacks, so we directly pass this
    # information to the constructor
    continuous_callbacks = tuple()
    discrete_callbacks = (callback..., diagnostic_callbacks)

    s = @timed_str begin
        all_callbacks =
            SciMLBase.CallbackSet(continuous_callbacks, discrete_callbacks)
    end
    @info "Prepared SciMLBase.CallbackSet callbacks: $s"
    steps_cycle_non_diag = n_steps_per_cycle_per_cb(all_callbacks, sim_info.dt)
    steps_cycle_diag =
        n_steps_per_cycle_per_cb_diagnostic(diagnostics_functions)
    steps_cycle = lcm([steps_cycle_non_diag..., steps_cycle_diag...])
    @info "n_steps_per_cycle_per_cb (non diagnostics): $steps_cycle_non_diag"
    @info "n_steps_per_cycle_per_cb_diagnostic: $steps_cycle_diag"
    @info "n_steps_per_cycle (non diagnostics): $steps_cycle"

    tspan = (t_start, sim_info.t_end)
    s = @timed_str begin
        integrator_args, integrator_kwargs = args_integrator(
            config.parsed_args,
            Y,
            p,
            tspan,
            ode_algo,
            all_callbacks,
        )
    end

    s = @timed_str begin
        integrator = SciMLBase.init(integrator_args...; integrator_kwargs...)
    end
    @info "init integrator: $s"
    reset_graceful_exit(output_dir)

    s = @timed_str init_diagnostics!(
        diagnostics_iterations,
        diagnostic_storage,
        diagnostic_accumulators,
        diagnostic_counters,
        output_dir,
        Y,
        p,
        t_start;
        warn_allocations = config.parsed_args["warn_allocations_diagnostics"],
    )
    @info "Init diagnostics: $s"

    return AtmosSimulation(
        job_id,
        output_dir,
        sim_info.start_date,
        sim_info.t_end,
        writers,
        integrator,
    )
end
