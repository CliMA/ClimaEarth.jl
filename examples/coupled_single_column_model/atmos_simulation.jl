using ClimaAtmos

function AtmosSimulation(float_type=Float64;
                         config_type        = "column",
                         t_end              = "10days",
                         #initial_condition  = "IsothermalProfile",
                         initial_condition  = "SimplePlume",
                         surface_setup      = "DefaultMoninObukhov",
                         z_elem             = 100,
                         z_max              = 10000.0,
                         z_stretch          = false,
                         moist              = "equil",
                         dt_save_to_sol     = "Inf",
                         precip_model       = "0M",
                         vert_diff          = "FriersonDiffusion",
                         implicit_diffusion = true,
                         dt                 = "15secs",
                         job_id             = "test",
                         toml               = ["parameters.toml"],
                         hyperdiff          = "none")

    args = Dict(
        "config"             => config_type,
        "initial_condition"  => initial_condition,
        "surface_setup"      => surface_setup,
        "z_elem"             => z_elem,
        "z_max"              => z_max,
        "z_stretch"          => z_stretch,
        "moist"              => moist,
        "dt_save_to_sol"     => dt_save_to_sol,
        "precip_model"       => precip_model,
        "vert_diff"          => vert_diff,
        "implicit_diffusion" => implicit_diffusion,
        "t_end"              => t_end,
        "dt"                 => dt,
        "job_id"             => job_id,
        "toml"               => toml,
        "hyperdiff"          => hyperdiff,
        "FLOAT_TYPE"         => string(float_type)
    )

    configuration = ClimaAtmos.AtmosConfig(args)

    for (k, v) in configuration.parsed_args
        @info string(k, " = ", v)
    end

    simulation = ClimaAtmos.get_simulation(configuration)

    return simulation
end

