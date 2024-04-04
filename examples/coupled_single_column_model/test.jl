using ClimaAtmos
configuration_dict = Dict("config" => "column",
                          "dt_save_to_sol" => "Inf",
                          "z_max" => 30000.0)
configuration = ClimaAtmos.AtmosConfig(configuration_dict)

