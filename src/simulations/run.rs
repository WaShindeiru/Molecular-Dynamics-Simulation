use std::fs;
use std::io;
use std::path::Path;

use crate::data::ValueUnits;
use crate::data::{ParticleConfig, SimulationConfig};
use crate::persistence::json::generator_config::GeneratorConfigFile;
use crate::persistence::json::particle_config::read_particle_config_from_json_file;
use crate::persistence::json::simulation_config::read_simulation_config_from_json_file;
use crate::sim_core::Engine;
use crate::simulations::generators::generate::Generator;
use crate::simulations::generators::generator_config::GeneratorConfig;
use crate::utils::logging;

pub fn run_with_particle_config(
  simulation_config: SimulationConfig,
  particle_config: ParticleConfig,
) -> io::Result<()> {
  let mut engine = Engine::from_configs(simulation_config, particle_config);
  engine.run();
  Ok(())
}

pub fn run_with_generator_config(
  simulation_config: SimulationConfig,
  generator_config: GeneratorConfig,
) -> io::Result<()> {
  let generator = generator_config.to_generator(&simulation_config);
  let particle_config = generator.generate();

  save_generator_config(&simulation_config, &generator_config)?;

  run_with_particle_config(simulation_config, particle_config)
}

pub fn run_from_paths(
  simulation_config_path: &str,
  particle_config_path: Option<&str>,
  generator_config_path: Option<&str>,
) -> io::Result<()> {
  let simulation_config = read_simulation_config_from_json_file(simulation_config_path)?;
  logging::init_logging(simulation_config.save_options.save_path.clone());

  match (particle_config_path, generator_config_path) {
    (Some(particle_path), None) => {
      let particle_config = read_particle_config_from_json_file(particle_path)?;
      run_with_particle_config(simulation_config, particle_config)
    }
    (None, Some(generator_path)) => {
      let generator_config =
        GeneratorConfigFile::from_json_file(generator_path)?.into_generator_config_unitless();
      run_with_generator_config(simulation_config, generator_config)
    }
    (Some(_), Some(_)) => Err(io::Error::new(
      io::ErrorKind::InvalidInput,
      "Exactly one of `particle_config_path` or `generator_config_path` must be provided (both were set).",
    )),
    (None, None) => Err(io::Error::new(
      io::ErrorKind::InvalidInput,
      "Exactly one of `particle_config_path` or `generator_config_path` must be provided (none were set).",
    )),
  }
}

fn save_generator_config(
  simulation_config: &SimulationConfig,
  generator_config: &GeneratorConfig,
) -> io::Result<()> {
  let save_path = simulation_config.save_options.save_path.trim();
  if save_path.is_empty() {
    return Err(io::Error::new(
      io::ErrorKind::InvalidInput,
      "Cannot save `GeneratorConfig`: `simulation_config.save_options.save_path` is empty.",
    ));
  }

  fs::create_dir_all(save_path)?;
  let generator_config_path = Path::new(save_path).join("generator_config.json");
  GeneratorConfigFile::new(generator_config.clone(), ValueUnits::Unitless)
    .to_value_units(ValueUnits::Si)
    .to_json_file(generator_config_path.to_string_lossy().as_ref())
}
