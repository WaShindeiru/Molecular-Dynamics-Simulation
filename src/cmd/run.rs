use std::io;
use std::path::{Path, PathBuf};

use crate::simulations::run::run_from_paths;

const PARAMETERS_FILE: &str = "parameters.json";
const PARTICLES_INITIAL_FILE: &str = "particles_initial.json";
const GENERATOR_CONFIG_FILE: &str = "generator_config.json";

#[derive(Debug, clap::Args)]
pub struct RunCommand {
  /// Path to simulation config JSON file (use with -p or -g).
  #[arg(
    short = 's',
    long,
    conflicts_with_all = ["sim_particles_dir", "sim_generator_dir"]
  )]
  pub simulation_config: Option<PathBuf>,

  /// Path to particles config JSON file (use with -s).
  #[arg(
    short = 'p',
    long,
    conflicts_with_all = ["sim_particles_dir", "sim_generator_dir", "generator_config"]
  )]
  pub particles_config: Option<PathBuf>,

  /// Path to generator config JSON file (use with -s).
  #[arg(
    short = 'g',
    long,
    conflicts_with_all = ["sim_particles_dir", "sim_generator_dir", "particles_config"]
  )]
  pub generator_config: Option<PathBuf>,

  /// Directory containing parameters.json and particles_initial.json.
  #[arg(
    long = "sp",
    conflicts_with_all = ["simulation_config", "particles_config", "generator_config", "sim_generator_dir"]
  )]
  pub sim_particles_dir: Option<PathBuf>,

  /// Directory containing parameters.json and generator_config.json.
  #[arg(
    long = "sg",
    conflicts_with_all = ["simulation_config", "particles_config", "generator_config", "sim_particles_dir"]
  )]
  pub sim_generator_dir: Option<PathBuf>,
}

pub fn run_simulation_command(command: RunCommand) -> io::Result<()> {
  let (simulation, particle, generator) = resolve_run_paths(&command)?;
  let simulation = simulation.to_string_lossy().into_owned();
  let particle = particle.map(|p| p.to_string_lossy().into_owned());
  let generator = generator.map(|g| g.to_string_lossy().into_owned());

  run_from_paths(
    simulation.as_str(),
    particle.as_deref(),
    generator.as_deref(),
  )
}

fn resolve_run_paths(
  command: &RunCommand,
) -> io::Result<(PathBuf, Option<PathBuf>, Option<PathBuf>)> {
  if let Some(dir) = &command.sim_particles_dir {
    return Ok((
      config_path(dir, PARAMETERS_FILE)?,
      Some(config_path(dir, PARTICLES_INITIAL_FILE)?),
      None,
    ));
  }

  if let Some(dir) = &command.sim_generator_dir {
    return Ok((
      config_path(dir, PARAMETERS_FILE)?,
      None,
      Some(config_path(dir, GENERATOR_CONFIG_FILE)?),
    ));
  }

  let simulation = command.simulation_config.as_ref().ok_or_else(|| {
    invalid_input(
      "Simulation config is required: use -s, --sp, or --sg.",
    )
  })?;

  match (
    command.particles_config.as_ref(),
    command.generator_config.as_ref(),
  ) {
    (Some(particles), None) => Ok((simulation.clone(), Some(particles.clone()), None)),
    (None, Some(generator)) => Ok((simulation.clone(), None, Some(generator.clone()))),
    (Some(_), Some(_)) => Err(invalid_input(
      "Exactly one of -p or -g must be provided (both were set).",
    )),
    (None, None) => Err(invalid_input(
      "Exactly one of -p or -g must be provided (none were set).",
    )),
  }
}

fn config_path(dir: &Path, file_name: &str) -> io::Result<PathBuf> {
  let path = dir.join(file_name);
  if path.is_file() {
    Ok(path)
  } else {
    Err(invalid_input(format!(
      "Expected config file `{}` in directory `{}`.",
      file_name,
      dir.display()
    )))
  }
}

fn invalid_input(message: impl Into<String>) -> io::Error {
  io::Error::new(io::ErrorKind::InvalidInput, message.into())
}
