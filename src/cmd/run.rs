use std::io;
use std::path::PathBuf;

use crate::simulations::run::run_from_paths;

#[derive(Debug, clap::Args)]
pub struct RunCommand {
  /// Path to simulation config JSON file.
  #[arg(short = 's', long)]
  pub simulation_config: PathBuf,

  #[command(flatten)]
  pub particle_or_generator: ParticleOrGenerator,
}

#[derive(Debug, clap::Args)]
#[group(required = true, multiple = false)]
pub struct ParticleOrGenerator {
  /// Path to particles config JSON file.
  #[arg(short = 'p', long)]
  pub particles_config: Option<PathBuf>,

  /// Path to generator config JSON file.
  #[arg(short = 'g', long)]
  pub generator_config: Option<PathBuf>,
}

pub fn run_simulation_command(command: RunCommand) -> io::Result<()> {
  let simulation = command.simulation_config.to_string_lossy().into_owned();
  let particle = command
    .particle_or_generator
    .particles_config
    .as_ref()
    .map(|p| p.to_string_lossy().into_owned());
  let generator = command
    .particle_or_generator
    .generator_config
    .as_ref()
    .map(|p| p.to_string_lossy().into_owned());

  run_from_paths(
    simulation.as_str(),
    particle.as_deref(),
    generator.as_deref(),
  )
}
