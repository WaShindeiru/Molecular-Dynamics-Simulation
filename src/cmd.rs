use std::io;
use std::path::PathBuf;

use clap::Parser;

use crate::simulations::run::run_from_paths;

/// Run molecular dynamics simulation from JSON config files.
#[derive(Debug, Parser)]
#[command(name = "carbon_nanotube", version, about)]
pub struct Cli {
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

pub fn run() -> io::Result<()> {
  let cli = Cli::parse();

  let simulation = cli.simulation_config.to_string_lossy().into_owned();
  let particle = cli
    .particle_or_generator
    .particles_config
    .as_ref()
    .map(|p| p.to_string_lossy().into_owned());
  let generator = cli
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
