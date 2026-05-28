use std::io;
use std::path::PathBuf;

use nalgebra::Vector3;

use crate::data::ValueUnits;
use crate::persistence::json::generator_config::GeneratorConfigFile;
use crate::persistence::transformation::{
  nanotube_txt_file_to_particles, particles_to_nanotube_generator_config,
  particles_to_velocity_nanotube_generator_config,
};
use crate::simulations::generators::core::generator_config::GeneratorConfig;
use crate::simulations::generators::core::generator_config::nanotube::NanotubeGeneratorParticleFile;

#[derive(Debug, clap::Args)]
pub struct ParseCommand {
  #[command(flatten)]
  pub input: ParseInput,

  #[command(flatten)]
  pub output: ParseOutput,
}

#[derive(Debug, clap::Args)]
#[group(required = true, multiple = false)]
pub struct ParseInput {
  /// Path to raw nanotube TXT definition.
  #[arg(short = 'r', long)]
  pub raw: Option<PathBuf>,

  /// Path to existing generator config JSON definition.
  #[arg(short = 'j', long)]
  pub json: Option<PathBuf>,
}

#[derive(Debug, clap::Args)]
#[group(required = true, multiple = true)]
pub struct ParseOutput {
  /// Output path for NanotubeGenerator config JSON.
  #[arg(short = 'b', long)]
  pub basic: Option<PathBuf>,

  /// Output path for VelocityNanotubeGenerator config JSON.
  #[arg(short = 'v', long)]
  pub velocity: Option<PathBuf>,
}

pub fn parse_command(command: ParseCommand) -> io::Result<()> {
  let particles = read_parse_particles(command.input)?;

  if let Some(path) = command.output.basic {
    let config = particles_to_nanotube_generator_config(
      particles.clone(),
      Vector3::zeros(),
      Vector3::zeros(),
      Vector3::zeros(),
    );
    GeneratorConfigFile::new(GeneratorConfig::Nanotube(config), ValueUnits::Unitless)
      .to_json_file(path.to_string_lossy().as_ref())?;
  }

  if let Some(path) = command.output.velocity {
    let config = particles_to_velocity_nanotube_generator_config(
      particles,
      Vector3::zeros(),
      vec![(0, Vector3::zeros())],
    );
    GeneratorConfigFile::new(
      GeneratorConfig::VelocityNanotube(config),
      ValueUnits::Unitless,
    )
    .to_json_file(path.to_string_lossy().as_ref())?;
  }

  Ok(())
}

fn read_parse_particles(input: ParseInput) -> io::Result<Vec<NanotubeGeneratorParticleFile>> {
  match (input.raw, input.json) {
    (Some(path), None) => nanotube_txt_file_to_particles(path),
    (None, Some(path)) => read_particles_from_generator_json(path),
    (Some(_), Some(_)) => Err(io::Error::new(
      io::ErrorKind::InvalidInput,
      "Exactly one of `--raw` or `--json` must be provided.",
    )),
    (None, None) => Err(io::Error::new(
      io::ErrorKind::InvalidInput,
      "One of `--raw` or `--json` must be provided.",
    )),
  }
}

fn read_particles_from_generator_json(
  path: PathBuf,
) -> io::Result<Vec<NanotubeGeneratorParticleFile>> {
  match GeneratorConfigFile::from_json_file(path.to_string_lossy().as_ref())?
    .into_generator_config_unitless()
  {
    GeneratorConfig::Nanotube(config) => Ok(config.particles),
    GeneratorConfig::VelocityNanotube(config) => Ok(config.particles),
    GeneratorConfig::Dense(_) => Err(io::Error::new(
      io::ErrorKind::InvalidData,
      "Dense generator config does not contain nanotube particles.",
    )),
  }
}
