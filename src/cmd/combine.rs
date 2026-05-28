use std::fs;
use std::io;
use std::path::PathBuf;

use crate::data::ValueUnits;
use crate::persistence::json::particle_config::ParticleConfigFile;

#[derive(Debug, clap::Args)]
pub struct CombineCommand {
  /// Path to the first initial particles JSON file.
  pub first_particles_config: PathBuf,

  /// Path to the second initial particles JSON file.
  pub second_particles_config: PathBuf,

  /// Output path for the merged initial particles JSON. Prints to stdout if omitted.
  #[arg(short, long)]
  pub output: Option<PathBuf>,
}

pub fn combine_command(command: CombineCommand) -> io::Result<()> {
  let first = read_particle_config_file(command.first_particles_config)?;
  let second = read_particle_config_file(command.second_particles_config)?;
  let merged = merge_particle_config_files(first, second).to_value_units(ValueUnits::Si);
  let json = serde_json::to_string_pretty(&merged)
    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

  if let Some(path) = command.output {
    fs::write(path, json)
  } else {
    println!("{json}");
    Ok(())
  }
}

fn read_particle_config_file(path: PathBuf) -> io::Result<ParticleConfigFile> {
  let content = fs::read_to_string(path)?;
  let config: ParticleConfigFile =
    serde_json::from_str(&content).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
  Ok(config.to_value_units(ValueUnits::Unitless))
}

fn merge_particle_config_files(
  mut first: ParticleConfigFile,
  mut second: ParticleConfigFile,
) -> ParticleConfigFile {
  let particle_id_offset = first
    .particles
    .iter()
    .map(|p| p.id)
    .max()
    .map_or(0, |id| id + 1);
  let velocity_manager_id_offset = first
    .velocity_managers
    .iter()
    .map(|vm| vm.id)
    .max()
    .map_or(0, |id| id + 1);

  for particle in second.particles.iter_mut() {
    particle.id += particle_id_offset;
    particle.velocity_manager_id = particle
      .velocity_manager_id
      .map(|id| id + velocity_manager_id_offset);
  }

  for velocity_manager in second.velocity_managers.iter_mut() {
    velocity_manager.id += velocity_manager_id_offset;
  }

  first.particles.extend(second.particles);
  first.velocity_managers.extend(second.velocity_managers);
  first.num_of_atoms += second.num_of_atoms;
  first.num_of_carbon_atoms += second.num_of_carbon_atoms;
  first.num_of_iron_atoms += second.num_of_iron_atoms;
  first.value_units = ValueUnits::Unitless;

  first
}
