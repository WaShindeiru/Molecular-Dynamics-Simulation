use std::fs;
use std::io;
use std::path::PathBuf;

use crate::persistence::json::particle_config::ParticleConfigFile;

#[derive(Debug, clap::Args)]
pub struct VelManReindexCommand {
  /// Path to the initial particles JSON file.
  pub particles_config: PathBuf,

  /// Id of the velocity manager to rename.
  pub old_id: usize,

  /// New id to assign to the velocity manager.
  pub new_id: usize,

  /// Output path for the reindexed particles JSON. Prints to stdout if omitted.
  #[arg(short, long)]
  pub output: Option<PathBuf>,
}

pub fn vel_man_reindex_command(command: VelManReindexCommand) -> io::Result<()> {
  let content = fs::read_to_string(command.particles_config)?;
  let mut config: ParticleConfigFile =
    serde_json::from_str(&content).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

  if command.old_id != command.new_id
    && config.velocity_managers.iter().any(|vm| vm.id == command.new_id)
  {
    return Err(io::Error::new(
      io::ErrorKind::InvalidData,
      format!("velocity manager id {} already in use", command.new_id),
    ));
  }

  let manager = config
    .velocity_managers
    .iter_mut()
    .find(|vm| vm.id == command.old_id)
    .ok_or_else(|| {
      io::Error::new(
        io::ErrorKind::InvalidData,
        format!("no velocity manager with id {} found", command.old_id),
      )
    })?;
  manager.id = command.new_id;

  for particle in config.particles.iter_mut() {
    if particle.velocity_manager_id == Some(command.old_id) {
      particle.velocity_manager_id = Some(command.new_id);
    }
  }

  let json = serde_json::to_string_pretty(&config)
    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

  if let Some(path) = command.output {
    fs::write(path, json)
  } else {
    println!("{json}");
    Ok(())
  }
}
