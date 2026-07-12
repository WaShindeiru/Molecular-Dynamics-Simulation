use std::fs;
use std::io;
use std::path::PathBuf;

use crate::persistence::json::particle_config::ParticleConfigFile;

#[derive(Debug, clap::Args)]
pub struct ReindexCommand {
  /// Path to the initial particles JSON file.
  pub particles_config: PathBuf,

  /// Output path for the reindexed particles JSON. Prints to stdout if omitted.
  #[arg(short, long)]
  pub output: Option<PathBuf>,
}

pub fn reindex_command(command: ReindexCommand) -> io::Result<()> {
  let content = fs::read_to_string(command.particles_config)?;
  let mut config: ParticleConfigFile =
    serde_json::from_str(&content).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

  for (new_id, particle) in config.particles.iter_mut().enumerate() {
    particle.id = new_id;
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
