use std::fs;
use std::io;
use std::path::PathBuf;

use crate::persistence::json::particle_config::ParticleConfigFile;

#[derive(Debug, clap::Args)]
pub struct InspectCommand {
  /// Path to the initial particles JSON file.
  pub particles_config: PathBuf,
}

pub fn inspect_command(command: InspectCommand) -> io::Result<()> {
  let content = fs::read_to_string(command.particles_config)?;
  let config: ParticleConfigFile =
    serde_json::from_str(&content).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
  let particle_config = config.try_into_runtime_unitless()?;

  let mut min = nalgebra::Vector3::from_element(f64::INFINITY);
  let mut max = nalgebra::Vector3::from_element(f64::NEG_INFINITY);

  for particle in &particle_config.atoms {
    let position = particle.get_position();
    min = min.zip_map(position, f64::min);
    max = max.zip_map(position, f64::max);
  }

  println!("x: min = {}, max = {}", min.x, max.x);
  println!("y: min = {}, max = {}", min.y, max.y);
  println!("z: min = {}, max = {}", min.z, max.z);

  Ok(())
}
