use std::fs;
use std::io;
use std::path::PathBuf;

use nalgebra::Vector3;

use crate::data::units::R_U;
use crate::data::ValueUnits;
use crate::persistence::json::particle_config::ParticleConfigFile;

#[derive(Debug, clap::Args)]
#[command(name = "move-minus")]
pub struct MoveMinusCommand {
  /// Path to the initial particles JSON file.
  pub particles_config: PathBuf,

  /// Displacement to subtract along the X axis in meters (SI).
  #[arg(short = 'x', long, default_value_t = 0.0)]
  pub x: f64,

  /// Displacement to subtract along the Y axis in meters (SI).
  #[arg(short = 'y', long, default_value_t = 0.0)]
  pub y: f64,

  /// Displacement to subtract along the Z axis in meters (SI).
  #[arg(short = 'z', long, default_value_t = 0.0)]
  pub z: f64,

  /// Output path for the moved particles JSON. Prints to stdout if omitted.
  #[arg(short, long)]
  pub output: Option<PathBuf>,
}

pub fn move_particles_minus_command(command: MoveMinusCommand) -> io::Result<()> {
  let content = fs::read_to_string(command.particles_config)?;
  let config: ParticleConfigFile =
    serde_json::from_str(&content).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

  let offset_unitless = -Vector3::new(command.x, command.y, command.z) / R_U;

  let moved = config
    .to_value_units(ValueUnits::Unitless)
    .translate(offset_unitless)
    .to_value_units(ValueUnits::Si);

  let json = serde_json::to_string_pretty(&moved)
    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

  if let Some(path) = command.output {
    fs::write(path, json)
  } else {
    println!("{json}");
    Ok(())
  }
}
