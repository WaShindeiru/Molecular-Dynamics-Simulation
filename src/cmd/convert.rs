use std::fs;
use std::io;
use std::path::PathBuf;

use crate::data::ValueUnits;
use crate::persistence::json::particle_config::ParticleConfigFile;

#[derive(Debug, clap::Args)]
pub struct ConvertCommand {
  /// Path to the initial particles JSON file.
  pub particles_config: PathBuf,

  /// Output path for the converted particles JSON. Prints to stdout if omitted.
  #[arg(short, long)]
  pub output: Option<PathBuf>,
}

pub fn convert_command(command: ConvertCommand) -> io::Result<()> {
  let content = fs::read_to_string(command.particles_config)?;
  let config: ParticleConfigFile =
    serde_json::from_str(&content).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

  let converted = config
    .to_value_units(ValueUnits::Unitless)
    .convert_to_custom_velocity_atoms()
    .to_value_units(ValueUnits::Si);

  let json = serde_json::to_string_pretty(&converted)
    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

  if let Some(path) = command.output {
    fs::write(path, json)
  } else {
    println!("{json}");
    Ok(())
  }
}
