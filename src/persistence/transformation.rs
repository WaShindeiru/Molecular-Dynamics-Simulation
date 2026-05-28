use std::path::Path;
use std::{fs, io};

use nalgebra::Vector3;

use crate::data::types::AtomType;
use crate::persistence::json::particle_config::Vector3Record;
use crate::persistence::json::velocity_manager_file::VelocityChangeEntry;
use crate::simulations::generators::core::generator_config::nanotube::{
  NanotubeGeneratorConfig, NanotubeGeneratorParticleFile,
};
use crate::simulations::generators::core::generator_config::velocity_nanotube::VelocityNanotubeGeneratorConfig;

pub fn nanotube_txt_file_to_particles<P: AsRef<Path>>(
  path: P,
) -> io::Result<Vec<NanotubeGeneratorParticleFile>> {
  let content = fs::read_to_string(path)?;
  nanotube_txt_definition_to_particles(&content)
}

pub fn nanotube_txt_definition_to_particles(
  content: &str,
) -> io::Result<Vec<NanotubeGeneratorParticleFile>> {
  let mut lines = content.lines();
  let atom_count = lines
    .next()
    .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Nanotube TXT file is empty"))?
    .trim()
    .parse::<usize>()
    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

  let mut particles = Vec::with_capacity(atom_count);

  for line in lines.take(atom_count) {
    let parts: Vec<_> = line.split_whitespace().collect();
    if parts.len() < 4 {
      return Err(io::Error::new(
        io::ErrorKind::InvalidData,
        format!("Invalid nanotube particle line: `{line}`"),
      ));
    }

    particles.push(NanotubeGeneratorParticleFile {
      particle_type: nanotube_atom_type(parts[0]),
      position: Vector3Record::from(Vector3::new(
        parse_coordinate(parts[1])?,
        parse_coordinate(parts[2])?,
        parse_coordinate(parts[3])?,
      )),
    });
  }

  if particles.len() != atom_count {
    return Err(io::Error::new(
      io::ErrorKind::InvalidData,
      format!(
        "Nanotube TXT file declared {atom_count} atoms but contained {} particle lines",
        particles.len()
      ),
    ));
  }

  Ok(particles)
}

pub fn particles_to_nanotube_generator_config(
  particles: Vec<NanotubeGeneratorParticleFile>,
  offset: Vector3<f64>,
  vel_mean: Vector3<f64>,
  vel_std_dev: Vector3<f64>,
) -> NanotubeGeneratorConfig {
  NanotubeGeneratorConfig::new(
    particles,
    Vector3Record::from(offset),
    Vector3Record::from(vel_mean),
    Vector3Record::from(vel_std_dev),
  )
}

pub fn particles_to_velocity_nanotube_generator_config(
  particles: Vec<NanotubeGeneratorParticleFile>,
  offset: Vector3<f64>,
  velocities: Vec<(usize, Vector3<f64>)>,
) -> VelocityNanotubeGeneratorConfig {
  VelocityNanotubeGeneratorConfig::new(
    particles,
    Vector3Record::from(offset),
    velocities
      .into_iter()
      .map(|(iteration, velocity)| VelocityChangeEntry::from_runtime(iteration, velocity))
      .collect(),
  )
}

fn nanotube_atom_type(atom: &str) -> AtomType {
  match atom {
    "Fe" => AtomType::Fe,
    _ => AtomType::C_nanotube,
  }
}

fn parse_coordinate(value: &str) -> io::Result<f64> {
  value
    .parse::<f64>()
    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}
