pub mod particle_type_file;

use std::{fs, io};

use nalgebra::Vector3;

use crate::data::constants::{ATOMIC_MASS_C, ATOMIC_MASS_FE};
use crate::data::types::AtomType;
use crate::data::units::{R_U, VELOCITY_U};
use crate::data::{ParticleConfig, ValueUnits};
use crate::particle::{Atom, Particle};
use crate::persistence::json::particle_config::particle_type_file::ParticleTypeFile;

#[derive(Clone, serde::Serialize, serde::Deserialize)]
pub struct ParticleConfigFile {
  #[serde(default)]
  pub value_units: ValueUnits,
  pub particles: Vec<ParticleInitialState>,
  pub num_of_atoms: usize,
  pub num_of_carbon_atoms: usize,
  pub num_of_iron_atoms: usize,
}

impl ParticleConfigFile {
  pub fn from_runtime(config: &ParticleConfig, target_units: ValueUnits) -> Self {
    let particles = config
      .atoms
      .iter()
      .map(ParticleInitialState::from_runtime)
      .collect();

    let unitless = Self {
      value_units: ValueUnits::Unitless,
      particles,
      num_of_atoms: config.num_of_atoms,
      num_of_carbon_atoms: config.num_of_carbon_atoms,
      num_of_iron_atoms: config.num_of_iron_atoms,
    };

    unitless.to_value_units(target_units)
  }

  pub fn try_into_runtime_unitless(self) -> io::Result<ParticleConfig> {
    let unitless = self.to_value_units(ValueUnits::Unitless);
    debug_assert_eq!(unitless.value_units, ValueUnits::Unitless);

    let atoms = unitless
      .particles
      .iter()
      .map(ParticleInitialState::try_to_runtime)
      .collect::<io::Result<Vec<_>>>()?;

    Ok(ParticleConfig::new(atoms))
  }

  pub fn to_value_units(&self, target: ValueUnits) -> Self {
    let source = self.value_units;
    let position_scale = ValueUnits::scale_between(source, target, R_U);
    let velocity_scale = ValueUnits::scale_between(source, target, VELOCITY_U);

    let particles = self
      .particles
      .iter()
      .cloned()
      .map(|particle| particle.to_value_units(position_scale, velocity_scale))
      .collect();

    Self {
      value_units: target,
      particles,
      num_of_atoms: self.num_of_atoms,
      num_of_carbon_atoms: self.num_of_carbon_atoms,
      num_of_iron_atoms: self.num_of_iron_atoms,
    }
  }
}

#[derive(Clone, serde::Serialize, serde::Deserialize)]
pub struct ParticleInitialState {
  id: usize,
  atom_type: AtomType,
  particle_type: ParticleTypeFile,
  position: Vector3Record,
  velocity: Vector3Record,
}

impl ParticleInitialState {
  fn from_runtime(particle: &Particle) -> Self {
    let position = particle.get_position();
    let velocity = particle.get_velocity();

    Self {
      id: particle.get_id(),
      atom_type: particle.get_type(),
      particle_type: ParticleTypeFile::from_runtime(particle),
      position: Vector3Record {
        x: position.x,
        y: position.y,
        z: position.z,
      },
      velocity: Vector3Record {
        x: velocity.x,
        y: velocity.y,
        z: velocity.z,
      },
    }
  }

  fn try_to_runtime(&self) -> io::Result<Particle> {
    let mass = match self.atom_type {
      AtomType::C => ATOMIC_MASS_C,
      AtomType::Fe => ATOMIC_MASS_FE,
    };

    Ok(match self.particle_type {
      ParticleTypeFile::Atom => Particle::Atom(Atom::new(
        self.id,
        self.atom_type,
        mass,
        self.position.to_runtime(),
        self.velocity.to_runtime(),
        Vector3::zeros(),
        Vector3::zeros(),
        0.0,
      )),
      ParticleTypeFile::CustomPathAtom => {
        return Err(io::Error::new(
          io::ErrorKind::InvalidData,
          "Cannot read CustomPathAtom from ParticleInitialState: path data is not stored in this JSON format.",
        ));
      }
    })
  }

  fn to_value_units(&self, position_scale: f64, velocity_scale: f64) -> Self {
    Self {
      id: self.id,
      atom_type: self.atom_type,
      particle_type: self.particle_type,
      position: self.position.scale(position_scale),
      velocity: self.velocity.scale(velocity_scale),
    }
  }
}

#[derive(Clone, serde::Serialize, serde::Deserialize)]
pub struct Vector3Record {
  x: f64,
  y: f64,
  z: f64,
}

impl Vector3Record {
  fn to_runtime(&self) -> Vector3<f64> {
    Vector3::new(self.x, self.y, self.z)
  }

  fn scale(&self, scale: f64) -> Self {
    Self {
      x: self.x * scale,
      y: self.y * scale,
      z: self.z * scale,
    }
  }
}

pub fn particle_config_to_initial_json_string(
  config: &ParticleConfig,
) -> Result<String, serde_json::Error> {
  serde_json::to_string_pretty(&ParticleConfigFile::from_runtime(config, ValueUnits::Si))
}

pub fn particle_config_to_initial_json_file(config: &ParticleConfig, path: &str) -> io::Result<()> {
  let json = particle_config_to_initial_json_string(config)
    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
  fs::write(path, json)
}

pub fn read_particle_config_from_json_str(content: &str) -> io::Result<ParticleConfig> {
  let config_file: ParticleConfigFile =
    serde_json::from_str(content).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
  config_file.try_into_runtime_unitless()
}

pub fn read_particle_config_from_json_file(path: &str) -> io::Result<ParticleConfig> {
  let content = fs::read_to_string(path)?;
  read_particle_config_from_json_str(&content)
}
