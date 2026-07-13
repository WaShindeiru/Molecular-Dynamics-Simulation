pub mod particle_type_file;

use std::{fs, io};

use nalgebra::Vector3;

use crate::data::constants::{ATOMIC_MASS_C, ATOMIC_MASS_FE};
use crate::data::types::AtomType;
use crate::data::units::{R_U, VELOCITY_U};
use crate::data::{ParticleConfig, ValueUnits};
use crate::particle::{Atom, Particle, CustomVelocityAtom};
use crate::persistence::json::particle_config::particle_type_file::ParticleTypeFile;
use crate::persistence::json::velocity_manager_file::VelocityManagerFile;

#[derive(Clone, serde::Serialize, serde::Deserialize)]
pub struct ParticleConfigFile {
  #[serde(default)]
  pub value_units: ValueUnits,
  pub particles: Vec<ParticleInitialState>,
  pub num_of_atoms: usize,
  pub num_of_carbon_atoms: usize,
  pub num_of_iron_atoms: usize,
  #[serde(default)]
  pub velocity_managers: Vec<VelocityManagerFile>,
}

impl ParticleConfigFile {
  pub fn from_runtime(config: &ParticleConfig, target_units: ValueUnits) -> Self {
    let particles = config
      .atoms
      .iter()
      .map(|particle| {
        let mut particle_file = ParticleInitialState::from_runtime(particle);
        if let Particle::CustomVelocityAtom(p) = particle {
          particle_file.velocity_manager_id = Some(p.get_particle_velocity_manager_id());
        }
        particle_file
      })
      .collect();

    let velocity_managers = config
      .velocity_schedules
      .iter()
      .map(VelocityManagerFile::from_schedule)
      .collect();

    let unitless = Self {
      value_units: ValueUnits::Unitless,
      particles,
      num_of_atoms: config.num_of_atoms,
      num_of_carbon_atoms: config.num_of_carbon_atoms,
      num_of_iron_atoms: config.num_of_iron_atoms,
      velocity_managers,
    };

    unitless.to_value_units(target_units)
  }

  pub fn try_into_runtime_unitless(self) -> io::Result<ParticleConfig> {
    let unitless = self.to_value_units(ValueUnits::Unitless);
    debug_assert_eq!(unitless.value_units, ValueUnits::Unitless);

    let velocity_schedules = unitless
      .velocity_managers
      .iter()
      .map(VelocityManagerFile::to_schedule)
      .collect();

    let atoms = unitless
      .particles
      .iter()
      .map(ParticleInitialState::try_to_runtime)
      .collect::<io::Result<Vec<_>>>()?;

    Ok(ParticleConfig::new_with_schedules(atoms, velocity_schedules))
  }

  pub fn convert_to_custom_velocity_atoms(mut self) -> Self {
    use crate::persistence::json::velocity_manager_file::{VelocityChangeEntry, VelocityManagerFile};

    const SHARED_VELOCITY_MANAGER_ID: usize = 0;

    self.velocity_managers = vec![VelocityManagerFile {
      id: SHARED_VELOCITY_MANAGER_ID,
      velocities: vec![VelocityChangeEntry::from_runtime(0, Vector3::zeros())],
    }];
    self.particles = self
      .particles
      .into_iter()
      .map(|p| p.into_custom_velocity_atom(SHARED_VELOCITY_MANAGER_ID))
      .collect();
    self
  }

  pub fn translate(mut self, offset: Vector3<f64>) -> Self {
    self.particles = self.particles.into_iter().map(|p| p.translate(offset)).collect();
    self
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

    let velocity_managers = self
      .velocity_managers
      .iter()
      .map(|vm| vm.scale_velocities(velocity_scale))
      .collect();

    Self {
      value_units: target,
      particles,
      num_of_atoms: self.num_of_atoms,
      num_of_carbon_atoms: self.num_of_carbon_atoms,
      num_of_iron_atoms: self.num_of_iron_atoms,
      velocity_managers,
    }
  }
}

#[derive(Clone, serde::Serialize, serde::Deserialize)]
pub struct ParticleInitialState {
  pub id: usize,
  atom_type: AtomType,
  particle_type: ParticleTypeFile,
  position: Vector3Record,
  velocity: Vector3Record,
  #[serde(skip_serializing_if = "Option::is_none")]
  #[serde(default)]
  pub velocity_manager_id: Option<usize>,
}

impl ParticleInitialState {
  pub fn new(
    id: usize,
    atom_type: AtomType,
    particle_type: ParticleTypeFile,
    position: Vector3<f64>,
    velocity: Vector3<f64>,
    velocity_manager_id: Option<usize>,
  ) -> Self {
    Self {
      id,
      atom_type,
      particle_type,
      position: Vector3Record::from(position),
      velocity: Vector3Record::from(velocity),
      velocity_manager_id,
    }
  }

  pub fn from_runtime(particle: &Particle) -> Self {
    Self::new(
      particle.get_id(),
      particle.get_type(),
      ParticleTypeFile::from_runtime(particle),
      *particle.get_position(),
      *particle.get_velocity(),
      None,
    )
  }

  pub fn try_to_runtime(&self) -> io::Result<Particle> {
    let mass = match self.atom_type {
      AtomType::C | AtomType::C_nanotube => ATOMIC_MASS_C,
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
      ParticleTypeFile::CustomVelocityAtom => {
        let particle_velocity_manager_id = self.velocity_manager_id.ok_or_else(|| {
          io::Error::new(
            io::ErrorKind::InvalidData,
            format!("CustomVelocityAtom id={} is missing velocity_manager_id", self.id),
          )
        })?;
        Particle::CustomVelocityAtom(CustomVelocityAtom::new(
          self.id,
          self.atom_type,
          mass,
          self.position.to_runtime(),
          particle_velocity_manager_id,
        ))
      }
    })
  }

  pub fn into_custom_velocity_atom(self, velocity_manager_id: usize) -> Self {
    Self {
      particle_type: ParticleTypeFile::CustomVelocityAtom,
      velocity: Vector3Record::from(Vector3::zeros()),
      velocity_manager_id: Some(velocity_manager_id),
      ..self
    }
  }

  pub fn translate(self, offset: Vector3<f64>) -> Self {
    Self {
      position: Vector3Record::from(self.position.to_runtime() + offset),
      ..self
    }
  }

  fn to_value_units(&self, position_scale: f64, velocity_scale: f64) -> Self {
    Self {
      id: self.id,
      atom_type: self.atom_type,
      particle_type: self.particle_type,
      position: self.position.scale(position_scale),
      velocity: self.velocity.scale(velocity_scale),
      velocity_manager_id: self.velocity_manager_id,
    }
  }
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct Vector3Record {
  x: f64,
  y: f64,
  z: f64,
}

impl From<Vector3<f64>> for Vector3Record {
  fn from(v: Vector3<f64>) -> Self {
    Self { x: v.x, y: v.y, z: v.z }
  }
}

impl Vector3Record {
  pub fn to_runtime(&self) -> Vector3<f64> {
    Vector3::new(self.x, self.y, self.z)
  }

  pub fn scale(&self, scale: f64) -> Self {
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
