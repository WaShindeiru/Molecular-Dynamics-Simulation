use nalgebra::Vector3;

use crate::data::types::AtomType;
use crate::data::units::{R_U, VELOCITY_U};
use crate::data::ValueUnits;
use crate::persistence::json::particle_config::particle_type_file::ParticleTypeFile;
use crate::persistence::json::particle_config::Vector3Record;

#[derive(serde::Serialize)]
pub struct VelocityParticleFile {
  pub value_units: ValueUnits,
  id: usize,
  iteration: usize,
  atom_type: AtomType,
  particle_type: ParticleTypeFile,
  position: Vector3Record,
  velocity: Vector3Record,
  velocity_magnitude: f64,
}

impl VelocityParticleFile {
  pub fn new(
    id: usize,
    iteration: usize,
    atom_type: AtomType,
    particle_type: ParticleTypeFile,
    position: Vector3<f64>,
    velocity: Vector3<f64>,
  ) -> Self {
    let velocity_magnitude = velocity.magnitude();
    Self {
      value_units: ValueUnits::Unitless,
      id,
      iteration,
      atom_type,
      particle_type,
      position: Vector3Record::from(position),
      velocity: Vector3Record::from(velocity),
      velocity_magnitude,
    }
  }

  pub fn to_value_units(&self, target: ValueUnits) -> Self {
    let source = self.value_units;
    let position_scale = ValueUnits::scale_between(source, target, R_U);
    let velocity_scale = ValueUnits::scale_between(source, target, VELOCITY_U);
    Self {
      value_units: target,
      id: self.id,
      iteration: self.iteration,
      atom_type: self.atom_type,
      particle_type: self.particle_type,
      position: self.position.scale(position_scale),
      velocity: self.velocity.scale(velocity_scale),
      velocity_magnitude: self.velocity_magnitude * velocity_scale,
    }
  }
}
