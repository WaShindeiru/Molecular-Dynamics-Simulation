use nalgebra::Vector3;

use crate::data::types::AtomType;
use crate::particle::particle::ParticleKind;

#[derive(Clone)]
pub struct AtomDTO {
  pub id: usize,
  pub iteration: usize,
  pub kind: ParticleKind,
  pub atom_type: AtomType,
  pub position: Vector3<f64>,
  pub velocity: Vector3<f64>,

  pub velocity_manager_id: Option<usize>,
  pub control_velocity_manager_id: Option<usize>,

  pub kinetic_energy: f64,
  pub potential_energy: f64,
  pub thermostat_work: f64,
  pub potential_gravity_energy: f64,
  pub p_control_energy: f64,

  pub force: Vector3<f64>,
}
