use crate::data::types::AtomType;
use crate::particle::particle::ParticleKind;
use crate::persistence::dto::atom::AtomDTO;
use nalgebra::Vector3;

#[derive(Debug, PartialEq, Clone)]
pub struct CustomVelocityAtom {
  id: usize,
  iteration: usize,
  type_: AtomType,
  mass: f64,

  position: Vector3<f64>,
  velocity: Vector3<f64>,

  acceleration: Vector3<f64>,
  force: Vector3<f64>,

  kinetic_energy: f64,

  potential_energy: f64,
  previous_potential_energy: f64,
  phantom_energy: f64,

  thermostat_work: f64,
  potential_gravity_energy: f64,

  ignore_edge_conditions: bool,
}

impl CustomVelocityAtom {
  pub fn new(
    id: usize,
    type_: AtomType,
    mass: f64,
    position: Vector3<f64>,
  ) -> Self {
    CustomVelocityAtom {
      id,
      iteration: 0,
      type_,
      mass,
      position,
      velocity: Vector3::zeros(),
      force: Vector3::zeros(),
      acceleration: Vector3::zeros(),
      kinetic_energy: 0.0,
      potential_energy: 0.0,
      previous_potential_energy: 0.0,
      phantom_energy: 0.0,
      potential_gravity_energy: 0.0,
      thermostat_work: 0.0,
      ignore_edge_conditions: true,
    }
  }

  pub fn get_id(&self) -> usize {
    self.id
  }

  pub fn get_iteration(&self) -> usize {
    self.iteration
  }

  pub fn get_type(&self) -> AtomType {
    self.type_
  }

  pub fn get_position(&self) -> &Vector3<f64> {
    &self.position
  }

  pub fn get_velocity(&self) -> &Vector3<f64> {
    &self.velocity
  }

  pub fn get_acceleration(&self) -> &Vector3<f64> {
    &self.acceleration
  }

  pub fn get_mass(&self) -> f64 {
    self.mass
  }

  pub fn get_force(&self) -> &Vector3<f64> {
    &self.force
  }

  pub fn get_potential_energy(&self) -> f64 {
    self.potential_energy
  }

  pub fn get_previous_potential_energy(&self) -> f64 {
    self.previous_potential_energy
  }

  pub fn get_kinetic_energy(&self) -> f64 {
    self.kinetic_energy
  }

  pub fn get_potential_gravity_energy(&self) -> f64 {
    self.potential_gravity_energy
  }

  pub fn get_thermostat_work(&self) -> f64 {
    self.thermostat_work
  }

  pub fn ignore_edge_conditions(&self) -> bool {
    self.ignore_edge_conditions
  }

  pub fn set_velocity(&mut self, velocity_: Vector3<f64>) {
    self.velocity = velocity_;
    self.kinetic_energy = 0.0;
  }

  pub fn set_potential_energy(&mut self, potential_energy_: f64) {
    self.potential_energy = potential_energy_;
  }

  pub fn add_phantom_energy(&mut self, phantom_energy: f64) {
    self.phantom_energy += phantom_energy;
  }

  pub fn set_potential_gravity_energy(&mut self, potential_gravity_energy: f64) {
    self.potential_gravity_energy = 0.0;
  }

  pub fn set_force(&mut self, force_: Vector3<f64>) {
    self.force = force_;
  }

  pub fn set_acceleration(&mut self, acceleration_: Vector3<f64>) {
    self.acceleration = acceleration_;
  }

  pub fn set_thermostat_work(&mut self, thermostat_work: f64) {
    self.thermostat_work = 0.0;
  }

  pub fn set_iteration(&mut self, iteration_: usize) {
    self.iteration = iteration_;
  }

  pub fn update_position(&mut self, position_: Vector3<f64>) {
    self.position = position_;
  }

  pub fn reset_clone(&self) -> CustomVelocityAtom {
    let previous_potential_energy = self.potential_energy;

    CustomVelocityAtom {
      id: self.id,
      iteration: self.iteration,
      type_: self.type_,
      mass: self.mass,

      position: Vector3::zeros(),
      velocity: self.velocity,

      acceleration: Vector3::zeros(),
      force: Vector3::zeros(),

      kinetic_energy: 0.0,
      potential_energy: 0.0,
      previous_potential_energy: previous_potential_energy,
      phantom_energy: self.phantom_energy,
      thermostat_work: 0.0,
      potential_gravity_energy: 0.0,

      ignore_edge_conditions: self.ignore_edge_conditions,
    }
  }

  pub fn to_transfer_struct(&self) -> AtomDTO {
    AtomDTO {
      id: self.id,
      iteration: self.iteration,
      kind: ParticleKind::CustomVelocityAtom,
      atom_type: self.type_,
      position: self.position,
      velocity: self.velocity,
      kinetic_energy: self.kinetic_energy,
      potential_energy: self.potential_energy,
      phantom_energy: self.phantom_energy,
      thermostat_work: self.thermostat_work,
      potential_gravity_energy: self.potential_gravity_energy,
      force: self.force,
    }
  }
}
