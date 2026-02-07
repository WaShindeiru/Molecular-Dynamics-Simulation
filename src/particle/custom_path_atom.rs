use nalgebra::Vector3;
use crate::data::types::AtomType;
use crate::output::AtomDTO;

impl CustomPathAtom {
  pub fn get_id(&self) -> u64 {
    self.id
  }
  
  pub fn get_iteration(&self) -> usize {
    self.iteration
  }

  pub fn get_type(&self) -> &AtomType {
    &self.type_
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

  pub fn get_thermostat_work(&self) -> f64 { self.thermostat_work }

  pub fn set_velocity(&mut self, velocity_: Vector3<f64>) {
    self.velocity = velocity_;
    self.kinetic_energy = self.mass * velocity_.magnitude_squared() / 2.0;
  }

  pub fn set_potential_energy(&mut self, potential_energy_: f64) {
    self.potential_energy = potential_energy_;
  }

  pub fn set_force(&mut self, force_: Vector3<f64>) {
    self.force = force_;
  }

  pub fn set_acceleration(&mut self, acceleration_: Vector3<f64>) {
    self.acceleration = acceleration_;
  }

  pub fn set_thermostat_work(&mut self, thermostat_work: f64) {
    self.thermostat_work = thermostat_work
  }

  pub fn set_iteration(&mut self, iteration_: usize) {
    self.iteration = iteration_;
  }

  pub fn update_position(&mut self, position_: Vector3<f64>) {
    self.step += 1;
    if self.step < self.path.len() {
      self.position = self.path[self.step];
    }
  }

  pub fn to_transfer_struct(&self) -> AtomDTO {
    AtomDTO {
      id: self.id,
      iteration: self.iteration,
      atom_type: if self.type_ == AtomType::C { 0 } else { 1 },
      x: self.position.x,
      y: self.position.y,
      z: self.position.z,
      kinetic_energy: self.kinetic_energy,
      potential_energy: self.potential_energy,
      thermostat_work: self.thermostat_work,
      force_x: self.force.x,
      force_y: self.force.y,
      force_z: self.force.z,
    }
  }
}

#[derive(Debug, PartialEq, Clone)]
pub struct CustomPathAtom {
  id: u64,
  iteration: usize,
  type_: AtomType,
  mass: f64,

  position: Vector3<f64>,
  velocity: Vector3<f64>,

  acceleration: Vector3<f64>,
  force: Vector3<f64>,

  kinetic_energy: f64,
  potential_energy: f64,
  thermostat_work: f64,

  path: Vec<Vector3<f64>>,
  step: usize,
}

impl CustomPathAtom {
  pub fn new(
    id: u64,
    type_: AtomType,
    mass: f64,
    path: Vec<Vector3<f64>>,
  ) -> Self {
    CustomPathAtom {
      id,
      iteration: 0,
      type_,
      mass,
      position: path[0],
      velocity: Vector3::new(0., 0., 0.),
      force: Vector3::new(0., 0., 0.),
      acceleration: Vector3::new(0., 0., 0.),
      kinetic_energy: 0.,
      potential_energy: 0.,
      thermostat_work: 0.0,
      path,
      step: 0,
    }
  }

  pub fn new_custom_iteration(
    id: u64,
    iteration: usize,
    type_: AtomType,
    mass: f64,
    path: Vec<Vector3<f64>>,
  ) -> Self {
    CustomPathAtom {
      id,
      iteration,
      type_,
      mass,
      position: path[0],
      velocity: Vector3::new(0., 0., 0.),
      force: Vector3::new(0., 0., 0.),
      acceleration: Vector3::new(0., 0., 0.),
      kinetic_energy: 0.,
      potential_energy: 0.,
      thermostat_work: 0.0,
      path,
      step: 0,
    }
  }
}