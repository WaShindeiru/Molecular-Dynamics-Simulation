use nalgebra::Vector3;
use crate::data::types::AtomType;
use crate::output::AtomDTO;
use crate::particle::atom::ParticleOperations;

impl ParticleOperations for CustomPathAtom {
  fn get_id(&self) -> u64 {
    self.id
  }

  fn get_type(&self) -> &AtomType {
    &self.type_
  }

  fn get_position(&self) -> &Vector3<f64> {
    &self.position
  }

  fn get_velocity(&self) -> &Vector3<f64> {
    &self.velocity
  }

  fn get_acceleration(&self) -> &Vector3<f64> {
    &self.acceleration
  }

  fn get_mass(&self) -> f64 {
    self.mass
  }

  fn get_force(&self) -> &Vector3<f64> {
    &self.force
  }

  fn get_potential_energy(&self) -> f64 {
    self.potential_energy
  }

  fn set_velocity(&mut self, velocity_: Vector3<f64>) {
    self.velocity = velocity_;
    self.kinetic_energy = self.mass * velocity_.magnitude_squared() / 2.0;
  }

  fn set_potential_energy(&mut self, potential_energy_: f64) {
    self.potential_energy = potential_energy_;
  }

  fn set_force(&mut self, force_: Vector3<f64>) {
    self.force = force_;
  }

  fn set_acceleration(&mut self, acceleration_: Vector3<f64>) {
    self.acceleration = acceleration_;
  }

  fn update_position(&mut self, position_: Vector3<f64>) {
    self.step += 1;
    if self.step < self.path.len() {
      self.position = self.path[self.step];
    }
  }

  fn to_transfer_struct(&self) -> AtomDTO {
    AtomDTO {
      id: self.id,
      atom_type: if self.type_ == AtomType::C { 0 } else { 1 },
      x: self.position.x,
      y: self.position.y,
      z: self.position.z,
      kinetic_energy: 0.,
      potential_energy: 0.,
    }
  }

  fn custom_clone(&self) -> Box<dyn ParticleOperations> {
    Box::new(self.clone())
  }
}

#[derive(Debug, PartialEq, Clone)]
pub struct CustomPathAtom {
  id: u64,
  type_: AtomType,
  mass: f64,

  position: Vector3<f64>,
  velocity: Vector3<f64>,

  acceleration: Vector3<f64>,
  force: Vector3<f64>,

  kinetic_energy: f64,
  potential_energy: f64,

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
      type_,
      mass,
      position: path[0],
      velocity: Vector3::new(0., 0., 0.),
      force: Vector3::new(0., 0., 0.),
      acceleration: Vector3::new(0., 0., 0.),
      kinetic_energy: 0.,
      potential_energy: 0.,
      path,
      step: 0,
    }
  }
}