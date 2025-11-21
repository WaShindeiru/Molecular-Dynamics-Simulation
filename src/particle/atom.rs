use std::sync::Mutex;
use nalgebra::base::Vector3;
use rand::distributions::{Distribution, Uniform};
use crate::data::constants::{ATOMIC_MASS_C, ATOMIC_MASS_FE};
use crate::data::types::AtomType;
use crate::output::{AtomDTO};

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Atom {
  id: u64,
  type_: AtomType,
  mass: f64,

  position: Vector3<f64>,
  velocity: Vector3<f64>,

  acceleration: Vector3<f64>,
  force: Vector3<f64>,

  kinetic_energy: f64,
  potential_energy: f64,
}

impl Atom {
  pub fn new(
    id: u64,
    type_: AtomType,
    mass: f64,
    position: Vector3<f64>,
    velocity: Vector3<f64>,
    force: Vector3<f64>,
    acceleration: Vector3<f64>,
    potential_energy: f64,
    ) -> Self {
    Atom {
      id,
      type_,
      mass,
      position,
      velocity,
      force,
      acceleration,
      kinetic_energy: mass * velocity.magnitude_squared() / 2.0,
      potential_energy,
    }
  }
  
  pub fn get_id(&self) -> u64 {
    self.id
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
  
  pub fn get_force(&self) -> &Vector3<f64> {
    &self.force
  }
  
  pub fn get_acceleration(&self) -> &Vector3<f64> {
    &self.acceleration
  }

  pub fn get_mass(&self) -> f64 {
    self.mass
  }

  pub fn set_velocity(&mut self, velocity_: Vector3<f64>) {
    self.velocity = velocity_;
    self.kinetic_energy = self.mass * velocity_.magnitude_squared() / 2.0;
  }

  pub fn set_potential_energy(&mut self, potential_energy_: f64) {
    self.potential_energy = potential_energy_;
  }

  pub fn get_potential_energy(&self) -> f64 {
    self.potential_energy
  }

  pub fn set_position(&mut self, position_: Vector3<f64>) {
    self.position = position_;
  }

  pub fn to_transfer_struct(&self) -> AtomDTO {
    AtomDTO {
      id: self.id,
      atom_type: if self.type_ == AtomType::C { 0 } else { 1 },
      x: self.position.x,
      y: self.position.y,
      z: self.position.z,
      kinetic_energy: self.kinetic_energy,
      potential_energy: self.potential_energy,
    }
  }
}

struct AtomFactory {
  counter: u64
}

impl AtomFactory {
  fn new() -> Self {
    AtomFactory{
      counter: 0,
    }
  }

  fn get_atom(&mut self, atom: AtomType, position_: Vector3<f64>, velocity_: Vector3<f64>) -> Atom {
    let result = match atom {
      AtomType::C => Atom{
        id: self.counter,
        type_: AtomType::C,
        mass: ATOMIC_MASS_C,
        position: position_,
        velocity: velocity_,
        force: Vector3::new(0.0, 0.0, 0.0),
        acceleration: Vector3::new(0.0, 0.0, 0.0),
        kinetic_energy: ATOMIC_MASS_C * velocity_.magnitude_squared() / 2.0,
        potential_energy: 0.0,
      },
      AtomType::Fe => Atom{
        type_: AtomType::Fe,
        id: self.counter,
        mass: ATOMIC_MASS_FE,
        position: position_,
        velocity: velocity_,
        force: Vector3::new(0.0, 0.0, 0.0),
        acceleration: Vector3::new(0.0, 0.0, 0.0),
        kinetic_energy: ATOMIC_MASS_FE * velocity_.magnitude_squared() / 2.0,
        potential_energy: 0.0,
      }
    };
    self.counter = self.counter + 1;

    result
  }
}

pub struct SafeAtomFactory {
  inner: Mutex<AtomFactory>
}

impl SafeAtomFactory {
  pub fn new() -> Self {
    Self {
      inner: Mutex::new(AtomFactory::new()),
    }
  }

  pub fn get_atom(&self, atom: AtomType, position_: Vector3<f64>, velocity_: Vector3<f64>) -> Atom {
    let mut factory = self.inner.lock().unwrap();
    factory.get_atom(atom, position_, velocity_)
  }

  pub fn get_atom_random(&self, atom: AtomType, lower_bound: f64, upper_bound: f64) -> Atom {
    let mut factory = self.inner.lock().unwrap();
    let mut rng = rand::thread_rng();
    let range = Uniform::new(lower_bound, upper_bound);
    let range_vel = Uniform::new(-1., 1.);

    let position_ = Vector3::new(
      range.sample(&mut rng),
      range.sample(&mut rng),
      range.sample(&mut rng),
    );
    let velocity_ = Vector3::new(
      range_vel.sample(&mut rng),
      range_vel.sample(&mut rng),
      range_vel.sample(&mut rng),
    );
    factory.get_atom(atom, position_, velocity_)
  }
}