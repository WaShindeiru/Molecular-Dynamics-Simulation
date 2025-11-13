use std::sync::Mutex;

use nalgebra::base::Vector3;
use crate::data::constants::{ATOMIC_MASS_C, ATOMIC_MASS_FE};
use crate::data::types::AtomType;
use crate::output::{change_length_unit, AtomDTO};
use crate::particle::atom_collection::AtomMetadata;

#[derive(Debug, PartialEq, Clone)]
pub struct Atom {
  id: u64,
  type_: AtomType,
  mass: f64,
  position: Vector3<f64>,
  velocity: Vector3<f64>,
  acceleration: Vector3<f64>,
  force: Vector3<f64>,
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
    ) -> Self {
    Atom {
      id,
      type_,
      mass,
      position,
      velocity,
      force,
      acceleration,
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
    }
  }
  
  pub fn to_atom_metadata(&self) -> AtomMetadata {
    AtomMetadata::new(
      self.id,
      &self.type_,
      self.mass,
      &self.position,
    )
  }
}

// impl PartialEq for Atom {
//     fn eq(&self, other: &Self) -> bool {
//         self.id == other.id
//     }
// }

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
      },
      AtomType::Fe => Atom{
        type_: AtomType::Fe,
        id: self.counter,
        mass: ATOMIC_MASS_FE,
        position: position_,
        velocity: velocity_,
        force: Vector3::new(0.0, 0.0, 0.0),
        acceleration: Vector3::new(0.0, 0.0, 0.0),
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

  // fn get_atom_default(&self, atom: AtomType, position_: Vector3<f64>, velocity_: Vector3<f64>) -> Atom {
  //     let mut factory = self.inner.lock().unwrap();
  //     factory.get_atom(atom, position_, velocity_)
  // }

  pub fn get_atom(&self, atom: AtomType, position_: Vector3<f64>, velocity_: Vector3<f64>) -> Atom {
    let mut factory = self.inner.lock().unwrap();
    factory.get_atom(atom, position_, velocity_)
  }
}