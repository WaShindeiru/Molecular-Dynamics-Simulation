use std::collections::HashMap;
use nalgebra::Vector3;
use crate::data::types::AtomType;

#[derive(Debug, PartialEq, Clone)]
pub struct AtomMetadata<'a> {
  pub id: u64,
  pub atom_type: & 'a AtomType,
  pub mass: f64,
  pub position: & 'a Vector3<f64>,
}

impl<'a> AtomMetadata <'a> {
  pub fn new(id: u64, atom_type: &'a AtomType, mass: f64, position: &'a Vector3<f64>) -> Self {
    AtomMetadata {
      id,
      atom_type,
      mass,
      position,
    }
  }

  pub fn get_id(&self) -> u64 {
    self.id
  }

  pub fn get_type(&self) -> &AtomType {
    &self.atom_type
  }

  pub fn get_mass(&self) -> f64 {
    self.mass
  }

  pub fn get_position(&self) -> &Vector3<f64> {
    &self.position
  }
}

pub trait AtomCollection {
  fn get_atom_by_id(&self, id: u64) -> Option<AtomMetadata>;
  fn get_all_atoms(&self) -> HashMap<u64, AtomMetadata>;
}