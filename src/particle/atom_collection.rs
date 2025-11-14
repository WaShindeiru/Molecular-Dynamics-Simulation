use std::collections::HashMap;
use nalgebra::Vector3;
use crate::data::types::AtomType;

pub trait AtomMetadata {
  fn get_id(&self) -> u64;
  fn get_type(&self) -> &AtomType;
  fn get_mass(&self) -> f64;
  fn get_position(&self) -> &Vector3<f64>;
}

pub trait AtomCollection {
  fn get_atom_by_id(&self, id: u64) -> Option<&dyn AtomMetadata>;
  fn get_all_atoms(&self) -> HashMap<u64, Box<dyn AtomMetadata>>;
}