use std::collections::HashMap;
use nalgebra::Vector3;
use crate::data::types::AtomType;
use crate::particle::{Atom, AtomCollection, SimpleAtomContainer};
use crate::particle::atom_collection::AtomMetadata;

pub struct AtomData {
  id: u64,
  type_: AtomType,
  mass: f64,
  position: Vector3<f64>,
}

impl AtomData {
  pub fn new(id: u64, type_: AtomType, mass: f64, position: Vector3<f64>) -> Self {
    AtomData {
      id,
      type_,
      mass,
      position,
    }
  }

  pub fn get_id(&self) -> u64 {
    self.id
  }

  pub fn get_type(&self) -> &AtomType {
    &self.type_
  }

  pub fn get_mass(&self) -> f64 {
    self.mass
  }

  pub fn get_position(&self) -> &Vector3<f64> {
    &self.position
  }

  pub fn to_atom_metadata(&self) -> AtomMetadata {
    AtomMetadata::new(
      self.get_id(),
      self.get_type(),
      self.get_mass(),
      self.get_position(),
    )
  }
}

pub struct AtomDataContainer {
  atom_map: HashMap<u64, AtomData>,
}

impl AtomDataContainer {
  pub fn new() -> Self {
    AtomDataContainer {
      atom_map: HashMap::new(),
    }
  }

  pub fn add_atom(&mut self, atom: AtomData) {
    let id = atom.get_id();
    self.atom_map.insert(id, atom);
  }

  pub fn len(&self) -> usize {
    self.atom_map.len()
  }

  pub fn get_map(&self) -> &HashMap<u64, AtomData> {
    &self.atom_map
  }

  pub fn get_atom_data(&self, id: u64) -> Option<&AtomData> {
    self.atom_map.get(&id)
  }
}

impl AtomCollection for AtomDataContainer {
  fn get_atom_by_id(&self, id: u64) -> Option<AtomMetadata> {
    self.atom_map.get(&id).map(|atom| atom.to_atom_metadata())
  }

  fn get_all_atoms(&self) -> HashMap<u64, AtomMetadata> {
    let mut metadata_map: HashMap<u64, AtomMetadata> = HashMap::with_capacity(self.atom_map.len());

    for (id, atom) in &self.atom_map {
      metadata_map.insert(*id, atom.to_atom_metadata());
    }

    metadata_map
  }
}

pub struct AtomForceData {
  id: u64,
  velocity: Vector3<f64>,
  acceleration: Vector3<f64>,
  force: Vector3<f64>,
}

impl AtomForceData {
  pub fn new(id: u64, velocity: Vector3<f64>, acceleration: Vector3<f64>, force: Vector3<f64>) -> Self {
    AtomForceData {
      id,
      velocity,
      acceleration,
      force,
    }
  }

  pub fn get_id(&self) -> u64 {
    self.id
  }

  pub fn get_velocity(&self) -> &Vector3<f64> {
    &self.velocity
  }

  pub fn get_acceleration(&self) -> &Vector3<f64> {
    &self.acceleration
  }

  pub fn get_force(&self) -> &Vector3<f64> {
    &self.force
  }
}

pub struct AtomForceContainer {
  force_map: HashMap<u64, AtomForceData>,
}

impl AtomForceContainer {
  pub fn new() -> Self {
    AtomForceContainer {
      force_map: HashMap::new(),
    }
  }

  pub fn add_atom_force(&mut self, atom_force: AtomForceData) {
    let id = atom_force.get_id();
    self.force_map.insert(id, atom_force);
  }

  pub fn len(&self) -> usize {
    self.force_map.len()
  }

  pub fn get_atom_force(&self, id: u64) -> Option<&AtomForceData> {
    self.force_map.get(&id)
  }

  pub fn get_map(&self) -> &HashMap<u64, AtomForceData> {
    &self.force_map
  }
}

pub fn new_atom_container_from_parts(atom_data_container: AtomDataContainer, mut atom_force_container: AtomForceContainer) -> SimpleAtomContainer {
  let mut atom_map: HashMap<u64, Atom> = HashMap::with_capacity(atom_data_container.len());

  for (id, atom_data) in atom_data_container.atom_map {
    let atom_force = atom_force_container.force_map.remove(&id).unwrap();

    let atom = Atom::new(
      atom_data.id,
      atom_data.type_,
      atom_data.mass,
      atom_data.position,
      atom_force.velocity,
      atom_force.force,
      atom_force.acceleration,
    );

    atom_map.insert(id, atom);
  }

  SimpleAtomContainer::new_from_map(atom_map)
}
