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
  
  pub fn new_from_atom(atom: &Atom) -> Self {
    AtomData {
      id: atom.get_id(),
      type_: atom.get_type().clone(),
      mass: atom.get_mass(),
      position: atom.get_position().clone(),
    }
  }
}

impl AtomMetadata for AtomData {
  fn get_id(&self) -> u64 {
    self.id
  }

  fn get_type(&self) -> &AtomType {
    &self.type_
  }

  fn get_mass(&self) -> f64 {
    self.mass
  }

  fn get_position(&self) -> &Vector3<f64> {
    &self.position
  }
}

pub struct AtomDataContainer {
  atom_map: HashMap<u64, Box<AtomData>>,
}

impl AtomDataContainer {
  pub fn new() -> Self {
    AtomDataContainer {
      atom_map: HashMap::new(),
    }
  }

  pub fn add_atom(&mut self, atom: AtomData) {
    let id = atom.get_id();
    self.atom_map.insert(id, Box::new(atom));
  }

  pub fn len(&self) -> usize {
    self.atom_map.len()
  }

  pub fn get_map(&self) -> &HashMap<u64, Box<AtomData>> {
    &self.atom_map
  }

  pub fn get_atom_data(&self, id: u64) -> Option<&AtomData> {
    self.atom_map.get(&id).map(|atom| atom.as_ref())
  }
}

impl AtomCollection for AtomDataContainer {
  fn get_atom_by_id(&self, id: u64) -> Option<&dyn AtomMetadata> {
    self.atom_map.get(&id).map(|atom| atom.as_ref() as &dyn AtomMetadata)
  }

  fn get_all_atoms(&self) -> HashMap<u64, Box<dyn AtomMetadata>> {
    let mut metadata_map: HashMap<u64, Box<dyn AtomMetadata>> = HashMap::with_capacity(self.atom_map.len());

    for (id, atom_data) in &self.atom_map {
      metadata_map.insert(*id, Box::new(AtomData {
        id: atom_data.id,
        type_: atom_data.type_.clone(),
        mass: atom_data.mass,
        position: atom_data.position.clone(),
      }) as Box<dyn AtomMetadata>);
    }

    metadata_map
  }
}

pub struct AtomForceData {
  id: u64,
  velocity: Vector3<f64>,
  acceleration: Vector3<f64>,
  force: Vector3<f64>,
  
  potential_energy: f64,
}

impl AtomForceData {
  pub fn new(id: u64, velocity: Vector3<f64>, acceleration: Vector3<f64>, force: Vector3<f64>,
             potential_energy: f64) -> Self {
    AtomForceData {
      id,
      velocity,
      acceleration,
      force,
      potential_energy,
    }
  }
  
  pub fn new_from_atom(atom: &Atom) -> Self {
    AtomForceData {
      id: atom.get_id(),
      velocity: atom.get_velocity().clone(),
      acceleration: atom.get_acceleration().clone(),
      force: atom.get_force().clone(),
      potential_energy: atom.get_potential_energy(),
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
  
  pub fn get_potential_energy(&self) -> f64 {
    self.potential_energy
  }
  
  pub fn set_force(&mut self, force: Vector3<f64>) {
    self.force = force;
  }
  
  pub fn set_acceleration(&mut self, acceleration: Vector3<f64>) {
    self.acceleration = acceleration;
  }
  
  pub fn set_potential_energy(&mut self, potential_energy: f64) {
    self.potential_energy = potential_energy;
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
  
  pub fn get_atom_force_mut(&mut self, id: u64) -> Option<&mut AtomForceData> {
    self.force_map.get_mut(&id)
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
      atom_force.potential_energy,
    );

    atom_map.insert(id, atom);
  }

  SimpleAtomContainer::new_from_map(atom_map)
}
