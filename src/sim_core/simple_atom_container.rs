use std::collections::HashMap;
use crate::particle::atom_collection::{AtomCollection, AtomMetadata};
use crate::particle::atom::*;
use crate::sim_core::atom_wrapper::{AtomData, AtomDataContainer, AtomForceContainer, AtomForceData};

pub struct SimpleAtomContainer {
  atom_map: HashMap<u64, Box<dyn ParticleOperations>>,
}

impl SimpleAtomContainer {
  pub fn new() -> Self {
    SimpleAtomContainer {
      atom_map: HashMap::new(),
    }
  }

  pub fn new_from_atoms(atoms: Vec<Box<dyn ParticleOperations>>) -> Self {
    let mut atom_map: HashMap<u64, Box<dyn ParticleOperations>> = HashMap::with_capacity(atoms.len());

    for atom in atoms {
      let id = atom.get_id();
      atom_map.insert(id, atom);
    }

    SimpleAtomContainer {
      atom_map,
    }
  }

  pub fn new_from_map(atom_map: HashMap<u64, Box<dyn ParticleOperations>>) -> Self {
    SimpleAtomContainer {
      atom_map,
    }
  }

  pub fn new_fixed_cap(capacity: usize) -> Self {
    SimpleAtomContainer {
      atom_map: HashMap::with_capacity(capacity),
    }
  }

  pub fn add_atom(&mut self, atom: Box<dyn ParticleOperations>) {
    let id = atom.get_id();
    self.atom_map.insert(id, atom);
  }

  pub fn get_map(&self) -> &HashMap<u64, Box<dyn ParticleOperations>> {
    &self.atom_map
  }

  pub fn len(&self) -> usize {
    self.atom_map.len()
  }

  pub fn get_atom_by_index(&self, index: u64) -> Option<&Box<dyn ParticleOperations>> { self.atom_map.get(&index) }

  pub fn get_atom_by_id(&self, id: u64) -> Option<&Box<dyn ParticleOperations>> {
    self.atom_map.get(&id)
  }
  
  pub fn create_parts(&self) -> (AtomDataContainer, AtomForceContainer) {
    let mut data_container = AtomDataContainer::new();
    let mut force_container = AtomForceContainer::new();

    for (id, atom) in &self.atom_map {
      data_container.add_atom(AtomData::new_from_atom(atom));
      force_container.add_atom_force(AtomForceData::new_from_atom(atom));
    }

    (data_container, force_container)
  }
}

// impl AtomCollection for SimpleAtomContainer {
// 
//   fn get_all_atoms(&self) -> HashMap<u64, AtomMetadata> {
//     let mut metadata_map: HashMap<u64, AtomMetadata> = HashMap::with_capacity(self.atom_map.len());
// 
//     for (id, atom) in &self.atom_map {
//       metadata_map.insert(*id, atom.to_atom_metadata());
//     }
// 
//     metadata_map
//   }
// }