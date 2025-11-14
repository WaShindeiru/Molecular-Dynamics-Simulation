use std::collections::HashMap;
use crate::particle::atom_collection::{AtomCollection, AtomMetadata};
use crate::particle::atom::*;

pub struct SimpleAtomContainer {
  atom_map: HashMap<u64, Atom>,
}

impl SimpleAtomContainer {
  pub fn new() -> Self {
    SimpleAtomContainer {
      atom_map: HashMap::new(),
    }
  }

  pub fn new_from_atoms(atoms: Vec<Atom>) -> Self {
    let mut atom_map: HashMap<u64, Atom> = HashMap::with_capacity(atoms.len());

    for atom in atoms {
      let id = atom.get_id();
      atom_map.insert(id, atom);
    }

    SimpleAtomContainer {
      atom_map,
    }
  }

  pub fn new_from_map(atom_map: HashMap<u64, Atom>) -> Self {
    SimpleAtomContainer {
      atom_map,
    }
  }

  pub fn new_fixed_cap(capacity: usize) -> Self {
    SimpleAtomContainer {
      atom_map: HashMap::with_capacity(capacity),
    }
  }

  pub fn add_atom(&mut self, atom: Atom) {
    let id = atom.get_id();
    self.atom_map.insert(id, atom);
  }

  pub fn get_map(&self) -> &HashMap<u64, Atom> {
    &self.atom_map
  }

  pub fn len(&self) -> usize {
    self.atom_map.len()
  }

  pub fn get_atom_by_index(&self, index: u64) -> Option<&Atom> { self.atom_map.get(&index) }

  fn get_atom_by_id(&self, id: u64) -> Option<&Atom> {
    self.atom_map.get(&id)
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