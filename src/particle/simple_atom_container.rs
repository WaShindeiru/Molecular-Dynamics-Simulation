use std::collections::HashMap;
use crate::particle::atom_collection::AtomCollection;
use super::atom::*;

pub struct SimpleAtomContainer {
  atoms: Vec<Atom>,
  atom_map: HashMap<u64, usize>,
}

impl SimpleAtomContainer {
  pub fn new() -> Self {
    SimpleAtomContainer {
      atoms: Vec::new(),
      atom_map: HashMap::new(),
    }
  }

  pub fn new_from_atoms(atoms: Vec<Atom>) -> Self {
    let mut atom_map = HashMap::with_capacity(atoms.len());
    for (index, atom) in atoms.iter().enumerate() {
      atom_map.insert(atom.get_id(), index);
    }

    SimpleAtomContainer {
      atoms,
      atom_map,
    }
  }

  pub fn new_fixed_cap(capacity: usize) -> Self {
    SimpleAtomContainer {
      atoms: Vec::with_capacity(capacity),
      atom_map: HashMap::with_capacity(capacity),
    }
  }

  pub fn add_atom(&mut self, atom: Atom) {
    let id = atom.get_id();
    self.atoms.push(atom);
    self.atom_map.insert(id, self.atoms.len() - 1);
  }

  pub fn len(&self) -> usize {
    self.atoms.len()
  }

  pub fn get_all_atoms_mutable(&mut self) -> &mut Vec<Atom> {
    &mut self.atoms
  }

  pub fn get_atom_by_index(&self, index: usize) -> Option<&Atom> { self.atoms.get(index) }
}

impl AtomCollection for SimpleAtomContainer {
  fn get_atom_by_id(&self, id: u64) -> Option<&Atom> {
    self.atom_map.get(&id).and_then(|&index| self.atoms.get(index))
  }

  fn get_all_atoms(&self) -> &Vec<Atom> {
    &self.atoms
  }
}