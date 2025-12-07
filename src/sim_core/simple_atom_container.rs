use std::collections::HashMap;
use crate::particle::atom_collection::{AtomCollection, AtomMetadata};
use crate::particle::Particle;
use crate::sim_core::atom_wrapper::{AtomData, AtomDataContainer, AtomForceContainer, AtomForceData};

pub struct SimpleAtomContainer {
  atoms: Vec<Particle>,
  potential_energy: f64,
}

impl SimpleAtomContainer {
  pub fn new() -> Self {
    SimpleAtomContainer {
      atoms: Vec::new(),
      potential_energy: 0.,
    }
  }

  pub fn new_from_atoms(atoms: Vec<Particle>) -> Self {
    SimpleAtomContainer {
      atoms,
      potential_energy: 0.
    }
  }

  pub fn new_fixed_cap(capacity: usize) -> Self {
    SimpleAtomContainer {
      atoms: Vec::with_capacity(capacity),
      potential_energy: 0.
    }
  }

  pub fn add_atom(&mut self, atom: Particle) {
    let id = atom.get_id();
    self.atoms.push(atom);
  }

  pub fn get_atoms(&self) -> &Vec<Particle> {
    &self.atoms
  }

  pub fn len(&self) -> usize {
    self.atoms.len()
  }

  pub fn get_atom(&self, id: usize) -> Option<&Particle> {
    self.atoms.get(id)
  }
  
  pub fn create_parts(&self) -> (AtomDataContainer, AtomForceContainer) {
    let mut data_container = AtomDataContainer::new();
    let mut force_container = AtomForceContainer::new();

    for atom in &self.atoms {
      data_container.add_atom(AtomData::new_from_atom(atom));
      force_container.add_atom_force(AtomForceData::new_from_atom(atom));
    }

    (data_container, force_container)
  }

  pub fn set_potential_energy(&mut self, potential_energy: f64) {
    self.potential_energy = potential_energy;
  }
  
  pub fn get_potential_energy(&self) -> f64 {
    self.potential_energy
  }
}