use std::collections::HashMap;
use nalgebra::Vector3;
use crate::particle::atom::AtomForce;
use crate::particle::{AtomCollection, SimpleAtomContainer};

pub struct SimpleForceContainer {
  forces: Vec<AtomForce>,
  atom_map: HashMap<u64, usize>,
}

impl SimpleForceContainer {
  pub fn new() -> Self {
    SimpleForceContainer {
      forces: Vec::new(),
      atom_map: HashMap::new(),
    }
  }
  
  pub fn new_from_atom_container(atom_cont: &SimpleAtomContainer) -> Self {
    let mut atom_map: HashMap<u64, usize> = HashMap::with_capacity(atom_cont.len());
    let mut forces: Vec<AtomForce> = Vec::with_capacity(atom_cont.len());
    
    for atom in atom_cont.get_all_atoms().iter() {
      let atom_force = AtomForce::from_atom(atom);
      let id = atom_force.get_id();
      forces.push(atom_force);
      atom_map.insert(id, forces.len() - 1);
    }
    
    SimpleForceContainer{
      atom_map,
      forces,
    }
  }

  pub fn set_force_for_atom(&mut self, atom_id: u64, force: Vector3<f64>) {
    let atom_index = *self.atom_map.get(&atom_id).unwrap();
    let old_atom_force = self.forces.get_mut(atom_index).unwrap();

    let acceleration = force / old_atom_force.get_mass();

    old_atom_force.set_force(force);
    old_atom_force.set_acceleration(acceleration);
  }

  pub fn get_atom_force_by_id(&self, atom_id: u64) -> &AtomForce {
    let atom_index = *self.atom_map.get(&atom_id).unwrap();
    self.forces.get(atom_index).unwrap()
  }
  
  pub fn get_atom_force_by_index(&self, index: usize) -> Option<&AtomForce> {
    self.forces.get(index)
  }
}