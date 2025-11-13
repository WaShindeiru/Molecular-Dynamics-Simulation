use std::collections::HashMap;
use nalgebra::Vector3;
use crate::output::{AtomDTO, WorldDTO};
use crate::particle::{compute_force_i, Atom, SimpleAtomContainer};
use crate::sim_core::atom_wrapper::{new_atom_container_from_parts, AtomData, AtomDataContainer, AtomForceContainer, AtomForceData};

pub struct World {
  atoms: Vec<SimpleAtomContainer>,
  atom_count: usize,
  size: Vector3<f64>,
  current_iteration: usize,
}

impl World {
  pub fn new(size: Vector3<f64>) -> Self {
    World {
      atoms: Vec::new(),
      atom_count: 0,
      size,
      current_iteration: 0,
    }
  }

  pub fn new_from_atoms(atoms: Vec<Atom>, size: Vector3<f64>) -> Self {
    let atom_count = atoms.len();
    let atom_container = SimpleAtomContainer::new_from_atoms(atoms);
    World {
      atoms: vec![atom_container],
      atom_count,
      size,
      current_iteration: 0
    }
  }

  pub fn update(&mut self, time_step: f64, next_iteration: usize) {
    assert_eq!(self.atoms.len() - 1, self.current_iteration);
    let previous_atom_container = self.atoms.get(self.current_iteration).unwrap();

    let mut half_velocity_cache: HashMap<u64, Vector3<f64>> = HashMap::with_capacity(self.atom_count);

    let mut atom_data_container: AtomDataContainer = AtomDataContainer::new();
    let mut force_container: AtomForceContainer = AtomForceContainer::new();

    for (i_index, atom_i) in previous_atom_container.get_map() {
      let half_velocity_i = atom_i.get_velocity() + atom_i.get_acceleration() * (time_step * 0.5);
      half_velocity_cache.insert(*i_index, half_velocity_i);

      let next_position = atom_i.get_position() + half_velocity_i * time_step;
      let next_atom_data = AtomData::new(
        atom_i.get_id(),
        atom_i.get_type().clone(),
        atom_i.get_mass(),
        next_position,
      );

      atom_data_container.add_atom(next_atom_data);
    }

    for (i_index, atom_data_i) in atom_data_container.get_map().iter() {
      let atom_force_i = compute_force_i(&atom_data_container, &atom_data_i.to_atom_metadata());

      let atom_acceleration_i: Vector3<f64> = atom_force_i / atom_data_i.get_mass();
      let atom_velocity_i: Vector3<f64> = half_velocity_cache.get(i_index).unwrap() + atom_acceleration_i * (time_step * 0.5);

      let atom_force_data = AtomForceData::new(
        *i_index,
        atom_force_i,
        atom_velocity_i,
        atom_acceleration_i,
      );

      force_container.add_atom_force(atom_force_data);
    }

    let next_iteration_atom_container = new_atom_container_from_parts(atom_data_container, force_container);

    self.current_iteration += 1;
    assert_eq!(self.current_iteration, next_iteration);

    self.atoms.push(next_iteration_atom_container);

    assert_eq!(self.atoms.len() - 1, self.current_iteration);
  }

  pub fn to_transfer_struct(&self) -> WorldDTO {
    let mut all_atoms_dto: Vec<Vec<AtomDTO>> = Vec::with_capacity(self.atoms.len());

    for iteration in 0..self.atoms.len() {
      let atom_container = self.atoms.get(iteration).unwrap();
      let mut atoms_dto: Vec<AtomDTO> = Vec::with_capacity(atom_container.len());

      for (id, atom) in atom_container.get_map().iter() {
        atoms_dto.push(atom.to_transfer_struct());
      }

      all_atoms_dto.push(atoms_dto);
    }
    
    WorldDTO {
      num_of_atoms: all_atoms_dto.len(),
      atoms: all_atoms_dto,
      box_x: self.size.x,
      box_y: self.size.y,
      box_z: self.size.z,
    }
  }
}
