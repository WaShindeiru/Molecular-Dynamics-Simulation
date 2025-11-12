use nalgebra::Vector3;
use crate::output::{change_length_unit, AtomDTO, WorldDTO};
use crate::sim_core::simple_force_container::SimpleForceContainer;
use crate::particle::{compute_force_i, Atom, AtomCollection, SimpleAtomContainer};
use crate::particle::atom::AtomForce;

pub struct World {
  atoms: Vec<SimpleAtomContainer>,
  atom_forces: Vec<SimpleForceContainer>,
  atom_count: usize,
  size: Vector3<f64>,
  integration_scheme: fn(&Atom, &AtomForce, f64) -> Atom,
  current_iteration: usize,
}


impl World {
  pub fn new(size: Vector3<f64>, integration_scheme: fn(&Atom, &AtomForce, f64) -> Atom) -> Self {
    World {
      atoms: Vec::new(),
      atom_forces: Vec::new(),
      atom_count: 0,
      size,
      integration_scheme,
      current_iteration: 0,
    }
  }

  pub fn new_from_atoms(atoms: Vec<Atom>, size: Vector3<f64>, integration_scheme: fn(&Atom, &AtomForce, f64) -> Atom) -> Self {
    let atom_count = atoms.len();
    let atom_container = SimpleAtomContainer::new_from_atoms(atoms);
    let force_container = SimpleForceContainer::new_from_atom_container(&atom_container);
    World {
      atoms: vec![atom_container],
      atom_forces: vec![force_container],
      atom_count,
      size,
      integration_scheme,
      current_iteration: 0
    }
  }

  pub fn update(&mut self, time_step: f64, next_iteration: usize) {
    let mut next_iteration_atom_container = SimpleAtomContainer::new_fixed_cap(self.atom_count);

    assert_eq!(self.atoms.len() - 1, self.current_iteration);
    assert_eq!(self.atom_forces.len() - 1, self.current_iteration);
    let previous_atom_container = self.atoms.get(self.current_iteration).unwrap();
    let previous_forces_container = self.atom_forces.get_mut(self.current_iteration).unwrap();

    for atom_i in previous_atom_container.get_all_atoms().iter() {
      let atom_i_id = atom_i.get_id();
      let atom_i_force = compute_force_i(previous_atom_container, atom_i);
      previous_forces_container.set_force_for_atom(atom_i_id, atom_i_force);
    }

    for i in 0..previous_atom_container.get_all_atoms().len() {
      let atom_i = previous_atom_container.get_all_atoms().get(i).unwrap();
      let atom_force_i = previous_forces_container.get_atom_force_for_atom_id(atom_i.get_id());

      let new_atom_i = (self.integration_scheme)(atom_i, atom_force_i, time_step);
      next_iteration_atom_container.add_atom(new_atom_i);
    }

    let next_iteration_force_container = SimpleForceContainer::new_from_atom_container(&next_iteration_atom_container);

    self.current_iteration += 1;
    assert_eq!(self.current_iteration, next_iteration);

    self.atoms.push(next_iteration_atom_container);
    self.atom_forces.push(next_iteration_force_container);

    assert_eq!(self.atoms.len() - 1, self.current_iteration);
    assert_eq!(self.atom_forces.len() - 1, self.current_iteration);
  }

  pub fn to_transfer_struct(&self) -> WorldDTO {
    let mut all_atoms_dto: Vec<Vec<AtomDTO>> = Vec::with_capacity(self.atoms.len());

    for iteration in 0..self.atoms.len() {
      let atom_container = self.atoms.get(iteration).unwrap();
      let mut atoms_dto: Vec<AtomDTO> = Vec::with_capacity(atom_container.len());

      for atom in atom_container.get_all_atoms().iter() {
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
