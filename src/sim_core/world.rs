use std::collections::HashMap;
use nalgebra::Vector3;
use crate::output::{AtomDTO, WorldDTO};
use crate::particle::{compute_force_i, compute_potential_energy_i, Atom, SimpleAtomContainer};
use crate::particle::atom::ParticleOperations;
use crate::particle::atom_collection::AtomMetadata;
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

  pub fn new_from_atoms(atoms: Vec<Box<dyn ParticleOperations>>, size: Vector3<f64>) -> Self {
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

    if self.current_iteration == 0 {
      let (atom_data_container, mut force_data_container) = self.atoms.get(0).unwrap().create_parts();

      for (id, atom_data_i) in atom_data_container.get_map().iter() {
        assert_eq!(*id, atom_data_i.get_id());
        let atom_force_i = compute_force_i(&atom_data_container, atom_data_i.as_ref());
        let atom_potential_energy_i = compute_potential_energy_i(&atom_data_container, atom_data_i.as_ref());

        let atom_acceleration_i: Vector3<f64> = atom_force_i / atom_data_i.get_mass();

        let atom_force_data = AtomForceData::new(
          atom_data_i.get_id(),
          atom_force_i,
          *force_data_container.get_atom_force(*id).unwrap().get_velocity(),
          atom_acceleration_i,
          atom_potential_energy_i,
        );

        force_data_container.add_atom_force(atom_force_data);
      }

      let previous_atom_container = self.atoms.get(self.current_iteration).unwrap();
      let initial_state = new_atom_container_from_parts(atom_data_container, 
                                                        force_data_container, previous_atom_container);
      *self.atoms.get_mut(0).unwrap() = initial_state;
    }

    let previous_atom_container = self.atoms.get(self.current_iteration).unwrap();

    let mut half_velocity_cache: HashMap<u64, Vector3<f64>> = HashMap::with_capacity(self.atom_count);

    let mut atom_data_container: AtomDataContainer = AtomDataContainer::new();
    let mut force_container: AtomForceContainer = AtomForceContainer::new();

    for (_, atom_i) in previous_atom_container.get_map() {
      let half_velocity_i: Vector3<f64> = atom_i.get_velocity() + atom_i.get_acceleration() * (time_step / 2.0);
      half_velocity_cache.insert(atom_i.get_id(), half_velocity_i);

      let next_position: Vector3<f64> = atom_i.get_position() + half_velocity_i * time_step;
      let next_atom_data = AtomData::new(
        atom_i.get_id(),
        *atom_i.get_type(),
        atom_i.get_mass(),
        next_position,
      );

      atom_data_container.add_atom(next_atom_data);
    }

    for (id, atom_data_i) in atom_data_container.get_map().iter() {
      assert_eq!(*id, atom_data_i.get_id());
      let atom_force_i = compute_force_i(&atom_data_container, atom_data_i.as_ref());
      let atom_potential_energy_i = compute_potential_energy_i(&atom_data_container, atom_data_i.as_ref());

      let atom_acceleration_i: Vector3<f64> = atom_force_i / atom_data_i.get_mass();
      let atom_velocity_i: Vector3<f64> = half_velocity_cache.get(&atom_data_i.get_id()).unwrap() + atom_acceleration_i * (time_step / 2.0);

      let atom_force_data = AtomForceData::new(
        atom_data_i.get_id(),
        atom_force_i,
        atom_velocity_i,
        atom_acceleration_i,
        atom_potential_energy_i,
      );

      force_container.add_atom_force(atom_force_data);
    }

    let next_iteration_atom_container = new_atom_container_from_parts(atom_data_container, 
                                                                      force_container, previous_atom_container);

    self.current_iteration += 1;
    assert_eq!(self.current_iteration, next_iteration);

    self.atoms.push(next_iteration_atom_container);

    assert_eq!(self.atoms.len() - 1, self.current_iteration);
  }

  pub fn update_simple_vel_verlet(&mut self, time_step: f64, next_iteration: usize) {
    assert_eq!(self.atoms.len() - 1, self.current_iteration);

    let previous_atom_container = self.atoms.get(self.current_iteration).unwrap();

    let mut atom_data_container: AtomDataContainer = AtomDataContainer::new();
    let mut force_container: AtomForceContainer = AtomForceContainer::new();

    for (_, atom_i) in previous_atom_container.get_map() {
      let next_position: Vector3<f64> = atom_i.get_position() + atom_i.get_velocity() * time_step + atom_i.get_acceleration() * (0.5 * time_step * time_step);
      let next_atom_data = AtomData::new(
        atom_i.get_id(),
        *atom_i.get_type(),
        atom_i.get_mass(),
        next_position,
      );

      atom_data_container.add_atom(next_atom_data);
    }

    for (id, atom_data_i) in atom_data_container.get_map().iter() {
      assert_eq!(*id, atom_data_i.get_id());
      let previous_atom_i = previous_atom_container.get_atom_by_id(*id).unwrap();
      let atom_force_i = compute_force_i(&atom_data_container, atom_data_i.as_ref());
      let atom_potential_energy_i = compute_potential_energy_i(&atom_data_container, atom_data_i.as_ref());

      let atom_acceleration_i: Vector3<f64> = atom_force_i / atom_data_i.get_mass();
      let atom_velocity_i: Vector3<f64> = previous_atom_i.get_velocity() + (previous_atom_i.get_acceleration() + atom_acceleration_i) * (0.5 * time_step);

      let atom_force_data = AtomForceData::new(
        atom_data_i.get_id(),
        atom_force_i,
        atom_velocity_i,
        atom_acceleration_i,
        atom_potential_energy_i,
      );

      force_container.add_atom_force(atom_force_data);
    }

    let next_iteration_atom_container = new_atom_container_from_parts(atom_data_container, 
                                                                      force_container, previous_atom_container);

    self.current_iteration += 1;
    assert_eq!(self.current_iteration, next_iteration);

    self.atoms.push(next_iteration_atom_container);

    assert_eq!(self.atoms.len() - 1, self.current_iteration);
  }

  pub fn update_semi_implicit_euler(&mut self, time_step: f64, next_iteration: usize) {
    let mut next_iteration_atom_container = SimpleAtomContainer::new_fixed_cap(self.atom_count);

    assert_eq!(self.atoms.len() - 1, self.current_iteration);
    let previous_atom_container = self.atoms.get(self.current_iteration).unwrap();

    let (previous_atom_data_container, mut previous_force_data_container) = previous_atom_container.create_parts();

    for (_, atom_data_i) in previous_atom_data_container.get_map().iter() {
      let atom_force_i = compute_force_i(&previous_atom_data_container, atom_data_i.as_ref());
      let atom_potential_energy_i = compute_potential_energy_i(&previous_atom_data_container, atom_data_i.as_ref());

      let atom_force_data_i = previous_force_data_container.get_atom_force_mut(atom_data_i.get_id()).unwrap();
      atom_force_data_i.set_force(atom_force_i);
      atom_force_data_i.set_potential_energy(atom_potential_energy_i);
      let atom_acceleration_i: Vector3<f64> = atom_force_i / atom_data_i.get_mass();
      atom_force_data_i.set_acceleration(atom_acceleration_i);
    }

    for (id, atom_data_i) in previous_atom_data_container.get_map().iter() {
      assert_eq!(*id, atom_data_i.get_id());
      let previous_atom_i = previous_atom_container.get_atom_by_id(atom_data_i.get_id()).unwrap();
      let atom_force_data_i = previous_force_data_container.get_atom_force(atom_data_i.get_id()).unwrap();

      let new_velocity: Vector3<f64> = atom_force_data_i.get_velocity() + atom_force_data_i.get_acceleration() * time_step;
      let new_position: Vector3<f64> = atom_data_i.get_position() + new_velocity * time_step;
      
      let mut new_atom = previous_atom_i.custom_clone();
      new_atom.update_position(new_position);
      new_atom.set_velocity(new_velocity);
      new_atom.set_force(Vector3::new(0., 0., 0.));
      new_atom.set_acceleration(Vector3::new(0., 0., 0.));
      new_atom.set_potential_energy(0.);

      next_iteration_atom_container.add_atom(new_atom);
    }

    *self.atoms.get_mut(self.current_iteration).unwrap() 
      = new_atom_container_from_parts(previous_atom_data_container, previous_force_data_container, 
                                      previous_atom_container);

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
