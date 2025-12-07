use nalgebra::Vector3;
use crate::output::{AtomDTO, WorldDTO};
use crate::particle::{potential, SimpleAtomContainer};
use crate::particle::atom::ParticleOperations;

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

  pub fn update_verlet(&mut self, time_step: f64, next_iteration: usize) {
    let mut next_iteration_atom_container = SimpleAtomContainer::new_fixed_cap(self.atom_count);

    assert_eq!(self.atoms.len() - 1, self.current_iteration);
    let previous_atom_container = self.atoms.get(self.current_iteration).unwrap();

    let mut half_velocity_cache: Vec<Vector3<f64>> = vec![Vector3::new(0., 0., 0.); self.atom_count];
    let mut new_position_atoms: Vec<Box<dyn ParticleOperations>> = Vec::with_capacity(self.atom_count);
    
    let max_key = previous_atom_container.get_map().keys().max().unwrap();
    for i in 0..=*max_key {
      let atom_i = previous_atom_container.get_atom_by_id(i).unwrap();
      let half_velocity_i: Vector3<f64> = atom_i.get_velocity() + atom_i.get_acceleration() * (time_step / 2.0);
      half_velocity_cache[i as usize] = half_velocity_i;
      
      let next_position: Vector3<f64> = atom_i.get_position() + half_velocity_i * time_step;
      let mut new_atom_data = atom_i.custom_clone();
      new_atom_data.update_position(next_position);
      new_position_atoms.push(new_atom_data);
    }
    
    let fpinfo = potential::compute_forces_potential(&new_position_atoms);
    let potential_energy = fpinfo.potential_energy;
    let forces = fpinfo.fp;
    
    for (i, particle_i) in new_position_atoms.iter().enumerate() {
      assert_eq!(i, particle_i.get_id() as usize);
      
      let new_force = forces.get(i).unwrap().force;
      let new_acceleration = new_force / particle_i.get_mass();
      let new_velocity = half_velocity_cache.get(i).unwrap() + 0.5 * new_acceleration * time_step;
      
      let mut new_atom = particle_i.custom_clone();
      new_atom.set_force(new_force);
      new_atom.set_acceleration(new_acceleration);
      new_atom.set_velocity(new_velocity);
      
      next_iteration_atom_container.add_atom(new_atom);
    }
    
    next_iteration_atom_container.set_potential_energy(potential_energy);

    self.current_iteration += 1;
    assert_eq!(self.current_iteration, next_iteration);

    self.atoms.push(next_iteration_atom_container);

    assert_eq!(self.atoms.len() - 1, self.current_iteration);
  }

  // pub fn update_old(&mut self, time_step: f64, next_iteration: usize) {
  //   assert_eq!(self.atoms.len() - 1, self.current_iteration);
  //
  //   if self.current_iteration == 0 {
  //     let (atom_data_container, mut force_data_container) = self.atoms.get(0).unwrap().create_parts();
  //
  //     for (id, atom_data_i) in atom_data_container.get_map().iter() {
  //       assert_eq!(*id, atom_data_i.get_id());
  //       let atom_force_i = compute_force_i(&atom_data_container, atom_data_i.as_ref());
  //       let atom_potential_energy_i = compute_potential_energy_i(&atom_data_container, atom_data_i.as_ref());
  //
  //       let atom_acceleration_i: Vector3<f64> = atom_force_i / atom_data_i.get_mass();
  //
  //       let atom_force_data = AtomForceData::new(
  //         atom_data_i.get_id(),
  //         atom_force_i,
  //         *force_data_container.get_atom_force(*id).unwrap().get_velocity(),
  //         atom_acceleration_i,
  //         atom_potential_energy_i,
  //       );
  //
  //       force_data_container.add_atom_force(atom_force_data);
  //     }
  //
  //     let previous_atom_container = self.atoms.get(self.current_iteration).unwrap();
  //     let initial_state = new_atom_container_from_parts(atom_data_container,
  //                                                       force_data_container, previous_atom_container);
  //     *self.atoms.get_mut(0).unwrap() = initial_state;
  //   }
  //
  //   let previous_atom_container = self.atoms.get(self.current_iteration).unwrap();
  //
  //   let mut half_velocity_cache: HashMap<u64, Vector3<f64>> = HashMap::with_capacity(self.atom_count);
  //
  //   let mut atom_data_container: AtomDataContainer = AtomDataContainer::new();
  //   let mut force_container: AtomForceContainer = AtomForceContainer::new();
  //
  //   for (_, atom_i) in previous_atom_container.get_map() {
  //     let half_velocity_i: Vector3<f64> = atom_i.get_velocity() + atom_i.get_acceleration() * (time_step / 2.0);
  //     half_velocity_cache.insert(atom_i.get_id(), half_velocity_i);
  //
  //     let next_position: Vector3<f64> = atom_i.get_position() + half_velocity_i * time_step;
  //     let next_atom_data = AtomData::new(
  //       atom_i.get_id(),
  //       *atom_i.get_type(),
  //       atom_i.get_mass(),
  //       next_position,
  //     );
  //
  //     atom_data_container.add_atom(next_atom_data);
  //   }
  //
  //   for (id, atom_data_i) in atom_data_container.get_map().iter() {
  //     assert_eq!(*id, atom_data_i.get_id());
  //     let atom_force_i = compute_force_i(&atom_data_container, atom_data_i.as_ref());
  //     let atom_potential_energy_i = compute_potential_energy_i(&atom_data_container, atom_data_i.as_ref());
  //
  //     let atom_acceleration_i: Vector3<f64> = atom_force_i / atom_data_i.get_mass();
  //     let atom_velocity_i: Vector3<f64> = half_velocity_cache.get(&atom_data_i.get_id()).unwrap() + atom_acceleration_i * (time_step / 2.0);
  //
  //     let atom_force_data = AtomForceData::new(
  //       atom_data_i.get_id(),
  //       atom_force_i,
  //       atom_velocity_i,
  //       atom_acceleration_i,
  //       atom_potential_energy_i,
  //     );
  //
  //     force_container.add_atom_force(atom_force_data);
  //   }
  //
  //   let next_iteration_atom_container = new_atom_container_from_parts(atom_data_container,
  //                                                                     force_container, previous_atom_container);
  //
  //   self.current_iteration += 1;
  //   assert_eq!(self.current_iteration, next_iteration);
  //
  //   self.atoms.push(next_iteration_atom_container);
  //
  //   assert_eq!(self.atoms.len() - 1, self.current_iteration);
  // }
  //
  // pub fn update_simple_vel_verlet(&mut self, time_step: f64, next_iteration: usize) {
  //   assert_eq!(self.atoms.len() - 1, self.current_iteration);
  //
  //   let previous_atom_container = self.atoms.get(self.current_iteration).unwrap();
  //
  //   let mut atom_data_container: AtomDataContainer = AtomDataContainer::new();
  //   let mut force_container: AtomForceContainer = AtomForceContainer::new();
  //
  //   for (_, atom_i) in previous_atom_container.get_map() {
  //     let next_position: Vector3<f64> = atom_i.get_position() + atom_i.get_velocity() * time_step + atom_i.get_acceleration() * (0.5 * time_step * time_step);
  //     let next_atom_data = AtomData::new(
  //       atom_i.get_id(),
  //       *atom_i.get_type(),
  //       atom_i.get_mass(),
  //       next_position,
  //     );
  //
  //     atom_data_container.add_atom(next_atom_data);
  //   }
  //
  //   for (id, atom_data_i) in atom_data_container.get_map().iter() {
  //     assert_eq!(*id, atom_data_i.get_id());
  //     let previous_atom_i = previous_atom_container.get_atom_by_id(*id).unwrap();
  //     let atom_force_i = compute_force_i(&atom_data_container, atom_data_i.as_ref());
  //     let atom_potential_energy_i = compute_potential_energy_i(&atom_data_container, atom_data_i.as_ref());
  //
  //     let atom_acceleration_i: Vector3<f64> = atom_force_i / atom_data_i.get_mass();
  //     let atom_velocity_i: Vector3<f64> = previous_atom_i.get_velocity() + (previous_atom_i.get_acceleration() + atom_acceleration_i) * (0.5 * time_step);
  //
  //     let atom_force_data = AtomForceData::new(
  //       atom_data_i.get_id(),
  //       atom_force_i,
  //       atom_velocity_i,
  //       atom_acceleration_i,
  //       atom_potential_energy_i,
  //     );
  //
  //     force_container.add_atom_force(atom_force_data);
  //   }
  //
  //   let next_iteration_atom_container = new_atom_container_from_parts(atom_data_container,
  //                                                                     force_container, previous_atom_container);
  //
  //   self.current_iteration += 1;
  //   assert_eq!(self.current_iteration, next_iteration);
  //
  //   self.atoms.push(next_iteration_atom_container);
  //
  //   assert_eq!(self.atoms.len() - 1, self.current_iteration);
  // }

  pub fn update_semi_implicit_euler(&mut self, time_step: f64, next_iteration: usize) {
    let mut next_iteration_atom_container = SimpleAtomContainer::new_fixed_cap(self.atom_count);

    assert_eq!(self.atoms.len() - 1, self.current_iteration);
    let previous_atom_container = self.atoms.get(self.current_iteration).unwrap();

    let mut previous_atoms: Vec<Box<dyn ParticleOperations>> = Vec::with_capacity(self.atom_count);

    let max_key = previous_atom_container.get_map().keys().max().unwrap();
    for i in 0..=*max_key {
      previous_atoms.push(previous_atom_container.get_atom_by_id(i).unwrap().custom_clone());
    }

    let fpinfo = potential::compute_forces_potential(&previous_atoms);
    let potential_energy = fpinfo.potential_energy;
    let forces = fpinfo.fp;

    for (i_, particle_i) in previous_atom_container.get_map().iter() {
      let i = *i_;
      assert_eq!(i, particle_i.get_id());

      let new_force = forces.get(i as usize).unwrap().force;
      let new_acceleration = new_force / particle_i.get_mass();
      let new_velocity = particle_i.get_velocity() + new_acceleration * time_step;
      let new_position = particle_i.get_position() + new_velocity * time_step;

      let mut new_atom = particle_i.custom_clone();
      new_atom.set_force(new_force);
      new_atom.set_acceleration(new_acceleration);
      new_atom.set_velocity(new_velocity);
      new_atom.update_position(new_position);

      next_iteration_atom_container.add_atom(new_atom);
    }

    next_iteration_atom_container.set_potential_energy(potential_energy);

    self.current_iteration += 1;
    assert_eq!(self.current_iteration, next_iteration);

    self.atoms.push(next_iteration_atom_container);

    assert_eq!(self.atoms.len() - 1, self.current_iteration);
  }

  // pub fn update_semi_implicit_euler_old(&mut self, time_step: f64, next_iteration: usize) {
  //   let mut next_iteration_atom_container = SimpleAtomContainer::new_fixed_cap(self.atom_count);
  //
  //   assert_eq!(self.atoms.len() - 1, self.current_iteration);
  //   let previous_atom_container = self.atoms.get(self.current_iteration).unwrap();
  //
  //   let (previous_atom_data_container, mut previous_force_data_container) = previous_atom_container.create_parts();
  //
  //   for (_, atom_data_i) in previous_atom_data_container.get_map().iter() {
  //     let atom_force_i = compute_force_i(&previous_atom_data_container, atom_data_i.as_ref());
  //     // let atom_potential_energy_i = compute_potential_energy_i(&previous_atom_data_container, atom_data_i.as_ref());
  //
  //     let atom_force_data_i = previous_force_data_container.get_atom_force_mut(atom_data_i.get_id()).unwrap();
  //     atom_force_data_i.set_force(atom_force_i);
  //     // atom_force_data_i.set_potential_energy(atom_potential_energy_i);
  //     let atom_acceleration_i: Vector3<f64> = atom_force_i / atom_data_i.get_mass();
  //     atom_force_data_i.set_acceleration(atom_acceleration_i);
  //   }
  //
  //   for (id, atom_data_i) in previous_atom_data_container.get_map().iter() {
  //     assert_eq!(*id, atom_data_i.get_id());
  //     let previous_atom_i = previous_atom_container.get_atom_by_id(atom_data_i.get_id()).unwrap();
  //     let atom_force_data_i = previous_force_data_container.get_atom_force(atom_data_i.get_id()).unwrap();
  //
  //     let new_velocity: Vector3<f64> = atom_force_data_i.get_velocity() + atom_force_data_i.get_acceleration() * time_step;
  //     let new_position: Vector3<f64> = atom_data_i.get_position() + new_velocity * time_step;
  //
  //     let mut new_atom = previous_atom_i.custom_clone();
  //     new_atom.update_position(new_position);
  //     new_atom.set_velocity(new_velocity);
  //     new_atom.set_force(Vector3::new(0., 0., 0.));
  //     new_atom.set_acceleration(Vector3::new(0., 0., 0.));
  //     new_atom.set_potential_energy(0.);
  //
  //     next_iteration_atom_container.add_atom(new_atom);
  //   }
  //
  //   *self.atoms.get_mut(self.current_iteration).unwrap()
  //     = new_atom_container_from_parts(previous_atom_data_container, previous_force_data_container,
  //                                     previous_atom_container);
  //
  //   self.current_iteration += 1;
  //   assert_eq!(self.current_iteration, next_iteration);
  //
  //   self.atoms.push(next_iteration_atom_container);
  //   assert_eq!(self.atoms.len() - 1, self.current_iteration);
  // }
  //
  pub fn to_transfer_struct(&self) -> WorldDTO {
    let mut all_atoms_dto: Vec<Vec<AtomDTO>> = Vec::with_capacity(self.atoms.len());
    let mut potential_energies: Vec<f64> = Vec::with_capacity(self.atoms.len());

    for iteration in 0..self.atoms.len() {
      let atom_container = self.atoms.get(iteration).unwrap();
      let mut atoms_dto: Vec<AtomDTO> = Vec::with_capacity(atom_container.len());

      for (id, atom) in atom_container.get_map().iter() {
        atoms_dto.push(atom.to_transfer_struct());
      }

      all_atoms_dto.push(atoms_dto);
      potential_energies.push(atom_container.get_potential_energy());
    }

    WorldDTO {
      num_of_atoms: all_atoms_dto.len(),
      atoms: all_atoms_dto,
      potential_energy: potential_energies,
      box_x: self.size.x,
      box_y: self.size.y,
      box_z: self.size.z,
    }
  }
}
