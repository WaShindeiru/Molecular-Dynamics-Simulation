use nalgebra::Vector3;
use crate::particle::{potential, Particle, SimpleAtomContainer};
use super::World;

impl World {
  
  // TODO: Fix boundary conditions check, what's wrong about them?
  pub fn update_verlet(&mut self, time_step: f64, next_iteration: usize) {
    let mut next_iteration_atom_container = SimpleAtomContainer::new_fixed_cap(self.atom_count);

    assert_eq!(self.atoms.len() - 1, self.current_iteration);
    let previous_atom_container = self.atoms.get(self.current_iteration).unwrap();

    let mut half_velocity_cache: Vec<Vector3<f64>> = vec![Vector3::new(0., 0., 0.); self.atom_count];
    let mut new_position_atoms: Vec<Particle> = Vec::with_capacity(self.atom_count);

    for (i, atom_i) in previous_atom_container.get_atoms().iter().enumerate() {
      assert_eq!(i, atom_i.get_id() as usize);
      let half_velocity_i: Vector3<f64> = atom_i.get_velocity() + atom_i.get_acceleration() * (time_step / 2.0);
      half_velocity_cache[i] = half_velocity_i;

      let next_position: Vector3<f64> = atom_i.get_position() + half_velocity_i * time_step;
      let mut new_atom_data = atom_i.custom_clone();
      new_atom_data.update_position(next_position);

      // new_atom_data = self.apply_boundary_constraint(new_atom_data);

      new_position_atoms.push(new_atom_data);
    }

    let fpinfo = potential::compute_forces_potential(&new_position_atoms);
    let potential_energy = fpinfo.potential_energy;
    let forces = fpinfo.fp;

    for (i, particle_i) in new_position_atoms.iter().enumerate() {
      assert_eq!(i, particle_i.get_id() as usize);

      let new_force = forces.get(i).unwrap().force;
      let new_potential_energy = forces.get(i).unwrap().potential_energy;
      let new_acceleration = new_force / particle_i.get_mass();
      let new_velocity = half_velocity_cache.get(i).unwrap() + 0.5 * new_acceleration * time_step;

      let mut new_atom = particle_i.custom_clone();
      new_atom.set_force(new_force);
      new_atom.set_potential_energy(new_potential_energy);
      new_atom.set_acceleration(new_acceleration);
      new_atom.set_velocity(new_velocity);

      new_atom = self.apply_boundary_constraint(new_atom);

      next_iteration_atom_container.add_atom(new_atom);
    }

    next_iteration_atom_container.set_potential_energy(potential_energy);

    self.current_iteration += 1;
    assert_eq!(self.current_iteration, next_iteration);

    self.atoms.push(next_iteration_atom_container);

    assert_eq!(self.atoms.len() - 1, self.current_iteration);
  }
}