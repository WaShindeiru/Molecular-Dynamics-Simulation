use crate::particle::{potential, SimpleAtomContainer};
use crate::sim_core::world::simple_world::SimpleWorld;

impl SimpleWorld {

  pub fn update_semi_implicit_euler(&mut self, time_step: f64, next_iteration: usize) {
    let mut next_iteration_atom_container = SimpleAtomContainer::new_fixed_cap(self.atom_count);

    assert_eq!(self.atoms.len() - 1, self.current_index);
    let previous_atom_container = self.atoms.get(self.current_index).unwrap();

    let fpinfo = potential::compute_forces_potential(&previous_atom_container.get_atoms());
    let potential_energy = fpinfo.potential_energy;
    let forces = fpinfo.fp;

    for (i, particle_i) in previous_atom_container.get_atoms().iter().enumerate() {
      assert_eq!(i, particle_i.get_id() as usize);

      let new_force = forces.get(i).unwrap().force;
      let new_acceleration = new_force / particle_i.get_mass();
      let new_velocity = particle_i.get_velocity() + new_acceleration * time_step;
      let new_position = particle_i.get_position() + new_velocity * time_step;

      let mut new_atom = particle_i.custom_clone();
      new_atom.set_force(new_force);
      new_atom.set_acceleration(new_acceleration);
      new_atom.set_velocity(new_velocity);
      new_atom.update_position(new_position);
      new_atom.set_iteration(next_iteration);

      next_iteration_atom_container.add_atom(new_atom);
    }

    next_iteration_atom_container.set_potential_energy(potential_energy);

    self.current_iteration += 1;
    assert_eq!(self.current_iteration, next_iteration);

    self.atoms.push(next_iteration_atom_container);
  }
}