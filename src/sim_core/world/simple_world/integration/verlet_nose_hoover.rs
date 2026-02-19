use std::collections::HashMap;
use nalgebra::Vector3;
use crate::data::units::K_B;
use crate::particle::{potential, Particle, SimpleAtomContainer};
use crate::sim_core::world::boundary_constraint::{apply_velocity_constraint, check_position_constraint, ParticleCompliance};
use crate::sim_core::world::simple_world::SimpleWorld;
use crate::utils::math::cos_from_vec;

fn compute_half_velocity_kinetic_energy(half_velocity_cache: &Vec<Vector3<f64>>, new_position_atoms: &Vec<Particle>) -> f64 {
  let mut kinetic_energy = 0.;

  for i in 0..half_velocity_cache.len() {
    let velocity = half_velocity_cache.get(i).unwrap();
    let mass = new_position_atoms.get(i).unwrap().get_mass();

    kinetic_energy += mass * velocity.magnitude().powi(2) / 2.0;
  }

  kinetic_energy
}

fn compute_new_thermostat_epsilon(thermostat_epsilon: f64, half_velocity_cache: &Vec<Vector3<f64>>,
                                  new_position_atoms: &Vec<Particle>, time_step: f64,
                                  q_effective_mass: f64, desired_temperature: f64) -> f64 {
  let mut new_thermostat_epsilon = thermostat_epsilon;

  let half_velocity_kinetic_energy = compute_half_velocity_kinetic_energy(half_velocity_cache, new_position_atoms);
  let num_of_particles = new_position_atoms.len();

  new_thermostat_epsilon += time_step / q_effective_mass
    * (half_velocity_kinetic_energy - 3. / 2. * K_B * desired_temperature * num_of_particles as f64);

  new_thermostat_epsilon
}

impl SimpleWorld {

  // TODO: Fix boundary conditions check
  pub fn update_verlet_nose_hoover(&mut self, time_step: f64, next_iteration: usize,
                                   desired_temperature: f64, q_effective_mass: f64) {
    let mut next_iteration_atom_container = SimpleAtomContainer::new_fixed_cap(self.atom_count);

    assert_eq!(self.atoms.len() - 1, self.current_index);
    let previous_atom_container = self.atoms.get(self.current_index).unwrap();
    let previous_thermostat_epsilon = *self.thermostat_epsilon.get(self.current_index).unwrap();

    let mut half_velocity_cache: Vec<Vector3<f64>> = vec![Vector3::new(0., 0., 0.);
                                                          self.atom_count];
    let mut new_position_atoms: Vec<Particle> = Vec::with_capacity(self.atom_count);
    let mut compliance_cache: HashMap<usize, ParticleCompliance> = HashMap::new();

    for (i, atom_i) in previous_atom_container.get_atoms().iter().enumerate() {
      assert_eq!(i, atom_i.get_id() as usize);
      let thermostat_difference = atom_i.get_acceleration() -
        previous_thermostat_epsilon * atom_i.get_velocity();
      let half_velocity_i: Vector3<f64> = atom_i.get_velocity() + thermostat_difference *
        (time_step / 2.0);
      half_velocity_cache[i] = half_velocity_i;

      let previous_position = atom_i.get_position();
      let next_position: Vector3<f64> = previous_position + half_velocity_i * time_step;
      let thermostat_work;

      if self.current_iteration == 0 {
        thermostat_work = 0.;
      } else {
        let thermostat_force = previous_thermostat_epsilon * atom_i.get_mass() *
          atom_i.get_velocity();
        let thermostat_path = next_position - previous_position;
        thermostat_work = thermostat_force.magnitude() * thermostat_path.magnitude()
          * cos_from_vec(&thermostat_force, &thermostat_path);
      }

      let mut new_atom_data = atom_i.custom_clone();
      new_atom_data.update_position(next_position);
      
      let (validated_position, compliance) = 
        check_position_constraint(new_atom_data.get_position().clone(), self.get_size());
      if !compliance.compliant {
        new_atom_data.update_position(validated_position);
      }
      compliance_cache.insert(new_atom_data.get_id() as usize, compliance);
      
      new_atom_data.set_thermostat_work(thermostat_work);

      new_position_atoms.push(new_atom_data);
    }

    let fpinfo = potential::compute_forces_potential(&new_position_atoms);
    let potential_energy = fpinfo.potential_energy;
    let forces = fpinfo.fp;

    let new_thermostat_epsilon =
      compute_new_thermostat_epsilon(previous_thermostat_epsilon, &half_velocity_cache,
                                     &new_position_atoms, time_step, q_effective_mass,
                                     desired_temperature);
    self.thermostat_epsilon.push(new_thermostat_epsilon);

    for (i, particle_i) in new_position_atoms.iter().enumerate() {
      assert_eq!(i, particle_i.get_id() as usize);

      let new_force = forces.get(i).unwrap().force;
      let new_potential_energy = forces.get(i).unwrap().potential_energy;
      let new_acceleration = new_force / particle_i.get_mass();

      let numerator = half_velocity_cache.get(i).unwrap() + 0.5 * new_acceleration * time_step;
      let denominator = 1.0 + 0.5 * time_step * new_thermostat_epsilon;

      let new_velocity = numerator / denominator;
      let compliance_i = compliance_cache.get(&i).unwrap();
      let validated_velocity = apply_velocity_constraint(&compliance_i, new_velocity);

      let mut new_atom = particle_i.custom_clone();
      new_atom.set_force(new_force);
      new_atom.set_potential_energy(new_potential_energy);
      new_atom.set_acceleration(new_acceleration);
      new_atom.set_velocity(validated_velocity);
      new_atom.set_iteration(next_iteration);

      new_atom = self.apply_boundary_constraint(new_atom);

      next_iteration_atom_container.add_atom(new_atom);
    }

    next_iteration_atom_container.set_potential_energy(potential_energy);

    self.current_iteration += 1;
    assert_eq!(self.current_iteration, next_iteration);

    self.atoms.push(next_iteration_atom_container);
  }
}