use std::collections::HashMap;
use std::sync::{Arc, RwLock};
use nalgebra::Vector3;
use crate::sim_core::world::boxed_world::box_task::ForceTaskResult;
use crate::sim_core::world::boxed_world::box_task::VelocityTaskResult;
use crate::particle::Particle;
use crate::sim_core::world::boxed_world::box_container::BoxContainer;
use crate::sim_core::world::boxed_world::box_task::BoxResult;
use crate::sim_core::world::boxed_world::integration::verlet_nose_hoover::{compute_forces_potential, verlet_noose_hoover_half_velocity_position, FPInfoBoxed};

pub fn handle_half_velocity_position_task(box_container: Arc<RwLock<BoxContainer>>, box_id: usize,
                                          time_step: f64,
                                          previous_thermostat_epsilon: f64, current_iteration: usize)
  -> VelocityTaskResult {

  let lock = box_container.read().unwrap();
  let previous_particles = lock.current_atoms_of_box(box_id);
  let atom_count = previous_particles.len();
  let (half_velocity_cache, new_position_atoms) =
    verlet_noose_hoover_half_velocity_position(previous_particles.values(), time_step,
                                               previous_thermostat_epsilon, atom_count, current_iteration);
  VelocityTaskResult {
    half_velocity_cache,
    new_position_atoms,
  }
}

pub fn handle_force_task(box_container: Arc<RwLock<BoxContainer>>, box_id: usize)
  -> ForceTaskResult {
  
  let particles_i: Vec<&Particle>;
  let particles_j: Vec<&Particle>;
  let info_boxed: FPInfoBoxed;
  let mut acceleration: HashMap<usize, Vector3<f64>> = HashMap::new();

  {
    let lock = box_container.read().unwrap();
    particles_i = lock.atoms_of_given_integration_box(box_id).collect();
    particles_j = lock.integration_particles_of_neighbour_boxes(box_id)
      .chain(particles_i.iter().copied())
      .collect();

    info_boxed = compute_forces_potential(particles_i.iter(), particles_j.iter());

    for particle_j_ in particles_j.iter() {
      let particle_j = *particle_j_;
      let j_id = particle_j.get_id() as usize;
      let new_force = info_boxed.fp.get(&j_id).unwrap().force;
      let new_acceleration: Vector3<f64> = new_force / particle_j.get_mass();

      acceleration.insert(j_id, new_acceleration);
    }
  }
  
  ForceTaskResult {
    box_id,
    info_boxed,
    acceleration,
  }
}