use std::collections::HashMap;
use std::sync::{Arc, RwLock};
use nalgebra::Vector3;
use crate::particle::Particle;
use crate::sim_core::world::boxed_world::box_container::BoxContainer;
use crate::sim_core::world::boxed_world::integration::verlet_nose_hoover::verlet_noose_hoover_half_velocity_position;

pub struct VelocityTaskResult {
  pub half_velocity_cache: HashMap<usize, Vector3<f64>>,
  pub new_position_atoms: HashMap<usize, Particle>,
}

pub fn handle_half_velocity_position_task(box_container: Arc<RwLock<BoxContainer>>, box_id: usize,
                                          time_step: f64,
                                          previous_thermostat_epsilon: f64, current_iteration: usize)
  -> VelocityTaskResult {

  let previous_particles = box_container.read().unwrap().current_atoms_of_box(box_id);
  let atom_count = previous_particles.len();
  let (half_velocity_cache, new_position_atoms) =
    verlet_noose_hoover_half_velocity_position(&previous_particles, time_step,
                                               previous_thermostat_epsilon, atom_count, current_iteration);
  VelocityTaskResult {
    half_velocity_cache,
    new_position_atoms,
  }
}