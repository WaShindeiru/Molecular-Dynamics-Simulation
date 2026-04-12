use std::collections::HashMap;
use std::sync::{Arc, RwLock};
use nalgebra::Vector3;
use crate::sim_core::world::boxed_world::box_task::ForceTaskResult;
use crate::sim_core::world::boxed_world::box_task::{VelocityTaskResult, VelocityParticle};
use crate::particle::Particle;
use crate::sim_core::world::boundary_constraint::EdgeCondition;
use crate::sim_core::world::boxed_world::box_container::BoxContainer;
use crate::sim_core::world::boxed_world::integration::verlet_nose_hoover::computation::{compute_forces_potential, verlet_noose_hoover_half_velocity_position, FPInfoBoxed, ForceComputationOperations};

pub fn handle_half_velocity_position_task(box_container: Arc<RwLock<BoxContainer>>, box_id: usize,
                                          time_step: f64,
                                          previous_thermostat_epsilon: f64, current_iteration: usize,
                                          container_size: Vector3<f64>, edge_condition: EdgeCondition)
  -> VelocityTaskResult {

  let lock = box_container.read().unwrap();
  let previous_particles = lock.current_atoms_of_box(box_id);
  let atom_count = previous_particles.len();
  let half_velocity_response =
    verlet_noose_hoover_half_velocity_position(previous_particles.values(), time_step,
                                               previous_thermostat_epsilon, atom_count, current_iteration,
                                               &container_size, edge_condition);

  let mut particles: HashMap<usize, VelocityParticle> = HashMap::new();

  for (id, particle) in previous_particles {
    let half_velocity = *half_velocity_response.half_velocity.get(id).unwrap();
    let new_position = *half_velocity_response.new_position.get(id).unwrap();
    let thermostat_work = *half_velocity_response.thermostat_work.get(id).unwrap();
    let compliance = *half_velocity_response.compliance.get(id).unwrap();

    let mut new_particle = particle.reset_clone();
    new_particle.update_position(new_position);
    new_particle.set_thermostat_work(thermostat_work);

    particles.insert(*id, VelocityParticle {
      particle: new_particle,
      half_velocity,
      compliance,
    });
  }

  VelocityTaskResult {
    particles,
  }
}

pub fn handle_force_task(box_container: Arc<RwLock<BoxContainer>>, box_id: usize)
  -> ForceTaskResult {
  
  let info_boxed: FPInfoBoxed;
  let mut acceleration: HashMap<usize, Vector3<f64>> = HashMap::new();

  {
    let lock = box_container.read().unwrap();

    let particles_i: Vec<_> = lock.atoms_for_force_computation_of_single_integration_box(box_id).collect();
    let particles_j: Vec<_> = lock.integration_particles_of_neighbour_boxes(box_id)
      .chain(lock.atoms_for_force_computation_of_single_integration_box(box_id))
      .collect();

    info_boxed = compute_forces_potential(particles_i.iter(), particles_j.iter());

    for particle in particles_j.iter() {
      let j_id = particle.get_id();
      let new_force = info_boxed.fp.get(&j_id).unwrap().force;
      let new_acceleration: Vector3<f64> = new_force / particle.get_mass();
      assert!(!new_acceleration.x.is_nan() && !new_acceleration.y.is_nan() && !new_acceleration.z.is_nan());

      acceleration.insert(j_id, new_acceleration);
    }
  }
  
  ForceTaskResult {
    box_id,
    info_boxed,
    acceleration,
  }
}