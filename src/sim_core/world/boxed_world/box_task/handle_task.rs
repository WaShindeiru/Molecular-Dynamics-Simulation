use crate::sim_core::world::boundary_constraint::EdgeCondition;
use crate::sim_core::world::boxed_world::box_container::BoxContainer;
use crate::sim_core::world::boxed_world::box_container::sim_box::SimulationBox;
use crate::sim_core::world::boxed_world::box_task::force_task_box_container::{
  ForceTaskBoxContainer, get_needed_box_id_periodic,
};
use crate::sim_core::world::boxed_world::box_task::{
  ForceTaskParticleData, ForceTaskResult, VelocityTaskParticleData, VelocityTaskResult,
};
use crate::sim_core::world::boxed_world::integration::verlet_nose_hoover::computation::{
  ForceComputationOperations, compute_forces_potential, verlet_noose_hoover_half_velocity_position,
};
use crate::sim_core::world::boxed_world::integration_cache::IntegrationCache;
use nalgebra::Vector3;
use std::collections::HashMap;
use std::sync::Arc;

pub fn handle_velocity_batch_task(
  task_id: usize,
  box_ids: &[usize],
  history: &BoxContainer<Arc<SimulationBox>>,
  time_step: f64,
  previous_thermostat_epsilon: f64,
  current_iteration: usize,
  container_size: Vector3<f64>,
  edge_condition: EdgeCondition,
) -> VelocityTaskResult {
  let atom_count = history.atom_count_of_boxes(box_ids);

  let half_velocity_response = verlet_noose_hoover_half_velocity_position(
    history.particles_of_boxes(box_ids),
    time_step,
    previous_thermostat_epsilon,
    atom_count,
    current_iteration,
    &container_size,
    edge_condition,
  );

  let particles: HashMap<usize, VelocityTaskParticleData> = half_velocity_response
    .half_velocity
    .keys()
    .map(|id_| {
      let id = *id_;
      let data = VelocityTaskParticleData {
        half_velocity: *half_velocity_response.half_velocity.get(&id).unwrap(),
        new_position: *half_velocity_response.new_position.get(&id).unwrap(),
        thermostat_work: *half_velocity_response.thermostat_work.get(&id).unwrap(),
        compliance: *half_velocity_response.compliance.get(&id).unwrap(),
      };
      (id, data)
    })
    .collect();

  VelocityTaskResult { task_id, particles }
}

pub fn handle_force_batch_task(
  task_id: usize,
  boundary_condition: EdgeCondition,
  box_ids: &[usize],
  integration_cache: &IntegrationCache,
) -> ForceTaskResult {
  assert!(boundary_condition == EdgeCondition::Periodic || boundary_condition == EdgeCondition::PeriodicAll, 
    "Only periodic and periodic_all boundary condition is supported for force batch task");

  let config = integration_cache.box_cache().config();
  let needed_ids = get_needed_box_id_periodic(&box_ids.to_vec(), config);

  let view = integration_cache.box_cache().view_select_boxes(&needed_ids);
  let force_container = ForceTaskBoxContainer::new(view, boundary_condition);

  let mut particles: HashMap<usize, ForceTaskParticleData> = HashMap::new();
  let mut potential_energy_total = 0.0f64;
  let mut optimization_considered_total = 0usize;
  let mut optimization_ignored_total = 0usize;

  for &box_id in box_ids {
    let particles_i: Vec<Box<dyn ForceComputationOperations>> = force_container
      .atoms_for_box(box_id)
      .expect("Primary box must be present in ForceTaskBoxContainer")
      .collect();

    let mut particles_j: Vec<Box<dyn ForceComputationOperations>> =
      force_container.neighbour_atoms_periodic(box_id).collect();
    particles_j.extend(force_container.atoms_for_box(box_id).unwrap());

    #[cfg(debug_assertions)]
    for particle in particles_i.iter() {
      let particle_id = particle.get_id();
      let particle_box_id = force_container.view().particle_box_id(particle_id);
      let correct_position = force_container
        .view()
        .get_box(particle_box_id)
        .expect("Logged particle box must be present in ForceTaskBoxContainer")
        .particle(particle_id)
        .get_position()
        .clone();
      let proxy_position = particle.get_position();

      log::debug!(
        "Force task {} box {} particles_i particle {} correct_position {:?} proxy_position {:?}",
        task_id,
        box_id,
        particle_id,
        correct_position,
        proxy_position
      );
    }

    #[cfg(debug_assertions)]
    for particle in particles_j.iter() {
      let particle_id = particle.get_id();
      let particle_box_id = force_container.view().particle_box_id(particle_id);
      let correct_position = force_container
        .view()
        .get_box(particle_box_id)
        .expect("Logged particle box must be present in ForceTaskBoxContainer")
        .particle(particle_id)
        .get_position()
        .clone();
      let proxy_position = particle.get_position();

      log::debug!(
        "Force task {} box {} particles_j particle {} correct_position {:?} proxy_position {:?}",
        task_id,
        box_id,
        particle_id,
        correct_position,
        proxy_position
      );
    }

    let info = compute_forces_potential(&particles_i, &particles_j);

    #[cfg(debug_assertions)]
    for particle in particles_j.iter() {
      let particle_id = particle.get_id();
      let fp = info.fp.get(&particle_id).unwrap();
      let particle_box_id = force_container.view().particle_box_id(particle_id);
      let real_position = force_container
        .view()
        .get_box(particle_box_id)
        .expect("Logged particle box must be present in ForceTaskBoxContainer")
        .particle(particle_id)
        .get_position()
        .clone();
      let proxy_position = particle.get_position();

      log::debug!(
        "Force task {} box {} computed particle {} real_position {:?} proxy_position {:?} force {:?} potential_energy {}",
        task_id,
        box_id,
        particle_id,
        real_position,
        proxy_position,
        fp.force,
        fp.potential_energy
      );
    }

    potential_energy_total += info.potential_energy;
    optimization_considered_total += info.optimization_considered;
    optimization_ignored_total += info.optimization_ignored;

    for particle in particles_j.iter() {
      let id = particle.get_id();
      let fp = info.fp.get(&id).unwrap();
      let particle_box_id = force_container.view().particle_box_id(id);
      particles
        .entry(id)
        .and_modify(|entry| {
          entry.force += fp.force;
          entry.potential_energy += fp.potential_energy;
        })
        .or_insert(ForceTaskParticleData {
          box_id: particle_box_id,
          force: fp.force,
          potential_energy: fp.potential_energy,
        });
    }
  }

  #[cfg(debug_assertions)]
  for (particle_id, particle_data) in particles.iter() {
    let position = integration_cache
      .box_cache()
      .get_box(particle_data.box_id)
      .particle(*particle_id)
      .get_position()
      .clone();

    log::debug!(
      "Force task {} result particle {} position {:?} force {:?} potential_energy {}",
      task_id,
      particle_id,
      position,
      particle_data.force,
      particle_data.potential_energy
    );
  }

  ForceTaskResult {
    task_id,
    potential_energy: potential_energy_total,
    optimization_considered: optimization_considered_total,
    optimization_ignored: optimization_ignored_total,
    particles,
  }
}
