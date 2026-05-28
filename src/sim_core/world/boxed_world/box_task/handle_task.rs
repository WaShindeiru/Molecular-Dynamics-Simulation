use nalgebra::Vector3;
use std::collections::HashMap;
use std::sync::Arc;

use crate::particle::Particle;
use crate::sim_core::world::boundary_constraint::{Compliance, EdgeCondition, ParticleCompliance};
use crate::sim_core::world::boxed_world::box_container::BoxContainer;
use crate::sim_core::world::boxed_world::box_container::sim_box::SimulationBox;
use crate::sim_core::world::boxed_world::box_task::force_task_box_container::{
  ForceTaskBoxContainer, get_needed_box_id_periodic,
};
use crate::sim_core::world::boxed_world::box_task::handle_task::handle_partial_velocity_step::{apply_velocity_constraint, handle_partial_velocity_step};
use crate::sim_core::world::boxed_world::box_task::{
  ForceTaskParticleData, ForceTaskResult, VelocityTaskParticleData, VelocityTaskResult,
};
use crate::sim_core::world::boxed_world::integration::verlet_nose_hoover::computation::{HalfVelocityResult, verlet_noose_hoover_half_velocity_position};
use crate::sim_core::world::computation::{ForceComputationOperations, compute_forces_potential};
use crate::sim_core::world::boxed_world::integration_cache::IntegrationCache;

mod partial_velocity_step;
mod handle_partial_velocity_step;

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
  // Separate CustomVelocityAtom particles from normal particles.
  let (normal_particles, custom_vel_particles): (Vec<_>, Vec<_>) = history
    .particles_of_boxes(box_ids)
    .partition(|p| !p.is_custom_velocity_atom());

  let atom_count = normal_particles.len();

  let half_velocity_response = verlet_noose_hoover_half_velocity_position(
    normal_particles,
    time_step,
    previous_thermostat_epsilon,
    atom_count,
    current_iteration,
    &container_size,
    edge_condition,
  );

  let HalfVelocityResult { half_velocity, new_position, thermostat_work, compliance } =
    half_velocity_response;

  let all_normal_particles: HashMap<usize, VelocityTaskParticleData> = compliance
    .into_iter()
    .map(|(id, compliance)| (id, VelocityTaskParticleData {
      half_velocity: half_velocity[&id],
      new_position: new_position[&id],
      thermostat_work: thermostat_work[&id],
      compliance,
    }))
    .collect();

  let (compliant, non_compliant): (HashMap<_, _>, HashMap<_, _>) = match edge_condition {
    EdgeCondition::Simple { .. } | EdgeCondition::Periodic { .. } => {
      all_normal_particles.into_iter().partition(|(_, d)| d.compliance.compliant)
    },
    EdgeCondition::PeriodicAll => (all_normal_particles, HashMap::new()),
  };

  let mut particles = apply_velocity_constraint(compliant, edge_condition);

  if !non_compliant.is_empty() {
    let corrected = handle_partial_velocity_step(
      history,
      non_compliant,
      container_size,
      edge_condition,
      previous_thermostat_epsilon,
      time_step,
    );
    particles.extend(corrected);
  }

  // CustomVelocityAtom: velocity is prescribed. The current velocity (set in the previous
  // iteration by apply_custom_velocities) is read directly from the history particle.
  // Position advances with simple Euler: pos += vel * dt.  No Verlet half-step, no
  // boundary checks (ignore_edge_conditions is always true for now).
  for particle_arc in custom_vel_particles {
    let id = particle_arc.get_id();
    let vel = *particle_arc.get_velocity();
    let new_position = particle_arc.get_position() + vel * time_step;
    particles.insert(id, VelocityTaskParticleData {
      half_velocity: vel,
      new_position,
      thermostat_work: 0.0,
      compliance: ParticleCompliance {
        compliant: true,
        x: Compliance::Compliant,
        y: Compliance::Compliant,
        z: Compliance::Compliant,
      },
    });
  }

  VelocityTaskResult { task_id, particles }
}

pub fn handle_force_batch_task(
  task_id: usize,
  boundary_condition: EdgeCondition,
  box_ids: &[usize],
  integration_cache: &IntegrationCache,
) -> ForceTaskResult {
  assert!(matches!(boundary_condition, EdgeCondition::Periodic { .. } | EdgeCondition::PeriodicAll), 
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
