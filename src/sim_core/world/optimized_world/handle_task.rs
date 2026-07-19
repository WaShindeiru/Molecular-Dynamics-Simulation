use std::collections::HashMap;

use nalgebra::Vector3;

use crate::sim_core::world::boundary_constraint::{Compliance, EdgeCondition, ParticleCompliance};
use crate::sim_core::world::boxed_world::box_task::handle_task::handle_partial_velocity_step::apply_velocity_constraint;
use crate::sim_core::world::boxed_world::box_task::{
  ForceTaskParticleData, ForceTaskResult, VelocityTaskParticleData, VelocityTaskResult,
};
use crate::sim_core::world::boxed_world::integration::verlet_nose_hoover::computation::HalfVelocityResult;
use crate::sim_core::world::cell::{FixedPositionParticle, LinkedCellContainer};
use crate::sim_core::world::computation::FP;
use crate::sim_core::world::optimized_world::computation::{
  compute_forces_potential, verlet_noose_hoover_half_velocity_position,
};
use crate::sim_core::world::optimized_world::handle_task::handle_partial_velocity_step::handle_partial_velocity_step;

pub(crate) mod handle_partial_velocity_step;

pub fn handle_velocity_batch_task(
  task_id: usize,
  cell_ids: &[usize],
  history: &LinkedCellContainer,
  time_step: f64,
  previous_thermostat_epsilon: f64,
  current_iteration: usize,
) -> VelocityTaskResult {
  let edge_condition = history.edge_condition();
  let container_size = history.config().world_size;

  let ids: Vec<usize> = cell_ids.iter().flat_map(|&cell_id| history.ids_in_cell(cell_id)).collect();

  let (normal_ids, custom_vel_ids): (Vec<usize>, Vec<usize>) = ids
    .into_iter()
    .partition(|&id| !history.particles().get(id).unwrap().is_custom_velocity_atom());

  let HalfVelocityResult { half_velocity, new_position, thermostat_work, compliance } =
    verlet_noose_hoover_half_velocity_position(
      history,
      &normal_ids,
      time_step,
      previous_thermostat_epsilon,
      current_iteration,
      &container_size,
      edge_condition,
    );

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
    }
    EdgeCondition::PeriodicAll => (all_normal_particles, HashMap::new()),
  };

  let mut result_particles = apply_velocity_constraint(compliant, edge_condition);

  if !non_compliant.is_empty() {
    if edge_condition.collision_split_enabled() {
      let corrected = handle_partial_velocity_step(
        history,
        non_compliant,
        container_size,
        edge_condition,
        previous_thermostat_epsilon,
        time_step,
      );
      result_particles.extend(corrected);
    } else {
      result_particles.extend(apply_velocity_constraint(non_compliant, edge_condition));
    }
  }

  for id in custom_vel_ids {
    let particle = history.particles().get(id).unwrap();
    let vel = *particle.get_velocity();
    let new_position = particle.get_position() + vel * time_step;
    result_particles.insert(id, VelocityTaskParticleData {
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

  VelocityTaskResult { task_id, particles: result_particles }
}

pub fn handle_force_batch_task(
  task_id: usize,
  cell_ids: &[usize],
  integration_cache: &LinkedCellContainer,
  fp: &mut Vec<FP>,
  gradients_cache: &mut Vec<Vector3<f64>>,
) -> ForceTaskResult {
  let mut particles: HashMap<usize, ForceTaskParticleData> = HashMap::new();
  let mut potential_energy_total = 0.0f64;

  for &cell_id in cell_ids {
    let particles_i: Vec<FixedPositionParticle> = integration_cache.atoms_for_cell_fixed(cell_id);

    let mut particles_j: Vec<FixedPositionParticle> =
      integration_cache.neighbour_atoms_periodic_fixed_positions(cell_id);
    particles_j.extend(integration_cache.atoms_for_cell_fixed(cell_id));

    let potential_energy =
      compute_forces_potential(&particles_i, &particles_j, integration_cache, fp, gradients_cache);
    potential_energy_total += potential_energy;

    for particle in particles_j.iter() {
      let id = particle.id;
      let fp_entry = &fp[id];
      particles
        .entry(id)
        .and_modify(|entry| {
          entry.force += fp_entry.force;
          entry.potential_energy += fp_entry.potential_energy;
        })
        .or_insert(ForceTaskParticleData {
          box_id: cell_id,
          force: fp_entry.force,
          potential_energy: fp_entry.potential_energy,
        });
    }
  }

  ForceTaskResult {
    task_id,
    potential_energy: potential_energy_total,
    particles,
  }
}
