use std::collections::HashMap;
use std::sync::Arc;

use nalgebra::Vector3;

use crate::particle::Particle;
use crate::sim_core::world::boundary_constraint::EdgeCondition;
use crate::sim_core::world::boundary_constraint::periodic::apply_velocity_constraint_periodic;
use crate::sim_core::world::boxed_world::box_container::{BoxContainer, sim_box::SimulationBox};
use crate::sim_core::world::boxed_world::box_task::VelocityTaskParticleData;
use crate::sim_core::world::boxed_world::box_task::handle_task::partial_velocity_step::PartialVelocityStep;
use crate::sim_core::world::boxed_world::integration::verlet_nose_hoover::computation::HalfVelocityResult;

pub fn apply_velocity_constraint(
  particles: HashMap<usize, VelocityTaskParticleData>,
  edge_condition: EdgeCondition,
) -> HashMap<usize, VelocityTaskParticleData> {
  
  match edge_condition {
    EdgeCondition::Simple { .. } => unimplemented!("not yet"),
    EdgeCondition::Periodic { .. } => {
      particles
        .into_iter()
        .map(|(id, mut d)| {
          d.half_velocity = apply_velocity_constraint_periodic(&d.compliance, &d.half_velocity);
          (id, d)
        })
        .collect()
    },
    EdgeCondition::PeriodicAll => particles,
  }
}

pub fn handle_partial_velocity_step(
  history: &BoxContainer<Arc<SimulationBox>>,
  non_compliant: HashMap<usize, VelocityTaskParticleData>,
  world_size: Vector3<f64>,
  edge_condition: EdgeCondition,
  thermostat_epsilon: f64,
  time_step: f64,
) -> HashMap<usize, VelocityTaskParticleData> {
  let trigger_small_subtask_size = match edge_condition {
    EdgeCondition::Simple { trigger_small_subtask_size, .. }
    | EdgeCondition::Periodic { trigger_small_subtask_size, .. } => trigger_small_subtask_size,
    EdgeCondition::PeriodicAll => panic!("can't call handle_partial_velocity_step for PeriodicAll"),
  };

  log::info!("doing handle_partial_velocity_step with {} subtasks.", trigger_small_subtask_size);

  let particles: HashMap<usize, Particle> = non_compliant.keys()
    .map(|id| (*id, Arc::unwrap_or_clone(history.get_particle(*id))))
    .collect();

  let compliance = non_compliant.iter()
    .map(|(id, d)| (*id, d.compliance))
    .collect();

  let HalfVelocityResult { half_velocity, new_position, thermostat_work, compliance } =
    PartialVelocityStep::new(
      world_size,
      edge_condition,
      particles,
      compliance,
      thermostat_epsilon,
      trigger_small_subtask_size,
      time_step,
    ).run();

  compliance
    .into_iter()
    .map(|(id, compliance)| (id, VelocityTaskParticleData {
      half_velocity: half_velocity[&id],
      new_position: new_position[&id],
      thermostat_work: thermostat_work[&id],
      compliance,
    }))
    .collect()
}
