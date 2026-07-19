use std::collections::HashMap;

use nalgebra::Vector3;

use crate::sim_core::world::boundary_constraint::EdgeCondition;
use crate::sim_core::world::boxed_world::box_task::VelocityTaskParticleData;
use crate::sim_core::world::boxed_world::box_task::handle_task::partial_velocity_step::PartialVelocityStep;
use crate::sim_core::world::boxed_world::integration::verlet_nose_hoover::computation::HalfVelocityResult;
use crate::sim_core::world::cell::{LinkedCellContainer};

pub fn handle_partial_velocity_step(
  history: &LinkedCellContainer,
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

  let initial_particles = non_compliant
    .keys()
    .map(|id| (*id, history.particles().get(*id).unwrap().clone()))
    .collect();
  let initial_compliance = non_compliant
    .iter()
    .map(|(id, d)| (*id, d.compliance))
    .collect();

  let HalfVelocityResult { half_velocity, new_position, thermostat_work, compliance } =
    PartialVelocityStep::new(
      world_size,
      edge_condition,
      initial_particles,
      initial_compliance,
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
