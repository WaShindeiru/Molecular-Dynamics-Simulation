use std::collections::HashMap;
use std::sync::Arc;
use nalgebra::Vector3;
use crate::sim_core::world::boundary_constraint::EdgeCondition;
use crate::sim_core::world::boxed_world::box_container::BoxContainer;
use crate::sim_core::world::boxed_world::box_container::sim_box::SimulationBox;
use crate::sim_core::world::boxed_world::box_task::{VelocityTaskParticleData, VelocityTaskResult};
use crate::sim_core::world::boxed_world::integration::verlet_nose_hoover::computation::verlet_noose_hoover_half_velocity_position;

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

  let particles: HashMap<usize, VelocityTaskParticleData> = half_velocity_response.half_velocity
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
