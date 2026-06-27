use nalgebra::Vector3;
use std::collections::HashMap;
use std::sync::Arc;


use crate::sim_core::world::boundary_constraint::{EdgeCondition, ParticleCompliance};

use crate::sim_core::world::boxed_world::box_container::BoxContainer;
use crate::sim_core::world::boxed_world::box_container::sim_box::SimulationBox;
use crate::sim_core::world::boxed_world::integration_cache::IntegrationCache;

mod force_task_box_container;
mod handle_task;
pub mod task_manager;

pub enum BoxTask {
  VelocityBatchTask {
    task_id: usize,
    box_ids: Arc<Vec<usize>>,
    history: Arc<BoxContainer<Arc<SimulationBox>>>,
    time_step: f64,
    previous_thermostat_epsilon: f64,
    current_iteration: usize,
    container_size: Vector3<f64>,
    edge_condition: EdgeCondition,
  },
  ForceBatchTask {
    task_id: usize,
    boundary_condition: EdgeCondition,
    box_ids: Arc<Vec<usize>>,
    integration_cache: Arc<IntegrationCache>,
  },
}

pub struct VelocityTaskParticleData {
  pub half_velocity: Vector3<f64>,
  pub new_position: Vector3<f64>,
  pub thermostat_work: f64,
  pub compliance: ParticleCompliance,
}

pub struct VelocityTaskResult {
  pub task_id: usize,
  pub particles: HashMap<usize, VelocityTaskParticleData>,
}

pub struct ForceTaskParticleData {
  pub box_id: usize,
  pub force: Vector3<f64>,
  pub potential_energy: f64,
}

pub struct ForceTaskResult {
  pub task_id: usize,
  pub potential_energy: f64,
  pub particles: HashMap<usize, ForceTaskParticleData>,
}

pub enum BoxResult {
  VelocityResult(VelocityTaskResult),
  ForceResult(ForceTaskResult),
}
