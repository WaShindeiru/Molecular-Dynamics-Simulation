use std::collections::HashMap;
use std::sync::{Arc, RwLock};
use nalgebra::Vector3;
use crate::particle::Particle;
use crate::sim_core::world::boundary_constraint::{EdgeCondition, ParticleCompliance};
use crate::sim_core::world::boxed_world::box_container::BoxContainer;
use crate::sim_core::world::boxed_world::box_container::sim_box::SimulationBox;
use crate::sim_core::world::boxed_world::history_manager::HistoryManager;
use crate::sim_core::world::boxed_world::history_manager::history_manager::HistoryManager;
use crate::sim_core::world::boxed_world::integration::verlet_nose_hoover::computation::FPInfoBoxed;
use crate::sim_core::world::boxed_world::integration_cache::IntegrationCache;

pub mod task_manager;
mod handle_task;

pub enum BoxTask {
  // VelocityTask {
  //   history_manager: Arc<RwLock<HistoryManager>>,
  //   box_id: usize,
  //   time_step: f64,
  //   previous_thermostat_epsilon: f64,
  //   current_iteration: usize,
  //   container_size: Vector3<f64>,
  //   edge_condition: EdgeCondition,
  // },
  // ForceTask {
  //   history_manager: Arc<RwLock<HistoryManager>>,
  //   box_id: usize,
  // },
  VelocityBatchTask {
    task_id: usize,
    box_ids: Vec<usize>,
    history: Arc<BoxContainer<Arc<SimulationBox>>>,
    time_step: f64,
    previous_thermostat_epsilon: f64,
    current_iteration: usize,
    container_size: Vector3<f64>,
    edge_condition: EdgeCondition,
  },
  ForceBatchTask {
    task_id: usize,
    box_ids: Vec<usize>,
    integration_cache: Arc<IntegrationCache>
  }
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
  pub acceleration: Vector3<f64>,
  pub potential_energy: f64,
}

pub struct ForceTaskResult {
  pub task_id: usize,
  pub potential_energy: f64,
  pub optimization_considered: usize,
  pub optimization_ignored: usize,
  pub particles: HashMap<usize, ForceTaskParticleData>
}

pub enum BoxResult {
  VelocityResult (VelocityTaskResult),
  ForceResult (ForceTaskResult),
}