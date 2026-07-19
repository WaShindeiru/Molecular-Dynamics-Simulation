use std::sync::Arc;

use crate::sim_core::world::boxed_world::box_task::{ForceTaskResult, VelocityTaskResult};
use crate::sim_core::world::cell::LinkedCellContainer;

pub enum OptimizedTask {
  VelocityBatchTask {
    task_id: usize,
    cell_ids: Arc<Vec<usize>>,
    history: Arc<LinkedCellContainer>,
    time_step: f64,
    previous_thermostat_epsilon: f64,
    current_iteration: usize,
  },
  ForceBatchTask {
    task_id: usize,
    cell_ids: Arc<Vec<usize>>,
    integration_cache: Arc<LinkedCellContainer>,
  },
}

pub enum OptimizedResult {
  VelocityResult(VelocityTaskResult),
  ForceResult(ForceTaskResult),
}
