mod threads;

use std::collections::HashMap;
use std::sync::Arc;
use std::sync::mpsc::{Receiver, RecvTimeoutError, Sender};
use std::thread::JoinHandle;
use std::time::Duration;

use crate::data::SimulationConfig;
use crate::sim_core::world::cell::{TaskSplitVariant, TaskSplitter};
use crate::sim_core::world::linked_cell_world::LinkedCellContainer;
use crate::sim_core::world::linked_cell_world::computation_collector::ComputationCollector;
use crate::sim_core::world::linked_cell_world::integration_cache::IntegrationCache;
use crate::sim_core::world::linked_cell_world::integration_cache::integration_cache_builder::IntegrationCacheBuilder;
use crate::sim_core::world::linked_cell_world::linked_cell_task::{LinkedCellResult, LinkedCellTask};
use threads::create_threads;

pub struct TaskManager {
  simulation_config: SimulationConfig,
  threads: Vec<JoinHandle<()>>,
  tx_task: Sender<LinkedCellTask>,
  rx_result: Receiver<LinkedCellResult>,
  num_workers: usize,
  task_worker_multiplier: f64,
  task_splitter: TaskSplitter,
  task_cell_mapping: Option<HashMap<usize, Arc<Vec<usize>>>>,
}

impl TaskManager {
  pub fn new(
    debug: bool,
    simulation_config: SimulationConfig,
    task_worker_multiplier: f64,
    split_variant: TaskSplitVariant,
    num_atoms: usize,
  ) -> Self {
    let (tx_task, rx_result, threads, num_workers) = create_threads(debug, num_atoms);
    TaskManager {
      simulation_config,
      threads,
      tx_task,
      rx_result,
      num_workers,
      task_worker_multiplier,
      task_splitter: TaskSplitter::new(split_variant),
      task_cell_mapping: None,
    }
  }

  pub fn num_workers(&self) -> usize {
    self.num_workers
  }

  fn clear(&mut self) {
    self.task_cell_mapping = None;
  }

  pub fn split_into_tasks_multiplier(&mut self, container: &LinkedCellContainer) {
    assert!(self.task_cell_mapping.is_none());
    let num_of_tasks = (self.num_workers as f64 * self.task_worker_multiplier).floor() as usize;
    self.task_cell_mapping = Some(self.task_splitter.split(num_of_tasks, container.config()));
  }

  pub fn task_cell_mapping(&self) -> Option<&HashMap<usize, Arc<Vec<usize>>>> {
    self.task_cell_mapping.as_ref()
  }

  pub fn half_velocity_step(
    &self,
    container: Arc<LinkedCellContainer>,
    thermostat_epsilon: f64,
    current_iteration: usize,
    time_step: f64,
  ) -> IntegrationCache {
    let mapping = self
      .task_cell_mapping
      .as_ref()
      .expect("split_into_tasks must be called before half_velocity_step");

    let num_tasks = mapping.len();
    let mut builder = IntegrationCacheBuilder::new(Arc::clone(&container));

    for (task_id, cell_ids) in mapping {
      let task = LinkedCellTask::VelocityBatchTask {
        task_id: *task_id,
        cell_ids: Arc::clone(cell_ids),
        history: Arc::clone(&container),
        time_step,
        previous_thermostat_epsilon: thermostat_epsilon,
        current_iteration,
      };
      self.tx_task.send(task).unwrap();
    }

    for _ in 0..num_tasks {
      match self.rx_result.recv_timeout(Duration::from_secs(20)) {
        Ok(LinkedCellResult::VelocityResult(result)) => {
          builder.add_velocity_results(result.particles);
        }
        Ok(_) => panic!("Expected VelocityResult, got wrong result type"),
        Err(RecvTimeoutError::Timeout) => panic!("Velocity step timed out after 20 seconds"),
        Err(RecvTimeoutError::Disconnected) => panic!("Worker channel disconnected"),
      }
    }

    builder.build().expect("Not all particles received velocity results")
  }

  pub fn force_step(&self, integration_cache: IntegrationCache) -> ComputationCollector {
    let mapping = self
      .task_cell_mapping
      .as_ref()
      .expect("split_into_tasks must be called before force_step");

    let num_tasks = mapping.len();
    let IntegrationCache { local_container, read_container, half_velocity_cache, particle_compliance } =
      integration_cache;

    let num_particles = half_velocity_cache.len();
    let mut collector = ComputationCollector::from_integration_cache(
      self.simulation_config.clone(),
      num_particles,
      half_velocity_cache,
      particle_compliance,
      local_container,
    );

    for (task_id, cell_ids) in mapping {
      let task = LinkedCellTask::ForceBatchTask {
        task_id: *task_id,
        cell_ids: Arc::clone(cell_ids),
        integration_cache: Arc::clone(&read_container),
      };
      self.tx_task.send(task).unwrap();
    }

    for _ in 0..num_tasks {
      match self.rx_result.recv_timeout(Duration::from_secs(60)) {
        Ok(LinkedCellResult::ForceResult(result)) => {
          collector.apply_force_results(&result.particles);
        }
        Ok(_) => panic!("Expected ForceResult, got wrong result type"),
        Err(RecvTimeoutError::Timeout) => panic!("Force step timed out after 60 seconds"),
        Err(RecvTimeoutError::Disconnected) => panic!("Worker channel disconnected"),
      }
    }

    #[cfg(debug_assertions)]
    collector.assert_zero_net_force();

    collector
  }
}
