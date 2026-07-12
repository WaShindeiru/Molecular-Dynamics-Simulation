use log::debug;
use std::collections::HashMap;
use std::sync::Arc;
use std::sync::mpsc::{Receiver, RecvTimeoutError, Sender};
use std::thread::JoinHandle;
use std::time::Duration;
use crate::perf_log;
use crate::data::SimulationConfig;
use crate::sim_core::world::boxed_world::box_container::BoxContainer;
use crate::sim_core::world::cell::box_container_config::BoxContainerConfig;
use crate::sim_core::world::boxed_world::box_container::sim_box::SimulationBox;
use crate::sim_core::world::boxed_world::box_task::task_manager::threads::create_threads;
use crate::sim_core::world::boxed_world::box_task::{BoxResult, BoxTask};
use crate::sim_core::world::boxed_world::computation_collector::ComputationCollector;
use crate::sim_core::world::boxed_world::integration_cache::IntegrationCache;
use crate::sim_core::world::boxed_world::integration_cache::integration_cache_builder::IntegrationCacheBuilder;
use crate::sim_core::world::cell::{TaskSplitVariant, TaskSplitter};

mod threads;

#[derive(Debug, Clone, Copy)]
pub struct TaskManagerConfig {
  pub debug: bool,
  pub task_worker_multiplier: f64,
  pub split: TaskSplitVariant,
}

pub struct TaskManager {
  simulation_config: SimulationConfig,
  container_config: BoxContainerConfig,

  debug: bool,
  threads: Vec<JoinHandle<()>>,
  tx_task: Sender<BoxTask>,
  rx_result: Receiver<BoxResult>,
  num_workers: usize,
  task_worker_multiplier: f64,
  task_splitter: TaskSplitter,

  task_box_mapping: Option<HashMap<usize, Arc<Vec<usize>>>>,
}

impl TaskManager {
  pub fn new(
    debug: bool,
    simulation_config: SimulationConfig,
    container_config: BoxContainerConfig,
    task_worker_multiplier: f64,
    split_variant: TaskSplitVariant,
    num_atoms: usize,
  ) -> Self {
    let (tx_task, rx_result, threads, num_workers) = create_threads(debug, num_atoms);

    TaskManager {
      debug,
      simulation_config,
      container_config,
      threads,
      tx_task,
      rx_result,
      num_workers,
      task_worker_multiplier,
      task_splitter: TaskSplitter::new(split_variant),
      task_box_mapping: Option::None,
    }
  }

  pub fn num_workers(&self) -> usize {
    self.num_workers
  }

  fn clear(&mut self) {
    self.task_box_mapping = None
  }

  pub fn split_into_tasks_multiplier(&mut self, config: &BoxContainerConfig) {
    assert!(self.task_box_mapping.is_none());
    let num_of_tasks = (self.num_workers as f64 * self.task_worker_multiplier).floor() as usize;
    self.task_box_mapping = Some(self.task_splitter.split(num_of_tasks, config));
  }

  pub fn task_box_mapping(&self) -> Option<&HashMap<usize, Arc<Vec<usize>>>> {
    self.task_box_mapping.as_ref()
  }

  pub fn half_velocity_step(
    &self,
    box_container: Arc<BoxContainer<Arc<SimulationBox>>>,
    thermostat_epsilon: f64,
    current_iteration: usize,
  ) -> Arc<IntegrationCache> {
    let mapping = self
      .task_box_mapping
      .as_ref()
      .expect("split_into_tasks must be called before half_velocity_step");

    let num_tasks = mapping.len();
    perf_log!("Before doing all_praticles_reset on box_container");
    let particles = box_container.all_particles_reset();
    perf_log!("After doing all_praticles_reset on box_container");
    let mut builder = IntegrationCacheBuilder::new(
      self.simulation_config.clone(),
      self.container_config,
      particles);

    for (task_id, box_ids) in mapping {
      let task = BoxTask::VelocityBatchTask {
        task_id: *task_id,
        box_ids: Arc::clone(box_ids),
        history: Arc::clone(&box_container),
        time_step: self.simulation_config.time_step,
        previous_thermostat_epsilon: thermostat_epsilon,
        current_iteration,
        container_size: self.container_config.world_size,
        edge_condition: self.simulation_config.edge_condition,
      };
      self.tx_task.send(task).unwrap();
      debug!("Sent VelocityBatchTask for task_id {}", task_id);
      perf_log!("Sent VelocityBatchTask for task_id {}", task_id);
    }

    for i in 0..num_tasks {
      match self.rx_result.recv_timeout(Duration::from_secs(20)) {
        Ok(BoxResult::VelocityResult(result)) => {
          perf_log!("Received VelocityResult for task_id {}", result.task_id);
          debug!("Received VelocityResult for task_id {}", result.task_id);
          builder.add_velocity_results(result.particles);
        }
        Ok(_) => panic!("Expected VelocityResult, got wrong result type"),
        Err(RecvTimeoutError::Timeout) => {
          panic!("Velocity step timed out after 20 seconds")
        }
        Err(RecvTimeoutError::Disconnected) => panic!("Worker channel disconnected"),
      }
    }

    Arc::new(
      builder
        .build()
        .expect("Not all particles received velocity results"),
    )
  }

  pub fn force_step(&self, integration_cache: Arc<IntegrationCache>) -> ComputationCollector {
    let mapping = self
      .task_box_mapping
      .as_ref()
      .expect("split_into_tasks must be called before force_step");

    let num_tasks = mapping.len();
    perf_log!("Before creating ComputationCollector from integration cache as part of Force step");
    let mut collector = ComputationCollector::from_integration_cache(
      self.simulation_config.clone(),
      Arc::clone(&integration_cache),
    );
    perf_log!("After creating ComputationCollector");

    for (task_id, box_ids) in mapping {
      let task = BoxTask::ForceBatchTask {
        task_id: *task_id,
        boundary_condition: self.simulation_config.edge_condition,
        box_ids: Arc::clone(box_ids),
        integration_cache: Arc::clone(&integration_cache),
      };
      self.tx_task.send(task).unwrap();
      perf_log!("Sent ForceBatchTask for task_id {}", task_id);
      debug!("Sent ForceBatchTask for task_id {}", task_id);
    }

    for _ in 0..num_tasks {
      match self.rx_result.recv_timeout(Duration::from_secs(20)) {
        Ok(BoxResult::ForceResult(result)) => {
          perf_log!("Received ForceResult for task_id {}", result.task_id);
          debug!("Received ForceResult for task_id {}", result.task_id);

          #[cfg(debug_assertions)]
          for (particle_id, particle_data) in result.particles.iter() {
            let position = integration_cache
              .box_cache()
              .get_box(particle_data.box_id)
              .particle(*particle_id)
              .get_position()
              .clone();

            debug!(
              "Received ForceResult task {} particle {} position {:?} force {:?} potential_energy {}",
              result.task_id,
              particle_id,
              position,
              particle_data.force,
              particle_data.potential_energy
            );
          }

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
