use log::debug;
use nalgebra::Vector3;
use std::collections::HashMap;
use std::sync::Arc;
use std::sync::mpsc::{Receiver, RecvTimeoutError, Sender};
use std::thread::JoinHandle;
use std::time::Duration;

use crate::data::SimulationConfig;
use crate::sim_core::world::boxed_world::box_container::BoxContainer;
use crate::sim_core::world::boxed_world::box_container::box_container_config::BoxContainerConfig;
use crate::sim_core::world::boxed_world::box_container::sim_box::{
  SimulationBox, get_id_simulation_box,
};
use crate::sim_core::world::boxed_world::box_task::task_manager::threads::create_threads;
use crate::sim_core::world::boxed_world::box_task::{BoxResult, BoxTask};
use crate::sim_core::world::boxed_world::computation_collector::ComputationCollector;
use crate::sim_core::world::boxed_world::integration_cache::IntegrationCache;
use crate::sim_core::world::boxed_world::integration_cache::integration_cache_builder::IntegrationCacheBuilder;

mod threads;

#[derive(Debug, Clone, Copy)]
pub struct TaskManagerConfig {
  pub debug: bool,
  pub task_worker_multiplier: f64,
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

  task_box_mapping: Option<HashMap<usize, Vec<usize>>>,
}

impl TaskManager {
  pub fn new(
    debug: bool,
    simulation_config: SimulationConfig,
    container_config: BoxContainerConfig,
    task_worker_multiplier: f64,
  ) -> Self {
    let (tx_task, rx_result, threads, num_workers) = create_threads(debug);

    TaskManager {
      debug,
      simulation_config,
      container_config,
      threads,
      tx_task,
      rx_result,
      num_workers,
      task_worker_multiplier,
      task_box_mapping: Option::None,
    }
  }

  pub fn num_workers(&self) -> usize {
    self.num_workers
  }

  fn clear(&mut self) {
    self.task_box_mapping = None
  }

  pub fn split_into_tasks(&mut self, num_of_task: usize, config: &BoxContainerConfig) {
    assert!(self.task_box_mapping.is_none());
    let nx = config.box_count_dim.x;
    let ny = config.box_count_dim.y;
    let nz = config.box_count_dim.z;

    let z_per_task = nz / num_of_task;
    let remainder = nz % num_of_task;

    let mut mapping: HashMap<usize, Vec<usize>> = HashMap::with_capacity(num_of_task);
    let mut z_start = 0;

    for task_id in 0..num_of_task {
      let z_count = z_per_task + if task_id < remainder { 1 } else { 0 };
      let z_end = z_start + z_count;

      let box_ids = (z_start..z_end)
        .flat_map(|z| {
          (0..ny).flat_map(move |y| {
            (0..nx)
              .map(move |x| get_id_simulation_box(&Vector3::new(x, y, z), &config.box_count_dim))
          })
        })
        .collect();

      mapping.insert(task_id, box_ids);
      z_start = z_end;
    }

    self.task_box_mapping = Some(mapping);
  }

  pub fn split_into_tasks_multiplier(&mut self, config: &BoxContainerConfig) {
    self.split_into_tasks(
      (self.num_workers as f64 * self.task_worker_multiplier).floor() as usize,
      config,
    );
  }

  pub fn task_box_mapping(&self) -> Option<&HashMap<usize, Vec<usize>>> {
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
    let particles = box_container.all_particles_reset();
    let mut builder = IntegrationCacheBuilder::new(self.container_config, particles);

    for (task_id, box_ids) in mapping {
      let task = BoxTask::VelocityBatchTask {
        task_id: *task_id,
        box_ids: box_ids.clone(),
        history: Arc::clone(&box_container),
        time_step: self.simulation_config.time_step,
        previous_thermostat_epsilon: thermostat_epsilon,
        current_iteration,
        container_size: self.container_config.world_size,
        edge_condition: self.simulation_config.edge_condition,
      };
      self.tx_task.send(task).unwrap();
      debug!("Sent VelocityBatchTask for task_id {}", task_id);
    }

    for _ in 0..num_tasks {
      match self.rx_result.recv_timeout(Duration::from_secs(20)) {
        Ok(BoxResult::VelocityResult(result)) => {
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
    let mut collector = ComputationCollector::from_integration_cache(
      self.simulation_config.clone(),
      Arc::clone(&integration_cache),
    );

    for (task_id, box_ids) in mapping {
      let task = BoxTask::ForceBatchTask {
        task_id: *task_id,
        boundary_condition: self.simulation_config.edge_condition,
        box_ids: box_ids.clone(),
        integration_cache: Arc::clone(&integration_cache),
      };
      self.tx_task.send(task).unwrap();
      debug!("Sent ForceBatchTask for task_id {}", task_id);
    }

    for _ in 0..num_tasks {
      match self.rx_result.recv_timeout(Duration::from_secs(20)) {
        Ok(BoxResult::ForceResult(result)) => {
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

    collector
  }
}
