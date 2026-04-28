use std::collections::HashMap;
use std::ops::Deref;
use std::sync::Arc;
use std::thread::JoinHandle;
use std::sync::mpsc::{Receiver, RecvTimeoutError, Sender};
use std::time::Duration;
use log::debug;
use nalgebra::Vector3;
use crate::data::SimulationConfig;
use crate::sim_core::world::boxed_world::box_container::BoxContainer;
use crate::sim_core::world::boxed_world::box_container::box_container_config::BoxContainerConfig;
use crate::sim_core::world::boxed_world::box_container::sim_box::{SimulationBox, get_id_simulation_box};
use crate::sim_core::world::boxed_world::box_task::{BoxResult, BoxTask};
use crate::sim_core::world::boxed_world::box_task::task_manager::threads::create_threads;
use crate::sim_core::world::boxed_world::integration_cache::IntegrationCache;
use crate::sim_core::world::boxed_world::integration_cache::integration_cache_builder::IntegrationCacheBuilder;

mod threads;

pub struct TaskManager {
  simulation_config: SimulationConfig,
  container_config: BoxContainerConfig,
  threads: Vec<JoinHandle<()>>,
  tx_task: Sender<BoxTask>,
  rx_result: Receiver<BoxResult>,
  num_workers: usize,
  task_box_mapping: Option<HashMap<usize, Vec<usize>>>
}

impl TaskManager {
  pub fn new(debug: bool, simulation_config: SimulationConfig, container_config: BoxContainerConfig) -> Self {
    let (tx_task, rx_result, threads, num_workers) = create_threads(debug);

    TaskManager {
      simulation_config,
      container_config,
      threads,
      tx_task,
      rx_result,
      num_workers,
      task_box_mapping: Option::None,
    }
  }

  pub fn num_workers(&self) -> usize {
    self.num_workers
  }

  fn clear(&mut self) {
    self.task_box_mapping = None
  }

  pub fn split_into_tasks<B, C>(&mut self, num_of_task: usize, container: C)
  where
    C: Deref<Target = BoxContainer<B>>,
  {
    assert!(self.task_box_mapping.is_none());

    let config = container.config();
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
        .flat_map(|z| (0..ny).flat_map(move |y| (0..nx).map(move |x| {
          get_id_simulation_box(&Vector3::new(x, y, z), &config.box_count_dim)
        })))
        .collect();

      mapping.insert(task_id, box_ids);
      z_start = z_end;
    }

    self.task_box_mapping = Some(mapping);
  }

  pub fn task_box_mapping(&self) -> Option<&HashMap<usize, Vec<usize>>> {
    self.task_box_mapping.as_ref()
  }

  pub fn half_velocity_step(
    &mut self,
    box_container: Arc<BoxContainer<Arc<SimulationBox>>>,
    thermostat_epsilon: f64,
    current_iteration: usize,
  ) -> IntegrationCache {
    let mapping = self.task_box_mapping.as_ref()
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
      match self.rx_result.recv_timeout(Duration::from_secs(60)) {
        Ok(BoxResult::VelocityResult(result)) => {
          debug!("Received VelocityResult for task_id {}", result.task_id);
          builder.add_velocity_results(result.particles);
        }
        Ok(_) => panic!("Expected VelocityResult, got wrong result type"),
        Err(RecvTimeoutError::Timeout) => panic!("Velocity step timed out after 20 seconds"),
        Err(RecvTimeoutError::Disconnected) => panic!("Worker channel disconnected"),
      }
    }

    builder.build().expect("Not all particles received velocity results")
  }
}