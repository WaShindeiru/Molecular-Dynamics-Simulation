use std::sync::{Arc, Mutex};

use crate::sim_core::world::WorldType;
use crate::sim_core::world::boundary_constraint::EdgeCondition;
use crate::sim_core::world::gravity::GravityManager;
use crate::sim_core::world::timestep::TimestepManager;
use crate::sim_core::world::thermostat::IntegrationAlgorithm;
use crate::sim_core::world::saver::SaveOptions;
use nalgebra::Vector3;

#[derive(Clone)]
pub struct SimulationConfig {
  pub world_size: Vector3<f64>,
  pub gravity_manager: Arc<Mutex<GravityManager>>,
  pub timestep_manager: Arc<Mutex<TimestepManager>>,
  pub num_of_iterations: usize,
  pub max_iteration_till_reset: usize,
  pub save_options: SaveOptions,
  pub integration_algorithm: IntegrationAlgorithm,
  pub world_type: WorldType,
  pub edge_condition: EdgeCondition,
}

impl SimulationConfig {
  pub fn new(
    world_size: Vector3<f64>,
    gravity_schedule: Vec<(usize, f64)>,
    timestep_schedule: Vec<(usize, f64)>,
    num_of_iterations: usize,
    max_iteration_till_reset: usize,
    save_options: SaveOptions,
    integration_algorithm: IntegrationAlgorithm,
    world_type: WorldType,
    edge_condition: EdgeCondition,
  ) -> Self {
    SimulationConfig {
      world_size,
      gravity_manager: Arc::new(Mutex::new(GravityManager::new(gravity_schedule))),
      timestep_manager: Arc::new(Mutex::new(TimestepManager::new(timestep_schedule))),
      num_of_iterations,
      max_iteration_till_reset,
      save_options,
      integration_algorithm,
      world_type,
      edge_condition,
    }
  }

  pub fn initial_gravity(&self) -> f64 {
    self
      .gravity_manager
      .lock()
      .expect("gravity manager lock poisoned")
      .initial_gravity()
  }

  pub fn initial_time_step(&self) -> f64 {
    self
      .timestep_manager
      .lock()
      .expect("timestep manager lock poisoned")
      .initial_timestep()
  }
}
