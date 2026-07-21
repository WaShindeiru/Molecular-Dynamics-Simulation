use std::sync::{Arc, Mutex};

use crate::sim_core::world::WorldType;
use crate::sim_core::world::boundary_constraint::EdgeCondition;
use crate::sim_core::world::gravity::GravityManager;
use crate::sim_core::world::thermostat::IntegrationAlgorithm;
use crate::sim_core::world::saver::SaveOptions;
use nalgebra::Vector3;

#[derive(Clone)]
pub struct SimulationConfig {
  pub world_size: Vector3<f64>,
  pub gravity_manager: Arc<Mutex<GravityManager>>,
  pub time_step: f64,
  pub num_of_iterations: usize,
  pub max_iteration_till_reset: usize,
  pub save_options: SaveOptions,
  pub integration_algorithm: IntegrationAlgorithm,
  pub world_type: WorldType,
  pub edge_condition: EdgeCondition,
  pub optimization: bool,
  /// P-controller gain for VelocityControlledParticle (OptimizedWorld only).
  pub alpha: f64,
}

impl SimulationConfig {
  pub fn new(
    world_size: Vector3<f64>,
    gravity_schedule: Vec<(usize, f64)>,
    time_step: f64,
    num_of_iterations: usize,
    max_iteration_till_reset: usize,
    save_options: SaveOptions,
    integration_algorithm: IntegrationAlgorithm,
    world_type: WorldType,
    edge_condition: EdgeCondition,
    optimization: bool,
    alpha: f64,
  ) -> Self {
    SimulationConfig {
      world_size,
      gravity_manager: Arc::new(Mutex::new(GravityManager::new(gravity_schedule))),
      time_step,
      num_of_iterations,
      max_iteration_till_reset,
      save_options,
      integration_algorithm,
      world_type,
      edge_condition,
      optimization,
      alpha,
    }
  }

  pub fn initial_gravity(&self) -> f64 {
    self
      .gravity_manager
      .lock()
      .expect("gravity manager lock poisoned")
      .initial_gravity()
  }
}
