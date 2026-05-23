use crate::sim_core::world::WorldType;
use crate::sim_core::world::boundary_constraint::EdgeCondition;
use crate::sim_core::world::thermostat::IntegrationAlgorithm;
use crate::sim_core::world::saver::SaveOptions;
use nalgebra::Vector3;

#[derive(Clone)]
pub struct SimulationConfig {
  pub world_size: Vector3<f64>,
  pub potential_gravity_max: f64,
  pub time_step: f64,
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
    potential_gravity_max: f64,
    time_step: f64,
    num_of_iterations: usize,
    max_iteration_till_reset: usize,
    save_options: SaveOptions,
    integration_algorithm: IntegrationAlgorithm,
    world_type: WorldType,
    edge_condition: EdgeCondition,
  ) -> Self {
    SimulationConfig {
      world_size,
      potential_gravity_max,
      time_step,
      num_of_iterations,
      max_iteration_till_reset,
      save_options,
      integration_algorithm,
      world_type,
      edge_condition,
    }
  }
}
