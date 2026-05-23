pub mod box_container;

use crate::persistence::dto::world::history::HistoryDTO;
use crate::sim_core::world::thermostat::IntegrationAlgorithm;
use nalgebra::Vector3;

pub struct BoxedWorldDTO {
  pub num_of_atoms: usize,
  pub size: Vector3<f64>,
  pub history: HistoryDTO,
  pub integration_algorithm: IntegrationAlgorithm,

  pub num_of_world_iterations: usize,
  pub number_of_resets: usize,
  pub max_iteration_till_reset: usize,

  pub laamps_frame_iteration_count: usize,
  pub energy_frame_iteration_count: usize,
}
