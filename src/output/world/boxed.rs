pub mod box_container;

use nalgebra::Vector3;
use crate::output::world::history::HistoryDTO;
use crate::sim_core::world::integration::IntegrationAlgorithm;

pub struct BoxedWorldDTO {
  pub num_of_atoms: usize,
  pub size: Vector3<f64>,
  pub history: HistoryDTO,
  pub integration_algorithm: IntegrationAlgorithm,

  pub num_of_world_iterations: usize,
  pub number_of_resets: usize,
  pub max_iteration_till_reset: usize,

  pub frame_iteration_count: usize,
}