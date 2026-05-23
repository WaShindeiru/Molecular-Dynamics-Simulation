use crate::persistence::dto::atom::AtomDTO;
use crate::sim_core::world::thermostat::IntegrationAlgorithm;

pub struct SimpleWorldDTO {
  pub num_of_atoms: usize,
  pub atoms: Vec<Vec<AtomDTO>>,
  pub potential_energy: Vec<f64>,
  pub thermostat_epsilon: Vec<f64>,
  pub box_x: f64,
  pub box_y: f64,
  pub box_z: f64,
  pub integration_algorithm: IntegrationAlgorithm,

  pub num_of_world_iterations: usize,
  pub number_of_resets: usize,
  pub max_iteration_till_reset: usize,

  pub laamps_frame_iteration_count: usize,
  pub energy_frame_iteration_count: usize,
}
