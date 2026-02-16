use nalgebra::Vector3;
use crate::data::InteractionType;
use crate::data::units::R_U;
use crate::sim_core::world::integration::IntegrationAlgorithm;

#[derive(Clone)]
pub struct AtomDTO {
  pub id: u64,
  pub iteration: usize,
  pub atom_type: u64,
  pub x: f64,
  pub y: f64,
  pub z: f64,
  
  pub kinetic_energy: f64,
  pub potential_energy: f64,
  pub thermostat_work: f64,

  pub force_x: f64,
  pub force_y: f64,
  pub force_z: f64,
}

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

  pub frame_iteration_count: usize,
}

pub struct BoxContainerDTO {
  pub box_type: InteractionType,
  pub atoms: Vec<Vec<AtomDTO>>,
  pub thermostat_epsilon: Vec<f64>,

  pub box_length: Vector3<f64>,
  pub box_count: usize,
  pub box_count_dim: Vector3<usize>,
}

pub struct BoxedWorldDTO {
  pub num_of_atoms: usize,
  pub size: Vector3<f64>,
  pub box_container: BoxContainerDTO,
  pub integration_algorithm: IntegrationAlgorithm,

  pub num_of_world_iterations: usize,
  pub number_of_resets: usize,
  pub max_iteration_till_reset: usize,

  pub frame_iteration_count: usize,
}

pub enum WorldDTO {
  SimpleWorldDTO(SimpleWorldDTO),
  BoxedWorldDTO(BoxedWorldDTO),
}

pub struct EngineDTO {
  pub num_of_iterations: usize,
  pub time_step: f64,
  pub world: WorldDTO,
}

pub fn change_length_unit(length: f64) -> f64 {
  length * R_U
}