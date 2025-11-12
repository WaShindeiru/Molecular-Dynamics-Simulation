use crate::data::units::R_U;

pub struct AtomDTO {
  pub id: u64,
  pub atom_type: u64,
  pub x: f64,
  pub y: f64,
  pub z: f64,
}

pub struct WorldDTO {
  pub num_of_atoms: usize,
  pub atoms: Vec<Vec<AtomDTO>>,
  pub box_x: f64,
  pub box_y: f64,
  pub box_z: f64,
}

pub struct EngineDTO {
  pub num_of_iterations: usize,
  pub time_step: f64,
  pub world: WorldDTO,
}

pub fn change_length_unit(length: f64) -> f64 {
  length * R_U
}