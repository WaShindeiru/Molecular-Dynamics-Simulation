use crate::output::world::WorldDTO;

pub struct EngineDTO {
  pub num_of_iterations: usize,
  pub time_step: f64,
  pub world: WorldDTO,
}