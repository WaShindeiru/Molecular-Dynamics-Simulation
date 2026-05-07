use std::io;
use crate::output::world::WorldDTO;

mod simple_world;
mod boxed_world;

#[derive(Clone, serde::Serialize, serde::Deserialize)]
pub struct SaveOptions {
  pub save: bool,
  pub save_path: String,
  pub save_laamps: bool,
  pub save_verbose: bool,
}

pub struct PartialWorldSaver {
  save_options: SaveOptions,

  thermostat_work_total: f64,
  laamps_frame_iteration_count_current_iteration: usize,
  energy_frame_iteration_count_current_iteration: usize,
}

impl PartialWorldSaver {
  pub fn new(save_options: SaveOptions) -> Self {
    PartialWorldSaver  {
      save_options,
      thermostat_work_total: 0.,
      laamps_frame_iteration_count_current_iteration: 0,
      energy_frame_iteration_count_current_iteration: 0,
    }
  }

  pub fn persist(&mut self, world: &WorldDTO) -> io::Result<()> {
    if self.save_options.save {
      match world {
        WorldDTO::SimpleWorldDTO(simple_world) => self.persist_simple_world(simple_world),
        WorldDTO::BoxedWorldDTO(boxed_world) => self.persist_boxed_world(boxed_world),
      }
    } else {
      Ok(())
    }
  }
}
