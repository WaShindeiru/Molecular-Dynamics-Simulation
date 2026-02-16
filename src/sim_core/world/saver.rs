use std::io;
use crate::output::{SimpleWorldDTO, WorldDTO};

mod simple_world;
mod boxed_world;

#[derive(Clone)]
pub struct SaveOptions {
  pub save: bool,
  pub save_path: String,
  pub save_laamps: bool,
  pub save_verbose: bool,
}

pub struct PartialWorldSaver {
  save_options: SaveOptions,

  thermostat_work_total: f64,
  frame_iteration_count_current_iteration: usize, // for laamps purposes save only
}

impl PartialWorldSaver {
  pub fn new(save_options: SaveOptions) -> Self {
    PartialWorldSaver  {
      save_options,
      thermostat_work_total: 0.,
      frame_iteration_count_current_iteration: 0,
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