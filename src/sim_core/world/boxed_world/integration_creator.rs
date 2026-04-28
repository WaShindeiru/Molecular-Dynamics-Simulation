use crate::sim_core::world::boxed_world::history_manager::sim_box::SimulationBox;
use crate::sim_core::world::boxed_world::cube::Cube;

pub struct IntegrationCreator {
  integration_box_modified: Cube<SimulationBox>,
}