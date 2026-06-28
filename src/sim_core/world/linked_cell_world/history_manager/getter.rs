use crate::data::SimulationConfig;
use crate::sim_core::world::boxed_world::box_container::box_container_config::BoxContainerConfig;
use crate::sim_core::world::linked_cell_world::history_manager::LinkedCellHistoryManager;

impl LinkedCellHistoryManager {
  pub fn simulation_config(&self) -> &SimulationConfig {
    &self.config
  }

  pub fn container_config(&self) -> &BoxContainerConfig {
    self.history.last().unwrap().config()
  }
}
