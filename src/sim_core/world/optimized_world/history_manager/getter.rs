use crate::data::SimulationConfig;
use crate::sim_core::world::cell::box_container_config::BoxContainerConfig;
use crate::sim_core::world::optimized_world::history_manager::OptimizedHistoryManager;

impl OptimizedHistoryManager {
  pub fn simulation_config(&self) -> &SimulationConfig {
    &self.config
  }

  pub fn container_config(&self) -> &BoxContainerConfig {
    self.history.last().unwrap().config()
  }
}
