use crate::data::SimulationConfig;
use crate::sim_core::world::boxed_world::box_container::box_container_config::BoxContainerConfig;
use crate::sim_core::world::boxed_world::history_manager::HistoryManager;

impl HistoryManager {
	pub fn simulation_config(&self) -> &SimulationConfig {
		&self.config
	}

	pub fn box_container_config(&self) -> &BoxContainerConfig {
		&self.box_container_config
	}
}