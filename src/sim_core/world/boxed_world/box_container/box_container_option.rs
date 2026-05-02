use std::sync::Arc;

use crate::sim_core::world::boxed_world::box_container::sim_box::{SimulationBox, get_coordinates_from_simulation_box_id};
use crate::sim_core::world::boxed_world::box_container::BoxContainer;

impl BoxContainer<Option<Arc<SimulationBox>>> {
  pub fn get_box(&self, box_id: usize) -> Option<Arc<SimulationBox>> {
		let coordinates = get_coordinates_from_simulation_box_id(box_id, &self.config.box_count_dim);
		self.simulation_boxes.get(coordinates.x, coordinates.y, coordinates.z).unwrap().clone()
	}
}