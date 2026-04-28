use std::collections::HashMap;
use std::sync::Arc;
use crate::particle::Particle;
use crate::sim_core::world::boxed_world::box_container::sim_box::get_coordinates_from_simulation_box_id;
use crate::sim_core::world::boxed_world::history_manager::HistoryManager;

impl HistoryManager {
	pub fn set_integration_box_cache(&mut self, cache: HashMap<usize, Particle>) {
		let mut box_id_cache: HashMap<usize, usize> = HashMap::with_capacity(cache.len());
		let mut new_boxes = self.boxes.last().unwrap().clone();

		for sim_box in new_boxes.iter_mut() {
			sim_box.clear_box();
		}

		for (i_id, particle_i) in cache {
			assert_eq!(i_id, particle_i.get_id() as usize);
			let box_id = self.assign_box_for_particle(&particle_i);
			box_id_cache.insert(particle_i.get_id() as usize, box_id);
			let coordinates = get_coordinates_from_simulation_box_id(box_id,
			                                                         &self.box_count_dim);
			new_boxes.get_mut(coordinates.x, coordinates.y, coordinates.z).unwrap()
				.add_particle(Arc::new(particle_i));
		}

		self.integration_box_cache = new_boxes;
		self.integration_box_id_cache = box_id_cache;
	}
}