use std::collections::HashMap;
use std::sync::Arc;
use crate::particle::Particle;
use crate::sim_core::world::boxed_world::box_container::box_container_config::BoxContainerConfig;
use crate::sim_core::world::boxed_world::box_container::BoxContainer;
use crate::sim_core::world::boxed_world::box_container::sim_box::{SimulationBox, get_coordinates_from_simulation_box_id, get_id_simulation_box};
use crate::utils::cube::Cube;

impl BoxContainer<SimulationBox> {
	pub fn new_local(config: BoxContainerConfig) -> Self {
		let boxes: Cube<SimulationBox> = Cube::new(
			config.box_count_dim.x,
			config.box_count_dim.y,
			config.box_count_dim.z,
		);

		BoxContainer {
			config,
			simulation_boxes: boxes,
			box_id_cache: HashMap::new(),
		}
	}

	pub fn get_box_mut(&mut self, box_id: usize) -> &mut SimulationBox {
		let coordinates = get_coordinates_from_simulation_box_id(box_id, &self.config.box_count_dim);
		self.simulation_boxes.get_mut(coordinates.x, coordinates.y, coordinates.z).unwrap()
	}

	pub fn add_particle(&mut self, particle: Arc<Particle>) {
		let particle_id = particle.get_id();
		let box_coordinates = self.config().box_coordinates_for_position(particle.get_position());
		self.simulation_boxes.get_mut(box_coordinates.x, box_coordinates.y, box_coordinates.z).unwrap().add_particle(particle);

		let box_id = get_id_simulation_box(&box_coordinates, &self.config().box_count_dim);
		self.box_id_cache.insert(particle_id, box_id);
	}

	pub fn get_box(&self, box_id: usize) -> &SimulationBox {
		let coordinates = get_coordinates_from_simulation_box_id(box_id, &self.config.box_count_dim);
		self.simulation_boxes.get(coordinates.x, coordinates.y, coordinates.z).unwrap()
	}

	pub fn simulation_boxes(&self) -> &Cube<SimulationBox> {
		&self.simulation_boxes
	}

	pub fn simulation_boxes_mut(&mut self) -> &mut Cube<SimulationBox> {
		&mut self.simulation_boxes
	}
	
	pub fn into_shared(self) -> BoxContainer<Arc<SimulationBox>> {
		let (x, y, z) = self.simulation_boxes.dimensions();
		let mut shared_boxes: Cube<Arc<SimulationBox>> = Cube::new(x, y, z);

		for ((xi, yi, zi), sim_box) in self.simulation_boxes.iter_with_coords() {
			shared_boxes.set(xi, yi, zi, Arc::new(sim_box.clone())).unwrap();
		}

		BoxContainer {
			config: self.config,
			simulation_boxes: shared_boxes,
			box_id_cache: self.box_id_cache,
		}
	}
}