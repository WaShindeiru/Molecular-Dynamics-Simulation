use std::collections::HashMap;
use std::ops::Deref;
use std::sync::Arc;
use crate::particle::Particle;
use crate::sim_core::world::boxed_world::box_container::box_container_config::BoxContainerConfig;
use crate::sim_core::world::boxed_world::box_container::BoxContainer;
use crate::sim_core::world::boxed_world::box_container::sim_box::{get_coordinates_from_simulation_box_id, SimulationBox};

impl<B> BoxContainer<B> {
	pub fn assign_box_id_for_particle(&self, particle: &Particle) -> usize {
		self.config.box_id_for_position(particle.get_position())
	}

	pub fn config(&self) -> &BoxContainerConfig {
		&self.config
	}

	pub fn particle_box_id(&self, particle_id: usize) -> usize {
		*self.box_id_cache.get(&particle_id).expect("Particle not found in box_id_cache")
	}
}

impl<B: Deref<Target = SimulationBox>> BoxContainer<B> {
	pub fn particles_of_box<'a>(&'a self, box_id: usize) -> impl Iterator<Item = Arc<Particle>> + 'a {
		let coordinates = get_coordinates_from_simulation_box_id(box_id, &self.config.box_count_dim);
		self.simulation_boxes
			.get(coordinates.x, coordinates.y, coordinates.z)
			.unwrap()
			.particles()
			.values()
			.cloned()
	}

	pub fn particles_of_boxes<'a>(&'a self, box_ids: &'a [usize]) -> impl Iterator<Item = Arc<Particle>> + 'a {
		box_ids.iter().flat_map(move |box_id_ref| {
			let box_id: usize = *box_id_ref;
			self.particles_of_box(box_id)
		})
	}

	pub fn all_particles_reset(&self) -> HashMap<usize, Particle> {
		self.simulation_boxes
			.iter()
			.flat_map(|sim_box| sim_box.particles().iter())
			.map(|(&id, particle)| (id, particle.reset_clone()))
			.collect()
	}

	pub fn atom_count_of_boxes(&self, box_ids: &[usize]) -> usize {
		box_ids.iter().map(|&box_id| {
			let coordinates = get_coordinates_from_simulation_box_id(box_id, &self.config.box_count_dim);
			self.simulation_boxes
				.get(coordinates.x, coordinates.y, coordinates.z)
				.unwrap()
				.len()
		}).sum()
	}
}