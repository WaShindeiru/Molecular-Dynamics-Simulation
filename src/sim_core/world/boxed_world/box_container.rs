use std::collections::HashMap;
use std::sync::Arc;
use crate::sim_core::world::boxed_world::box_container::box_container_config::BoxContainerConfig;
use crate::sim_core::world::boxed_world::box_container::sim_box::SimulationBox;
use crate::utils::cube::Cube;

pub mod sim_box;
pub mod box_container_config;
pub mod box_container_mutual;
pub mod box_container_arc;
pub mod box_container_value;
pub mod box_container_option;

pub struct BoxContainer<B = Arc<SimulationBox>> {
	config: BoxContainerConfig,
	simulation_boxes: Cube<B>,
	box_id_cache: HashMap<usize, usize>,
}
