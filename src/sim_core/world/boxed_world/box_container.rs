use crate::sim_core::world::boxed_world::box_container::box_container_config::BoxContainerConfig;
use crate::sim_core::world::boxed_world::box_container::sim_box::SimulationBox;
use crate::utils::cube::Cube;
use std::collections::HashMap;
use std::sync::Arc;

pub mod box_container_arc;
pub mod box_container_config;
pub mod box_container_mutual;
pub mod box_container_option;
pub mod box_container_value;
pub mod periodic;
pub mod sim_box;

pub use periodic::{face_neighbor_box_ids_periodic, get_needed_box_id_periodic};

pub struct BoxContainer<B = Arc<SimulationBox>> {
  config: BoxContainerConfig,
  simulation_boxes: Cube<B>,
  box_id_cache: HashMap<usize, usize>,
}
