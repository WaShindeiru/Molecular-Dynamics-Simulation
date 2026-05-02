pub mod atom;
pub mod world;
pub mod engine;

use nalgebra::Vector3;
use crate::data::units::R_U;
use crate::sim_core::world::boxed_world::box_container::box_container_config::BoxContainerConfig;
use crate::sim_core::world::integration::IntegrationAlgorithm;

pub fn change_length_unit(length: f64) -> f64 {
  length * R_U
}