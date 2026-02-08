use std::io;
use nalgebra::Vector3;
use crate::output::WorldDTO;
use crate::particle::Particle;
use crate::sim_core::world::integration::{IntegrationAlgorithm, IntegrationAlgorithmParams};
use crate::sim_core::world::saver::SaveOptions;

pub struct BoxedWorld {
  
}

impl BoxedWorld {
  pub fn new_from_atoms(
    atoms: Vec<Particle>,
    size: Vector3<f64>,
    max_iteration_till_reset: usize,
    frame_iteration_count: usize,
    integration_algorithm: IntegrationAlgorithm,
    save_options: SaveOptions
  ) -> Self {
    todo!("BoxedWorld::new_from_atoms not yet implemented")
  }

  pub fn save(&mut self) -> io::Result<()> {
    todo!("BoxedWorld::save not yet implemented")
  }

  pub fn update(&mut self, params: &IntegrationAlgorithmParams, time_step: f64, next_iteration: usize) {
    todo!("BoxedWorld::update not yet implemented")
  }

  pub fn reset_world(&mut self) {
    todo!("BoxedWorld::reset_world not yet implemented")
  }

  pub fn get_size(&self) -> &Vector3<f64> {
    todo!("BoxedWorld::get_size not yet implemented")
  }

  pub fn to_transfer_struct(&self) -> WorldDTO {
    todo!("BoxedWorld::to_transfer_struct not yet implemented")
  }
}
