use std::io;
use std::sync::mpsc::{Receiver, Sender};
use std::sync::{Arc, RwLock};
use std::thread::JoinHandle;
use log::info;
use nalgebra::Vector3;
use crate::data::InteractionType;
use crate::data::types::AtomType;
use crate::sim_core::world::boxed_world::box_container::BoxContainer;
use crate::output::WorldDTO;

use crate::sim_core::world::integration::{IntegrationAlgorithm, IntegrationAlgorithmParams};
use crate::sim_core::world::saver::SaveOptions;
use crate::particle::Particle;
use crate::sim_core::world::boxed_world::box_task::{BoxResult, BoxTask};
use box_task::threads::create_threads;

pub mod box_container;
pub mod cube;
pub mod box_task;
mod integration;

pub struct BoxedWorld {
  size: Vector3<f64>,
  box_container: Arc<RwLock<BoxContainer>>,
  num_of_atoms: usize,

  iteration: usize,
  
  max_iteration_till_reset: usize,
  reset_counter: usize,
  number_of_resets: usize,

  frame_iteration_count: usize,
  integration_algorithm: IntegrationAlgorithm,
  save_options: SaveOptions,

  threads: Vec<JoinHandle<()>>,
  tx_task: Sender<BoxTask>,
  rx_result: Receiver<BoxResult>,
}

impl BoxedWorld {
  pub fn new_from_atoms(
    atoms: Vec<Particle>,
    size: Vector3<f64>,
    max_iteration_till_reset: usize,
    frame_iteration_count: usize,
    integration_algorithm: IntegrationAlgorithm,
    save_options: SaveOptions,
  ) -> Self {

    let (tx_task, rx_result, threads) = create_threads();

    let box_type = {
      let mut fe = false;
      let mut c = false;

      for i in &atoms {
        if i.get_type() == AtomType::C {
          c = true;
        } else {
          fe = true;
        }
      }

      if fe && c {
        InteractionType::FeC
      } else if fe {
        InteractionType::FeFe
      } else {
        InteractionType::CC
      }
    };

    let temp = BoxContainer::new(atoms, size.clone(), box_type, max_iteration_till_reset);
    let box_container = Arc::new(RwLock::new(temp));

    // BoxedWorld {
    //   size,
    // 
    //   max_iteration_till_reset,
    //   reset_counter: 1,
    //   number_of_resets: 0,
    // 
    //   frame_iteration_count,
    //   integration_algorithm,
    //   save_options,
    // 
    //   threads,
    //   tx_task,
    //   rx_result,
    // 
    //   box_container,
    // }
    
    unimplemented!("aha");
  }

  pub fn save(&mut self) -> io::Result<()> {
    todo!("BoxedWorld::save not yet implemented")
  }

  pub fn update(&mut self, params: &IntegrationAlgorithmParams, time_step: f64, next_iteration: usize) {
    assert!(self.reset_counter <= self.max_iteration_till_reset);

    if self.reset_counter == self.max_iteration_till_reset {
      self.save().unwrap();
      self.reset_world();
    }

    match params {
      IntegrationAlgorithmParams::SemiImplicitEuler => unimplemented!("semi implicit euler for boxed world"),
      IntegrationAlgorithmParams::VelocityVerlet => unimplemented!("velocity verlet for boxed world!"),
      IntegrationAlgorithmParams::NoseHooverVerlet {
        desired_temperature,
        q_effective_mass } => {
          self.update_verlet_nose_hoover(time_step, next_iteration,
                                         *desired_temperature, *q_effective_mass);
      }
    }

    self.reset_counter += 1;
  }

  pub fn reset_world(&mut self) {
    info!("Resetting boxed world, reset number: {}", self.number_of_resets);

    let new_number_of_resets = self.number_of_resets + 1;

    self.box_container.write().unwrap().reset_container();

    self.number_of_resets = new_number_of_resets;
    self.reset_counter = 0;
  }

  pub fn get_size(&self) -> &Vector3<f64> {
    todo!("BoxedWorld::get_size not yet implemented")
  }

  pub fn to_transfer_struct(&self) -> WorldDTO {
    todo!("BoxedWorld::to_transfer_struct not yet implemented")
  }
}
