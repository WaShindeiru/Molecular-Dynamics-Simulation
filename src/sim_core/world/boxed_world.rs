use std::io;
use std::sync::mpsc::{Receiver, Sender};
use std::sync::{Arc, RwLock};
use std::thread::JoinHandle;
use log::info;
use nalgebra::Vector3;

use crate::data::InteractionType;
use crate::data::types::AtomType;
use crate::data::SimulationConfig;
use crate::sim_core::world::boxed_world::box_container::BoxContainer;
use crate::output::{BoxedWorldDTO, WorldDTO};
use crate::sim_core::world::integration::{new_integration_algorithm_state, IntegrationAlgorithm, IntegrationAlgorithmState};
use crate::sim_core::world::saver::{PartialWorldSaver, SaveOptions};
use crate::particle::Particle;
use crate::sim_core::world::boxed_world::box_task::{BoxResult, BoxTask};

use box_task::threads::create_threads;
use crate::sim_core::world::boundary_constraint::EdgeCondition;

pub mod box_container;
pub mod cube;
pub mod box_task;
pub mod integration;

pub struct BoxedWorld {
  config: SimulationConfig,
  box_container: Arc<RwLock<BoxContainer>>,

  iteration: usize,

  reset_counter: usize,
  number_of_resets: usize,

  frame_iteration_count: usize,
  integration_algorithm: IntegrationAlgorithm,
  integration_algorithm_state: IntegrationAlgorithmState,

  save_options: SaveOptions,
  world_saver: PartialWorldSaver,

  threads: Vec<JoinHandle<()>>,
  tx_task: Sender<BoxTask>,
  rx_result: Receiver<BoxResult>,
}

impl BoxedWorld {
  pub fn with_config(mut config: SimulationConfig) -> Self {
    let atoms = config.atoms.take().unwrap();
    let (tx_task, rx_result, threads) = create_threads(false);

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

    let frame_iteration_count = if !config.save_all_iterations {
      (config.one_frame_duration / config.time_step) as usize
    } else {
      1
    };

    BoxedWorld {
      config: config.clone(),
      box_container: Arc::new(RwLock::new(
        BoxContainer::with_config(config.clone(), Some(atoms))
        )
      ),

      iteration: 0,

      reset_counter: 1,
      number_of_resets: 0,

      frame_iteration_count,
      integration_algorithm_state: new_integration_algorithm_state(&config.integration_algorithm),
      integration_algorithm: config.integration_algorithm.clone(),

      save_options: config.save_options.clone(),
      world_saver: PartialWorldSaver::new(config.save_options.clone()),

      threads,
      tx_task,
      rx_result,
    }
  }

  // TODO: this is the same for both variants, maybe move world saver out of the world???
  pub fn save(&mut self) -> io::Result<()> {
    if self.save_options.save {
      let world = self.to_transfer_struct();
      self.world_saver.persist(&world)?;
    }
    Ok(())
  }

  pub fn update(&mut self, algorithm: &IntegrationAlgorithm, time_step: f64, next_iteration: usize) {
    assert!(self.reset_counter <= self.config.max_iteration_till_reset);

    if self.reset_counter == self.config.max_iteration_till_reset {
      self.save().unwrap();
      self.reset_world();
    }

    match algorithm {
      IntegrationAlgorithm::SemiImplicitEuler => unimplemented!("semi implicit euler for boxed world"),
      IntegrationAlgorithm::VelocityVerlet => unimplemented!("velocity verlet for boxed world!"),
      IntegrationAlgorithm::NoseHooverVerlet { .. } => {
          self.update_verlet_nose_hoover(time_step, next_iteration);
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
    &self.config.world_size
  }

  pub fn to_transfer_struct(&self) -> WorldDTO {
    let lower_index: usize;
    if self.number_of_resets > 0 {
      lower_index = 1;
    } else {
      lower_index = 0;
    }

    WorldDTO::BoxedWorldDTO(
      BoxedWorldDTO {
        num_of_atoms: self.config.num_of_atoms,
        size: self.config.world_size,
        box_container: self.box_container.read().unwrap().to_transfer_struct(lower_index),
        integration_algorithm: self.integration_algorithm.clone(),

        num_of_world_iterations: self.iteration,
        number_of_resets: self.number_of_resets,
        max_iteration_till_reset: self.config.max_iteration_till_reset,

        frame_iteration_count: self.frame_iteration_count,
      }
    )
  }

  pub fn get_particle_counts(&self) -> (usize, usize, usize) {
    self.box_container.read().unwrap().get_particle_counts()
  }
}
