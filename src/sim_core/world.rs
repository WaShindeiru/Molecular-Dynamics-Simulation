use std::io;
use nalgebra::Vector3;
use crate::output::WorldDTO;
use crate::particle::Particle;
use crate::sim_core::world::boxed_world::BoxedWorld;
use crate::sim_core::world::simple_world::SimpleWorld;
use crate::sim_core::world::integration::{IntegrationAlgorithm, IntegrationAlgorithmParams};
use crate::sim_core::world::saver::SaveOptions;

pub mod integration;
pub mod simple_world;
pub mod boxed_world;
pub mod saver;

fn get_index_for_iteration(current_iteration: usize, max_iteration_till_reset: usize, reset_count: usize) -> usize {
  let mut index = current_iteration;
  if reset_count > 0 {
    index = index - (max_iteration_till_reset - 1);
  }
  index = index - (max_iteration_till_reset) * (reset_count - 1);
  index
}

pub enum WorldType {
  SimpleWorld,
  BoxedWorld,
}

pub enum World {
  SimpleWorld(SimpleWorld),
  BoxedWorld(BoxedWorld),
}

impl World {
  pub fn new_from_atoms(
    atoms: Vec<Particle>,
    size: Vector3<f64>,
    max_iteration_till_reset: usize,
    frame_iteration_count: usize,
    integration_algorithm: IntegrationAlgorithm,
    save_options: SaveOptions,
    world_type: WorldType
  ) -> Self {
    match world_type {
      WorldType::BoxedWorld => {
        World::BoxedWorld(BoxedWorld::new_from_atoms(
          atoms,
          size,
          max_iteration_till_reset,
          frame_iteration_count,
          integration_algorithm,
          save_options
        ))
      }
      WorldType::SimpleWorld => {
        World::SimpleWorld(SimpleWorld::new_from_atoms(
          atoms,
          size,
          max_iteration_till_reset,
          frame_iteration_count,
          integration_algorithm,
          save_options
        ))
      }
    }
  }

  pub fn save(&mut self) -> io::Result<()> {
    match self {
      World::SimpleWorld(world) => world.save(),
      World::BoxedWorld(world) => world.save(),
    }
  }

  pub fn update(&mut self, params: &IntegrationAlgorithmParams, time_step: f64, next_iteration: usize) {
    match self {
      World::SimpleWorld(world) => world.update(params, time_step, next_iteration),
      World::BoxedWorld(world) => world.update(params, time_step, next_iteration),
    }
  }

  pub fn reset_world(&mut self) {
    match self {
      World::SimpleWorld(world) => world.reset_world(),
      World::BoxedWorld(world) => world.reset_world(),
    }
  }

  pub fn get_size(&self) -> &Vector3<f64> {
    match self {
      World::SimpleWorld(world) => world.get_size(),
      World::BoxedWorld(world) => world.get_size(),
    }
  }

  pub fn apply_boundary_constraint(&self, mut atom: Particle) -> Particle {
    let size = self.get_size();
    
    atom = match atom.get_position().x {
      x if x < 0.0 => {
        atom.set_velocity(Vector3::new(-atom.get_velocity().x, atom.get_velocity().y, atom.get_velocity().z));
        atom.update_position(Vector3::new(-atom.get_position().x, atom.get_position().y, atom.get_position().z));
        atom
      }
      x if x > size.x => {
        atom.set_velocity(Vector3::new(-atom.get_velocity().x, atom.get_velocity().y, atom.get_velocity().z));
        atom.update_position(Vector3::new(2. * size.x - atom.get_position().x, atom.get_position().y, atom.get_position().z));
        atom
      },
      _ => atom
    };

    atom = match atom.get_position().y {
      y if y < 0.0 => {
        atom.set_velocity(Vector3::new(atom.get_velocity().x, -atom.get_velocity().y, atom.get_velocity().z));
        atom.update_position(Vector3::new(atom.get_position().x, -atom.get_position().y, atom.get_position().z));
        atom
      }
      y if y > size.y => {
        atom.set_velocity(Vector3::new(atom.get_velocity().x, -atom.get_velocity().y, atom.get_velocity().z));
        atom.update_position(Vector3::new(atom.get_position().x, 2. * size.y - atom.get_position().y, atom.get_position().z));
        atom
      },
      _ => atom
    };

    atom = match atom.get_position().z {
      z if z < 0.0 => {
        atom.set_velocity(Vector3::new(atom.get_velocity().x, atom.get_velocity().y, -atom.get_velocity().z));
        atom.update_position(Vector3::new(atom.get_position().x, atom.get_position().y, -atom.get_position().z));
        atom
      }
      z if z > size.z => {
        atom.set_velocity(Vector3::new(atom.get_velocity().x, atom.get_velocity().y, -atom.get_velocity().z));
        atom.update_position(Vector3::new(atom.get_position().x, atom.get_position().y, 2. * size.z - atom.get_position().z));
        atom
      },
      _ => atom
    };

    atom
  }

  pub fn to_transfer_struct(&self) -> WorldDTO {
    match self {
      World::SimpleWorld(world) => world.to_transfer_struct(),
      World::BoxedWorld(world) => world.to_transfer_struct(),
    }
  }
}
