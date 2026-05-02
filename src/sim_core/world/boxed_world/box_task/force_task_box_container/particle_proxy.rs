use std::sync::Arc;
use nalgebra::Vector3;
use crate::data::types::AtomType;
use crate::particle::Particle;
use crate::sim_core::world::boxed_world::box_container::sim_box::{get_coordinates_from_simulation_box_id, SimBoxEdge};
use crate::sim_core::world::boxed_world::integration::verlet_nose_hoover::computation::ForceComputationOperations;

#[derive(Copy, Clone)]
pub enum AxisPlacement {
  Left,
  Normal
}

#[derive(Copy, Clone)]
pub struct ParticlePlacement {
  pub x: AxisPlacement,
  pub y: AxisPlacement,
}

pub struct ParticlePositionProxy {
  particle: Arc<Particle>,
  particle_placement: ParticlePlacement,
  world_size: Vector3<f64>,
}

pub fn new_particle_position_proxy(particle: Arc<Particle>, particle_placement: ParticlePlacement,
                                   world_size: Vector3<f64>) -> ParticlePositionProxy {
  ParticlePositionProxy {
    particle,
    particle_placement,
    world_size,
  }
}

impl ForceComputationOperations for ParticlePositionProxy {
  fn get_id(&self) -> usize {
    self.particle.get_id()
  }

  fn get_position(&self) -> Vector3<f64> {
    let mut position = self.particle.get_position().clone();

    position.x = match self.particle_placement.x {
      AxisPlacement::Left => position.x + self.world_size.x,
      AxisPlacement::Normal => position.x,
    };

    position.y = match self.particle_placement.y {
      AxisPlacement::Left => position.y + self.world_size.y,
      AxisPlacement::Normal => position.y,
    };

    position
  }

  fn get_type(&self) -> AtomType {
    self.particle.get_type()
  }

  fn get_mass(&self) -> f64 {
    self.particle.get_mass()
  }

  fn prototype_clone(&self) -> Box<dyn ForceComputationOperations> {
    self.particle.prototype_clone()
  }
}