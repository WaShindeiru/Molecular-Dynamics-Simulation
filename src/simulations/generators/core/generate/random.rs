use nalgebra::Vector3;
use rand::distr::{Distribution, Uniform};
use rand_distr::Normal;

use crate::data::ParticleConfig;
use crate::data::types::AtomType::{C, Fe};
use crate::particle::{Particle, SafeAtomFactory};
use crate::simulations::generators::core::generate::{Generator, GeneratorError};

pub struct RandomGenerator {
  num_particles: usize,
  potential_gravity_max: f64,
  world_size: Vector3<f64>,
  vel_mean: Vector3<f64>,
  vel_std_dev: Vector3<f64>,
}

impl RandomGenerator {
  pub fn new(
    num_particles: usize,
    potential_gravity_max: f64,
    world_size: Vector3<f64>,
    vel_mean: Vector3<f64>,
    vel_std_dev: Vector3<f64>,
  ) -> Self {
    Self {
      num_particles,
      potential_gravity_max,
      world_size,
      vel_mean,
      vel_std_dev,
    }
  }
}

impl Generator for RandomGenerator {
  fn generate(&self) -> Result<ParticleConfig, GeneratorError> {
    if self.num_particles == 0 {
      return Err(GeneratorError("num_particles must be greater than 0".to_string()));
    }
    if self.world_size.x <= 0.0 || self.world_size.y <= 0.0 || self.world_size.z <= 0.0 {
      return Err(GeneratorError(format!(
        "world_size components must be positive, got ({}, {}, {})",
        self.world_size.x, self.world_size.y, self.world_size.z,
      )));
    }

    let atom_factory = SafeAtomFactory::new(self.potential_gravity_max, self.world_size.z);

    let pos_x = Uniform::new(0.0, self.world_size.x)
      .map_err(|e| GeneratorError(format!("Invalid position distribution for x axis: {e}")))?;
    let pos_y = Uniform::new(0.0, self.world_size.y)
      .map_err(|e| GeneratorError(format!("Invalid position distribution for y axis: {e}")))?;
    let pos_z = Uniform::new(0.0, self.world_size.z)
      .map_err(|e| GeneratorError(format!("Invalid position distribution for z axis: {e}")))?;

    let type_range = Uniform::new(0.0, 1.0)
      .map_err(|e| GeneratorError(format!("Invalid atom-type distribution: {e}")))?;

    let vel_x = Normal::new(self.vel_mean.x, self.vel_std_dev.x)
      .map_err(|e| GeneratorError(format!("Invalid velocity distribution for x axis: {e}")))?;
    let vel_y = Normal::new(self.vel_mean.y, self.vel_std_dev.y)
      .map_err(|e| GeneratorError(format!("Invalid velocity distribution for y axis: {e}")))?;
    let vel_z = Normal::new(self.vel_mean.z, self.vel_std_dev.z)
      .map_err(|e| GeneratorError(format!("Invalid velocity distribution for z axis: {e}")))?;

    let mut rng = rand::rng();
    let mut atoms: Vec<Particle> = Vec::with_capacity(self.num_particles);

    for _ in 0..self.num_particles {
      let position = Vector3::new(
        pos_x.sample(&mut rng),
        pos_y.sample(&mut rng),
        pos_z.sample(&mut rng),
      );
      let velocity = Vector3::new(
        vel_x.sample(&mut rng),
        vel_y.sample(&mut rng),
        vel_z.sample(&mut rng),
      );
      let atom_type = if type_range.sample(&mut rng) < 0.75 { Fe } else { C };
      atoms.push(atom_factory.get_atom(atom_type, position, velocity));
    }

    Ok(ParticleConfig::new(atoms))
  }
}
