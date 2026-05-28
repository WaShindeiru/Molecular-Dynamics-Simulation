use nalgebra::Vector3;
use rand_distr::{Distribution, Normal};

use crate::data::ParticleConfig;
use crate::data::types::AtomType;
use crate::particle::{Particle, SafeAtomFactory};
use crate::simulations::generators::core::generate::{Generator, GeneratorError};

#[derive(Debug, Clone)]
pub struct NanotubeGeneratorParticle {
  pub position: Vector3<f64>,
  pub particle_type: AtomType,
}

impl NanotubeGeneratorParticle {
  pub fn new(position: Vector3<f64>, particle_type: AtomType) -> Self {
    Self { position, particle_type }
  }
}

pub struct NanotubeGenerator {
  particles: Vec<NanotubeGeneratorParticle>,
  potential_gravity_max: f64,
  world_size: Vector3<f64>,
  offset: Vector3<f64>,
  vel_mean: Vector3<f64>,
  vel_std_dev: Vector3<f64>,
}

impl NanotubeGenerator {
  pub fn new(
    particles: Vec<NanotubeGeneratorParticle>,
    potential_gravity_max: f64,
    world_size: Vector3<f64>,
    offset: Vector3<f64>,
    vel_mean: Vector3<f64>,
    vel_std_dev: Vector3<f64>,
  ) -> Self {
    Self {
      particles,
      potential_gravity_max,
      world_size,
      offset,
      vel_mean,
      vel_std_dev,
    }
  }
}

impl Generator for NanotubeGenerator {
  fn generate(&self) -> Result<ParticleConfig, GeneratorError> {
    let atom_factory = SafeAtomFactory::new(self.potential_gravity_max, self.world_size.z);

    let vel_x = Normal::new(self.vel_mean.x, self.vel_std_dev.x)
      .map_err(|e| GeneratorError(format!("Invalid velocity distribution for x axis: {e}")))?;
    let vel_y = Normal::new(self.vel_mean.y, self.vel_std_dev.y)
      .map_err(|e| GeneratorError(format!("Invalid velocity distribution for y axis: {e}")))?;
    let vel_z = Normal::new(self.vel_mean.z, self.vel_std_dev.z)
      .map_err(|e| GeneratorError(format!("Invalid velocity distribution for z axis: {e}")))?;

    let mut rng = rand::rng();
    let mut atoms: Vec<Particle> = Vec::with_capacity(self.particles.len());

    for (i, p) in self.particles.iter().enumerate() {
      let position = p.position + self.offset;

      if position.x < 0.0
        || position.x > self.world_size.x
        || position.y < 0.0
        || position.y > self.world_size.y
        || position.z < 0.0
        || position.z > self.world_size.z
      {
        return Err(GeneratorError(format!(
          "Particle {i} at position ({}, {}, {}) is out of world bounds ({}, {}, {})",
          position.x, position.y, position.z,
          self.world_size.x, self.world_size.y, self.world_size.z,
        )));
      }

      let velocity = Vector3::new(
        vel_x.sample(&mut rng),
        vel_y.sample(&mut rng),
        vel_z.sample(&mut rng),
      );

      atoms.push(atom_factory.get_atom(p.particle_type, position, velocity));
    }

    Ok(ParticleConfig::new(atoms))
  }
}
