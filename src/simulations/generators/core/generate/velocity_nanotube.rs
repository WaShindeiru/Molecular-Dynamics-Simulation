use nalgebra::Vector3;

use crate::data::config::particle_config::VelocityScheduleConfig;
use crate::data::ParticleConfig;
use crate::particle::{Particle, SafeAtomFactory};
use crate::simulations::generators::core::generate::nanotube::NanotubeGeneratorParticle;
use crate::simulations::generators::core::generate::{Generator, GeneratorError};

pub struct VelocityNanotubeGenerator {
  particles: Vec<NanotubeGeneratorParticle>,
  offset: Vector3<f64>,
  velocities: Vec<(usize, Vector3<f64>)>,
}

impl VelocityNanotubeGenerator {
  pub fn new(
    particles: Vec<NanotubeGeneratorParticle>,
    offset: Vector3<f64>,
    velocities: Vec<(usize, Vector3<f64>)>,
  ) -> Self {
    Self {
      particles,
      offset,
      velocities,
    }
  }
}

impl Generator for VelocityNanotubeGenerator {
  fn generate(&self) -> Result<ParticleConfig, GeneratorError> {
    if self.velocities.iter().all(|(iteration, _)| *iteration != 0) {
      return Err(GeneratorError(
        "VelocityNanotubeGenerator requires velocity entry for iteration 0".to_string(),
      ));
    }

    // All particles from this generator share the same velocity schedule, so they
    // share a single velocity manager instead of one duplicate manager each.
    let particle_velocity_manager_id = 0;
    let mut atoms: Vec<Particle> = Vec::with_capacity(self.particles.len());
    let atom_factory = SafeAtomFactory::new(0.0, 1.0);

    for p in self.particles.iter() {
      let position = p.position + self.offset;
      let particle =
        atom_factory.get_atom_custom_velocity(p.particle_type, position, particle_velocity_manager_id);
      atoms.push(particle);
    }

    let velocity_schedules = vec![VelocityScheduleConfig {
      particle_velocity_manager_id,
      velocities: self.velocities.clone(),
    }];

    Ok(ParticleConfig::new_with_schedules(atoms, velocity_schedules))
  }
}