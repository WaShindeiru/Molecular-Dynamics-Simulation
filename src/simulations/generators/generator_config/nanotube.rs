use crate::data::SimulationConfig;
use crate::data::types::AtomType;
use crate::data::units::{R_U, VELOCITY_U, ValueUnits};
use crate::persistence::json::particle_config::Vector3Record;
use crate::simulations::generators::generate::nanotube::{NanotubeGenerator, NanotubeGeneratorParticle};

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct NanotubeGeneratorParticleFile {
  pub position: Vector3Record,
  pub particle_type: AtomType,
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct NanotubeGeneratorConfig {
  pub particles: Vec<NanotubeGeneratorParticleFile>,
  pub offset: Vector3Record,
  pub vel_mean: Vector3Record,
  pub vel_std_dev: Vector3Record,
}

impl NanotubeGeneratorConfig {
  pub fn new(
    particles: Vec<NanotubeGeneratorParticleFile>,
    offset: Vector3Record,
    vel_mean: Vector3Record,
    vel_std_dev: Vector3Record,
  ) -> Self {
    Self { particles, offset, vel_mean, vel_std_dev }
  }

  pub fn to_generator(&self, simulation_config: &SimulationConfig) -> NanotubeGenerator {
    let particles = self.particles.iter().map(|p| {
      NanotubeGeneratorParticle::new(p.position.to_runtime(), p.particle_type)
    }).collect();

    NanotubeGenerator::new(
      particles,
      simulation_config.potential_gravity_max,
      simulation_config.world_size,
      self.offset.to_runtime(),
      self.vel_mean.to_runtime(),
      self.vel_std_dev.to_runtime(),
    )
  }

  pub fn to_value_units(&self, source: ValueUnits, target: ValueUnits) -> Self {
    let pos_scale = ValueUnits::scale_between(source, target, R_U);
    let vel_scale = ValueUnits::scale_between(source, target, VELOCITY_U);

    Self {
      particles: self.particles.iter().map(|p| NanotubeGeneratorParticleFile {
        position: p.position.scale(pos_scale),
        particle_type: p.particle_type,
      }).collect(),
      offset: self.offset.scale(pos_scale),
      vel_mean: self.vel_mean.scale(vel_scale),
      vel_std_dev: self.vel_std_dev.scale(vel_scale),
    }
  }
}
