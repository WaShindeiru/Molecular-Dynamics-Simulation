use crate::data::units::{R_U, VELOCITY_U, ValueUnits};
use crate::persistence::json::particle_config::Vector3Record;
use crate::persistence::json::velocity_manager_file::VelocityChangeEntry;
use crate::simulations::generators::core::generate::nanotube::NanotubeGeneratorParticle;
use crate::simulations::generators::core::generate::velocity_nanotube::VelocityNanotubeGenerator;
use crate::simulations::generators::core::generator_config::nanotube::NanotubeGeneratorParticleFile;

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct VelocityNanotubeGeneratorConfig {
  pub particles: Vec<NanotubeGeneratorParticleFile>,
  pub offset: Vector3Record,
  pub velocities: Vec<VelocityChangeEntry>,
}

impl VelocityNanotubeGeneratorConfig {
  pub fn new(
    particles: Vec<NanotubeGeneratorParticleFile>,
    offset: Vector3Record,
    velocities: Vec<VelocityChangeEntry>,
  ) -> Self {
    Self {
      particles,
      offset,
      velocities,
    }
  }

  pub fn to_generator(&self) -> VelocityNanotubeGenerator {
    let particles = self
      .particles
      .iter()
      .map(|p| NanotubeGeneratorParticle::new(p.position.to_runtime(), p.particle_type))
      .collect();

    VelocityNanotubeGenerator::new(
      particles,
      self.offset.to_runtime(),
      self.velocities.iter().map(|v| v.to_runtime()).collect(),
    )
  }

  pub fn to_value_units(&self, source: ValueUnits, target: ValueUnits) -> Self {
    let pos_scale = ValueUnits::scale_between(source, target, R_U);
    let vel_scale = ValueUnits::scale_between(source, target, VELOCITY_U);

    Self {
      particles: self
        .particles
        .iter()
        .map(|p| NanotubeGeneratorParticleFile {
          position: p.position.scale(pos_scale),
          particle_type: p.particle_type,
        })
        .collect(),
      offset: self.offset.scale(pos_scale),
      velocities: self.velocities.iter().map(|v| v.scale(vel_scale)).collect(),
    }
  }
}
