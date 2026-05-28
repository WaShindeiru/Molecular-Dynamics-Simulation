use std::fmt;

use crate::data::ParticleConfig;
use crate::data::SimulationConfig;
use crate::simulations::generators::core::generator_config::GeneratorConfig;
use crate::simulations::generators::core::generator_config::dense::DenseGeneratorConfig;
use dense::DenseGenerator;
use nanotube::NanotubeGenerator;
use velocity_nanotube::VelocityNanotubeGenerator;

pub mod dense;
pub mod nanotube;
pub mod velocity_nanotube;

#[derive(Debug)]
pub struct GeneratorError(pub String);

impl fmt::Display for GeneratorError {
  fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
    write!(f, "GeneratorError: {}", self.0)
  }
}

impl std::error::Error for GeneratorError {}

impl From<GeneratorError> for std::io::Error {
  fn from(e: GeneratorError) -> Self {
    std::io::Error::new(std::io::ErrorKind::InvalidData, e)
  }
}

pub enum GeneratorType {
  Dense(DenseGenerator),
  Nanotube(NanotubeGenerator),
  VelocityNanotube(VelocityNanotubeGenerator),
}

pub trait Generator {
  fn generate(&self) -> Result<ParticleConfig, GeneratorError>;
}

impl Generator for GeneratorType {
  fn generate(&self) -> Result<ParticleConfig, GeneratorError> {
    match self {
      GeneratorType::Dense(g) => g.generate(),
      GeneratorType::Nanotube(g) => g.generate(),
      GeneratorType::VelocityNanotube(g) => g.generate(),
    }
  }
}

impl GeneratorType {
  pub fn dense(
    potential_gravity_max: f64,
    world_size: nalgebra::Vector3<f64>,
    particle_distance: f64,
    offset: nalgebra::Vector3<f64>,
  ) -> Self {
    Self::Dense(DenseGenerator::new(
      potential_gravity_max,
      world_size,
      particle_distance,
      offset,
    ))
  }

  pub fn dense_from_config(
    generator_config: &DenseGeneratorConfig,
    simulation_config: &SimulationConfig,
  ) -> Self {
    Self::Dense(generator_config.to_generator(simulation_config))
  }

  pub fn nanotube(
    particles: Vec<nanotube::NanotubeGeneratorParticle>,
    potential_gravity_max: f64,
    world_size: nalgebra::Vector3<f64>,
    offset: nalgebra::Vector3<f64>,
    vel_mean: nalgebra::Vector3<f64>,
    vel_std_dev: nalgebra::Vector3<f64>,
  ) -> Self {
    Self::Nanotube(NanotubeGenerator::new(
      particles,
      potential_gravity_max,
      world_size,
      offset,
      vel_mean,
      vel_std_dev,
    ))
  }

  pub fn from_config(
    generator_config: &GeneratorConfig,
    simulation_config: &SimulationConfig,
  ) -> Self {
    generator_config.to_generator(simulation_config)
  }
}
