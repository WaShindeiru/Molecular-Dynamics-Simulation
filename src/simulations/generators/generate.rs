use crate::data::ParticleConfig;
use crate::data::SimulationConfig;
use crate::simulations::generators::generator_config::GeneratorConfig;
use crate::simulations::generators::generator_config::dense::DenseGeneratorConfig;
use dense::DenseGenerator;

pub mod dense;

pub enum GeneratorType {
  Dense(DenseGenerator),
}

pub trait Generator {
  fn generate(&self) -> ParticleConfig;
}

impl Generator for GeneratorType {
  fn generate(&self) -> ParticleConfig {
    match self {
      GeneratorType::Dense(g) => g.generate(),
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

  pub fn from_config(
    generator_config: &GeneratorConfig,
    simulation_config: &SimulationConfig,
  ) -> Self {
    generator_config.to_generator(simulation_config)
  }
}
