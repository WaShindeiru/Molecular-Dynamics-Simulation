use std::fmt;

pub mod verlet;
pub mod euler;
pub mod verlet_nose_hoover;

#[derive(Clone)]
pub enum IntegrationAlgorithm {
  SemiImplicitEuler,
  VelocityVerlet,
  NoseHooverVerlet,
}

impl fmt::Display for IntegrationAlgorithm {
  fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
    match self {
      IntegrationAlgorithm::SemiImplicitEuler => write!(f, "SemiImplicitEuler"),
      IntegrationAlgorithm::VelocityVerlet => write!(f, "VelocityVerlet"),
      IntegrationAlgorithm::NoseHooverVerlet => write!(f, "NoseHooverVerlet"),
    }
  }
}

pub enum IntegrationAlgorithmParams {
  SemiImplicitEuler,
  VelocityVerlet,
  NoseHooverVerlet {
    desired_temperature: f64,
    q_effective_mass: f64
  },
}

