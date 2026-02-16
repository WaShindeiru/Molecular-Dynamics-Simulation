use std::fmt;

#[derive(Debug, Clone, Copy)]
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

#[derive(Debug)]
pub enum IntegrationAlgorithmParams {
  SemiImplicitEuler,
  VelocityVerlet,
  NoseHooverVerlet {
    desired_temperature: f64,
    q_effective_mass: f64
  },
}

pub fn validate_integration_params(algorithm: &IntegrationAlgorithm, params: &IntegrationAlgorithmParams) -> bool {
  match (algorithm, params) {
    (IntegrationAlgorithm::SemiImplicitEuler, IntegrationAlgorithmParams::SemiImplicitEuler) => true,
    (IntegrationAlgorithm::VelocityVerlet, IntegrationAlgorithmParams::VelocityVerlet) => true,
    (IntegrationAlgorithm::NoseHooverVerlet, IntegrationAlgorithmParams::NoseHooverVerlet { .. }) => true,
    _ => false,
  }
}
