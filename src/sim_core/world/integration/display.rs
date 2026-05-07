use std::fmt;

use super::types::IntegrationAlgorithm;

impl fmt::Display for IntegrationAlgorithm {
  fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
    match self {
      IntegrationAlgorithm::SemiImplicitEuler => write!(f, "SemiImplicitEuler"),
      IntegrationAlgorithm::VelocityVerlet => write!(f, "VelocityVerlet"),
      IntegrationAlgorithm::NoseHooverVerlet {
        desired_temperature: _t,
        q_effective_mass: _q,
      } => {
        write!(f, "NoseHooverVerlet")
      },
    }
  }
}
