use std::fmt;

use super::types::IntegrationAlgorithm;

impl fmt::Display for IntegrationAlgorithm {
  fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
    match self {
      IntegrationAlgorithm::SemiImplicitEuler => write!(f, "SemiImplicitEuler"),
      IntegrationAlgorithm::VelocityVerlet => write!(f, "VelocityVerlet"),
      IntegrationAlgorithm::NoseHooverVerlet { .. } => {
        write!(f, "NoseHooverVerlet")
      }
    }
  }
}
