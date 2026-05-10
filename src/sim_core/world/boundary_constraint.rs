pub mod periodic;
pub mod simple;

#[derive(Debug, PartialEq, Copy, Clone)]
pub struct ParticleCompliance {
  pub compliant: bool,
  pub x: Compliance,
  pub y: Compliance,
  pub z: Compliance,
}

#[derive(Debug, PartialEq, Copy, Clone)]
pub enum Compliance {
  Compliant,
  ExceededLowerBoundary,
  ExceededHigherBoundary,
}

#[derive(Debug, PartialEq, Clone, Copy, serde::Serialize, serde::Deserialize)]
pub enum EdgeCondition {
  Simple,
  Periodic,
}
