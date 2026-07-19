use nalgebra::Vector3;

/// A particle id paired with a position already adjusted for periodic boundary conditions
/// (i.e. "fixed up" so pairwise distances to a reference cell can be computed directly,
/// without any wraparound logic at the call site).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FixedPositionParticle {
  pub id: usize,
  pub position: Vector3<f64>,
}
