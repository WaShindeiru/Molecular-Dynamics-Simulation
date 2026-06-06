#[derive(Debug, Clone, Copy)]
pub struct CorrectionParam {
  pub small_distance: SmallDistanceCorrection,
}

impl Default for CorrectionParam {
  fn default() -> Self {
    Self {
      small_distance: SmallDistanceCorrection::default(),
    }
  }
}

#[derive(Debug, Clone, Copy)]
pub struct SmallDistanceCorrection {
  pub enabled: bool,
  pub substep_count: usize,
  /// Unitless length; trigger when `r_ij < distance_threshold`.
  pub distance_threshold: f64,
}

impl Default for SmallDistanceCorrection {
  fn default() -> Self {
    Self {
      enabled: false,
      substep_count: 1,
      distance_threshold: 0.0,
    }
  }
}
