use crate::data::TimeIterationDistance;

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct GravityChangeEntry {
  pub distance: TimeIterationDistance,
  pub value: f64,
}

impl GravityChangeEntry {
  pub fn to_runtime(&self, time_step: f64) -> (usize, f64) {
    (self.distance.to_iteration(time_step), self.value)
  }

  pub fn from_runtime(iteration: usize, value: f64) -> Self {
    GravityChangeEntry {
      distance: TimeIterationDistance::Iteration { value: iteration },
      value,
    }
  }
}
