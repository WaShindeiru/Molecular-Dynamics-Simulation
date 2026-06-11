use crate::data::TimeIterationDistance;
use crate::data::units::TIME_U;
use crate::data::units::ValueUnits;

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct TimestepChangeEntry {
  pub distance: TimeIterationDistance,
  pub value: f64,
}

impl TimestepChangeEntry {
  pub fn to_runtime(&self, initial_timestep: f64) -> (usize, f64) {
    (self.distance.to_iteration(initial_timestep), self.value)
  }

  pub fn from_runtime(iteration: usize, value: f64) -> Self {
    TimestepChangeEntry {
      distance: TimeIterationDistance::Iteration { value: iteration },
      value,
    }
  }

  pub fn to_value_units(&self, source: ValueUnits, target: ValueUnits) -> Self {
    TimestepChangeEntry {
      distance: self.distance.to_value_units(source, target),
      value: self.value * ValueUnits::scale_between(source, target, TIME_U),
    }
  }
}
