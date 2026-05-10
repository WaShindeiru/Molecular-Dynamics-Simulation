use crate::data::units::{TIME_U, ValueUnits};
use crate::sim_core::world::saver::FrameSamplingConfig;

#[derive(Clone, serde::Serialize, serde::Deserialize)]
pub struct FrameSamplingConfigFile {
  pub one_frame_duration: f64,
  pub frame_iteration_count: usize,
}

impl FrameSamplingConfigFile {
  pub fn from_runtime(config: &FrameSamplingConfig) -> Self {
    Self {
      one_frame_duration: config.one_frame_duration,
      frame_iteration_count: config.frame_iteration_count,
    }
  }

  pub fn to_runtime(&self) -> FrameSamplingConfig {
    FrameSamplingConfig {
      one_frame_duration: self.one_frame_duration,
      frame_iteration_count: self.frame_iteration_count,
    }
  }

  pub fn to_value_units(&self, source: ValueUnits, target: ValueUnits) -> Self {
    Self {
      one_frame_duration: self.one_frame_duration
        * ValueUnits::scale_between(source, target, TIME_U),
      frame_iteration_count: self.frame_iteration_count,
    }
  }
}
