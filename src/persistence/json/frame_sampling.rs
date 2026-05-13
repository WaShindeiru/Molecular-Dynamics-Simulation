use crate::data::units::{TIME_U, ValueUnits};
use crate::sim_core::world::saver::FrameSamplingConfig;

#[derive(Clone, serde::Serialize, serde::Deserialize)]
pub struct FrameSamplingConfigFile {
  pub one_frame_duration: f64,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub frame_iteration_count: Option<usize>,
}

impl FrameSamplingConfigFile {
  pub fn from_runtime(config: &FrameSamplingConfig) -> Self {
    Self {
      one_frame_duration: config.one_frame_duration,
      frame_iteration_count: Some(config.frame_iteration_count),
    }
  }

  pub fn to_runtime(&self, time_step: f64, save_all_iterations: bool) -> FrameSamplingConfig {
    let frame_iteration_count = self.frame_iteration_count.unwrap_or_else(|| {
      if save_all_iterations {
        1
      } else {
        FrameSamplingConfig::iterations_from_duration(self.one_frame_duration, time_step)
      }
    });

    FrameSamplingConfig {
      one_frame_duration: self.one_frame_duration,
      frame_iteration_count,
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
