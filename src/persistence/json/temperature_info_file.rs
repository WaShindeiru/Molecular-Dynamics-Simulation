use crate::data::TimeIterationDistance;
use crate::data::units::{TEMPERATURE_U, ValueUnits};
use crate::sim_core::world::thermostat::{
  DEFAULT_ACCEPTANCE_TIME_UNITLESS, DEFAULT_TEMP_THRESHOLD_UNITLESS, TemperatureInfo,
};

fn default_acceptance_distance() -> TimeIterationDistance {
  TimeIterationDistance::Time {
    value: DEFAULT_ACCEPTANCE_TIME_UNITLESS,
  }
}

fn default_achieved_distance() -> TimeIterationDistance {
  TimeIterationDistance::Iteration { value: 0 }
}

fn default_temperature_threshold() -> f64 {
  DEFAULT_TEMP_THRESHOLD_UNITLESS
}

#[derive(Debug, Clone, Copy, serde::Serialize, serde::Deserialize)]
pub struct TemperatureInfoFile {
  pub desired_temperature: f64,
  #[serde(default = "default_acceptance_distance")]
  pub acceptance_distance: TimeIterationDistance,
  #[serde(alias = "distance", default = "default_achieved_distance")]
  pub achieved_distance: TimeIterationDistance,
  #[serde(default = "default_temperature_threshold")]
  pub threshold: f64,
  pub save: bool,
}

impl TemperatureInfoFile {
  pub fn from_runtime(info: &TemperatureInfo) -> Self {
    TemperatureInfoFile {
      desired_temperature: info.desired_temperature,
      acceptance_distance: info.acceptance_distance,
      achieved_distance: info.achieved_distance,
      threshold: info.threshold,
      save: info.save,
    }
  }

  pub fn to_runtime(&self) -> TemperatureInfo {
    TemperatureInfo {
      desired_temperature: self.desired_temperature,
      acceptance_distance: self.acceptance_distance,
      achieved_distance: self.achieved_distance,
      threshold: self.threshold,
      save: self.save,
    }
  }

  pub fn to_value_units(&self, source: ValueUnits, target: ValueUnits) -> Self {
    self.to_runtime().to_value_units(source, target).into()
  }
}

impl From<TemperatureInfo> for TemperatureInfoFile {
  fn from(info: TemperatureInfo) -> Self {
    TemperatureInfoFile::from_runtime(&info)
  }
}

impl From<TemperatureInfoFile> for TemperatureInfo {
  fn from(file: TemperatureInfoFile) -> Self {
    file.to_runtime()
  }
}
