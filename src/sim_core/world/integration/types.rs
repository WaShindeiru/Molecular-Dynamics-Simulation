use crate::data::units::TEMPERATURE_U;
use crate::data::units::TIME_U;

pub const DEFAULT_TEMP_THRESHOLD_UNITLESS: f64 = 30. / TEMPERATURE_U;
pub const DEFAULT_ACCEPTANCE_TIME_UNITLESS: f64 = 2000. * 1e-18 / TIME_U;

fn default_acceptance_distance() -> TimeIterationDistance {
  TimeIterationDistance::Time { value: DEFAULT_ACCEPTANCE_TIME_UNITLESS }
}

fn default_achieved_distance() -> TimeIterationDistance {
  TimeIterationDistance::Iteration { value: 0 }
}

fn default_temperature_threshold() -> f64 {
  DEFAULT_TEMP_THRESHOLD_UNITLESS
}

#[derive(Debug, Clone, Copy, serde::Serialize, serde::Deserialize)]
#[serde(tag = "type")]
pub enum TimeIterationDistance {
  Time { value: f64 },
  Iteration { value: usize },
}

#[derive(Debug, Clone, Copy, serde::Serialize, serde::Deserialize)]
pub struct TemperatureInfo {
  pub desired_temperature: f64,
  #[serde(default = "default_acceptance_distance")]
  pub acceptance_distance: TimeIterationDistance,
  #[serde(alias = "distance", default = "default_achieved_distance")]
  pub achieved_distance: TimeIterationDistance,
  #[serde(default = "default_temperature_threshold")]
  pub threshold: f64,
}

impl TemperatureInfo {
  pub fn new(desired_temperature: f64, achieved_distance: TimeIterationDistance) -> Self {
    TemperatureInfo {
      desired_temperature,
      acceptance_distance: default_acceptance_distance(),
      achieved_distance,
      threshold: default_temperature_threshold(),
    }
  }

  pub fn with_params(
    desired_temperature: f64,
    acceptance_distance: TimeIterationDistance,
    achieved_distance: TimeIterationDistance,
    threshold: f64,
  ) -> Self {
    TemperatureInfo {
      desired_temperature,
      acceptance_distance,
      achieved_distance,
      threshold,
    }
  }
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
#[serde(tag = "type")]
pub enum IntegrationAlgorithm {
  SemiImplicitEuler,
  VelocityVerlet,
  NoseHooverVerlet {
    desired_temperature: Vec<TemperatureInfo>,
    q_effective_mass: f64,
  },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NoseHooverStage {
  StabilizeTemperature,
  TemperatureAchieved,
  LastEntry,
}

#[derive(Debug, Clone, Copy)]
pub struct TemperatureIteration {
  pub iteration: usize,
  pub temperature: f64,
}

#[derive(Debug, Clone, Copy, Default)]
pub struct TemperatureHistoryEntry {
  pub temperature_started: Option<TemperatureIteration>,
  pub temperature_achieved: Option<TemperatureIteration>,
  pub temperature_switched: Option<TemperatureIteration>,
}

#[derive(Debug, Clone, Copy)]
pub enum IntegrationStateUpdateResponse {
  SemiImplicitEuler,
  VelocityVerlet,
  NoseHooverVerlet {
    updated: bool,
    temperature: f64,
  }
}
