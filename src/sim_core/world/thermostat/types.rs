use crate::data::units::TEMPERATURE_U;
use crate::data::units::TIME_U;
use crate::data::units::ValueUnits;
use crate::particle::Particle;

pub const DEFAULT_TEMP_THRESHOLD_UNITLESS: f64 = 30. / TEMPERATURE_U;
pub const DEFAULT_ACCEPTANCE_TIME_UNITLESS: f64 = 2000. * 1e-18 / TIME_U;

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
#[serde(tag = "type")]
pub enum TimeIterationDistance {
  Time { value: f64 },
  Iteration { value: usize },
}

impl TimeIterationDistance {
  pub fn to_value_units(&self, source: ValueUnits, target: ValueUnits) -> Self {
    match self {
      TimeIterationDistance::Time { value } => TimeIterationDistance::Time {
        value: value * ValueUnits::scale_between(source, target, TIME_U),
      },
      TimeIterationDistance::Iteration { value } => {
        TimeIterationDistance::Iteration { value: *value }
      }
    }
  }
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
  pub save: bool,
}

impl TemperatureInfo {
  pub fn new(desired_temperature: f64, achieved_distance: TimeIterationDistance) -> Self {
    TemperatureInfo {
      desired_temperature,
      acceptance_distance: default_acceptance_distance(),
      achieved_distance,
      threshold: default_temperature_threshold(),
      save: false,
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
      save: false,
    }
  }

  pub fn to_value_units(&self, source: ValueUnits, target: ValueUnits) -> Self {
    TemperatureInfo {
      desired_temperature: self.desired_temperature
        * ValueUnits::scale_between(source, target, TEMPERATURE_U),
      acceptance_distance: self.acceptance_distance.to_value_units(source, target),
      achieved_distance: self.achieved_distance.to_value_units(source, target),
      threshold: self.threshold * ValueUnits::scale_between(source, target, TEMPERATURE_U),
      save: self.save,
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

impl IntegrationAlgorithm {
  pub fn to_value_units(&self, source: ValueUnits, target: ValueUnits) -> Self {
    match self {
      IntegrationAlgorithm::SemiImplicitEuler => IntegrationAlgorithm::SemiImplicitEuler,
      IntegrationAlgorithm::VelocityVerlet => IntegrationAlgorithm::VelocityVerlet,
      IntegrationAlgorithm::NoseHooverVerlet {
        desired_temperature,
        q_effective_mass,
      } => IntegrationAlgorithm::NoseHooverVerlet {
        desired_temperature: desired_temperature
          .iter()
          .map(|t| t.to_value_units(source, target))
          .collect(),
        q_effective_mass: *q_effective_mass,
      },
    }
  }
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

#[derive(Debug, Clone, Default)]
pub struct TemperatureHistoryEntry {
  pub temperature_started: Option<TemperatureIteration>,
  pub temperature_achieved: Option<TemperatureIteration>,
  pub temperature_switched: Option<TemperatureIteration>,
  pub particles: Option<Vec<Particle>>,
}

#[derive(Debug, Clone, Copy)]
pub enum IntegrationStateUpdateResponse {
  SemiImplicitEuler,
  VelocityVerlet,
  NoseHooverVerlet { updated: bool, temperature: f64 },
}
