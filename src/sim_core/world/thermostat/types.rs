use crate::data::TimeIterationDistance;
use crate::data::units::TEMPERATURE_U;
use crate::data::units::TIME_U;
use crate::data::units::ValueUnits;
use crate::particle::Particle;

pub mod simple_temperature_info_generator;

pub const DEFAULT_TEMP_THRESHOLD_UNITLESS: f64 = 30. / TEMPERATURE_U;
pub const DEFAULT_ACCEPTANCE_TIME_UNITLESS: f64 = 2000. * 1e-18 / TIME_U;

pub trait TemperatureInfoGenerator {
  fn generate(&self) -> Result<Vec<TemperatureInfo>, String>;
}

/// A single [`TemperatureInfo`] entry or a generator that expands to many.
pub enum TemperatureInfoSource {
  Direct(TemperatureInfo),
  Generated(Box<dyn TemperatureInfoGenerator>),
}

impl From<TemperatureInfo> for TemperatureInfoSource {
  fn from(info: TemperatureInfo) -> Self {
    TemperatureInfoSource::Direct(info)
  }
}

impl<G: TemperatureInfoGenerator + 'static> From<G> for TemperatureInfoSource {
  fn from(generator: G) -> Self {
    TemperatureInfoSource::Generated(Box::new(generator))
  }
}

/// Expands each source in order and concatenates the results.
pub fn collect_temperature_infos(
  sources: &[TemperatureInfoSource],
) -> Result<Vec<TemperatureInfo>, String> {
  let mut result = Vec::new();

  for source in sources {
    match source {
      TemperatureInfoSource::Direct(info) => result.push(*info),
      TemperatureInfoSource::Generated(generator) => {
        result.extend(generator.generate()?);
      }
    }
  }

  Ok(result)
}

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

#[derive(Debug, Clone, Copy)]
pub struct TemperatureInfo {
  pub desired_temperature: f64,
  pub acceptance_distance: TimeIterationDistance,
  pub achieved_distance: TimeIterationDistance,
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

#[derive(Debug, Clone)]
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

#[cfg(test)]
mod collect_tests {
  use crate::data::TimeIterationDistance;

  use super::simple_temperature_info_generator::SimpleTemperatureInfoGenerator;
  use super::{collect_temperature_infos, TemperatureInfo, TemperatureInfoSource};

  #[test]
  fn preserves_order_across_direct_entries_and_generators() {
    let first = TemperatureInfo::new(5000., TimeIterationDistance::Iteration { value: 10 });
    let second = TemperatureInfo::new(4000., TimeIterationDistance::Iteration { value: 20 });
    let generator = SimpleTemperatureInfoGenerator::new()
      .start_temperature(3000.)
      .end_temperature(3000.)
      .acceptance_distance(TimeIterationDistance::Iteration { value: 30 })
      .achieved_distance(TimeIterationDistance::Iteration { value: 30 });
    let third = TemperatureInfo::new(1000., TimeIterationDistance::Iteration { value: 40 });

    let sources = [
      TemperatureInfoSource::from(first),
      TemperatureInfoSource::from(second),
      TemperatureInfoSource::from(generator),
      TemperatureInfoSource::from(third),
    ];

    let infos = collect_temperature_infos(&sources).unwrap();
    let temps: Vec<f64> = infos.iter().map(|i| i.desired_temperature).collect();

    assert_eq!(temps, vec![5000., 4000., 3000., 1000.]);
  }

  #[test]
  fn preserves_order_across_direct_entries_and_generators_v2() {
    let first = TemperatureInfo::new(5000., TimeIterationDistance::Iteration { value: 10 });
    let second = TemperatureInfo::new(4000., TimeIterationDistance::Iteration { value: 20 });
    let generator = SimpleTemperatureInfoGenerator::new()
      .start_temperature(1000.)
      .end_temperature(3000.)
      .acceptance_distance(TimeIterationDistance::Iteration { value: 0 })
      .achieved_distance(TimeIterationDistance::Iteration { value: 500 })
      .temperature_step(500.0);
    let third = TemperatureInfo::new(1000., TimeIterationDistance::Iteration { value: 40 });

    let sources = [
      TemperatureInfoSource::from(first),
      TemperatureInfoSource::from(second),
      TemperatureInfoSource::from(generator),
      TemperatureInfoSource::from(third),
    ];

    let infos = collect_temperature_infos(&sources).unwrap();
    let temps: Vec<f64> = infos.iter().map(|i| i.desired_temperature).collect();

    assert_eq!(temps, vec![5000., 4000., 1000., 1500., 2000., 2500., 3000., 1000.]);
  }

  #[test]
  fn expands_generator_sequence_in_place() {
    let generator = SimpleTemperatureInfoGenerator::new()
      .start_temperature(300.)
      .end_temperature(500.)
      .temperature_step(100.)
      .acceptance_distance(TimeIterationDistance::Iteration { value: 0 })
      .achieved_distance(TimeIterationDistance::Iteration { value: 0 });

    let sources = [TemperatureInfoSource::from(generator)];
    let infos = collect_temperature_infos(&sources).unwrap();
    let temps: Vec<f64> = infos.iter().map(|i| i.desired_temperature).collect();

    assert_eq!(temps, vec![300., 400., 500.]);
  }

  #[test]
  fn propagates_generator_errors() {
    let generator = SimpleTemperatureInfoGenerator::new();
    let sources = [TemperatureInfoSource::from(generator)];

    let err = collect_temperature_infos(&sources).unwrap_err();
    assert_eq!(err, "start_temperature must be set");
  }
}
