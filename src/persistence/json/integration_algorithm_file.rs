use crate::data::units::ValueUnits;
use crate::sim_core::world::thermostat::IntegrationAlgorithm;

use super::temperature_info_source_file::{
  collect_temperature_infos_from_file, TemperatureInfoSourceFile,
};

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
#[serde(tag = "type")]
pub enum IntegrationAlgorithmFile {
  SemiImplicitEuler,
  VelocityVerlet,
  NoseHooverVerlet {
    desired_temperature: Vec<TemperatureInfoSourceFile>,
    q_effective_mass: f64,
  },
}

impl IntegrationAlgorithmFile {
  pub fn from_runtime(algorithm: &IntegrationAlgorithm) -> Self {
    match algorithm {
      IntegrationAlgorithm::SemiImplicitEuler => IntegrationAlgorithmFile::SemiImplicitEuler,
      IntegrationAlgorithm::VelocityVerlet => IntegrationAlgorithmFile::VelocityVerlet,
      IntegrationAlgorithm::NoseHooverVerlet {
        desired_temperature,
        q_effective_mass,
      } => IntegrationAlgorithmFile::NoseHooverVerlet {
        desired_temperature: desired_temperature
          .iter()
          .map(TemperatureInfoSourceFile::from_runtime)
          .collect(),
        q_effective_mass: *q_effective_mass,
      },
    }
  }

  pub fn to_runtime(&self) -> Result<IntegrationAlgorithm, String> {
    match self {
      IntegrationAlgorithmFile::SemiImplicitEuler => Ok(IntegrationAlgorithm::SemiImplicitEuler),
      IntegrationAlgorithmFile::VelocityVerlet => Ok(IntegrationAlgorithm::VelocityVerlet),
      IntegrationAlgorithmFile::NoseHooverVerlet {
        desired_temperature,
        q_effective_mass,
      } => Ok(IntegrationAlgorithm::NoseHooverVerlet {
        desired_temperature: collect_temperature_infos_from_file(desired_temperature)?,
        q_effective_mass: *q_effective_mass,
      }),
    }
  }

  pub fn to_value_units(&self, source: ValueUnits, target: ValueUnits) -> Self {
    match self {
      IntegrationAlgorithmFile::SemiImplicitEuler => IntegrationAlgorithmFile::SemiImplicitEuler,
      IntegrationAlgorithmFile::VelocityVerlet => IntegrationAlgorithmFile::VelocityVerlet,
      IntegrationAlgorithmFile::NoseHooverVerlet {
        desired_temperature,
        q_effective_mass,
      } => IntegrationAlgorithmFile::NoseHooverVerlet {
        desired_temperature: desired_temperature
          .iter()
          .map(|entry| match entry {
            TemperatureInfoSourceFile::Direct(file) => {
              TemperatureInfoSourceFile::Direct(file.to_value_units(source, target))
            }
            TemperatureInfoSourceFile::Simple(file) => {
              TemperatureInfoSourceFile::Simple(file.to_value_units(source, target))
            }
          })
          .collect(),
        q_effective_mass: *q_effective_mass,
      },
    }
  }
}

#[cfg(test)]
mod tests {
  use crate::data::TimeIterationDistance;

  use super::*;

  #[test]
  fn deserializes_mixed_temperature_sources_preserving_order() {
    let json = r#"{
      "type": "NoseHooverVerlet",
      "q_effective_mass": 1.0,
      "desired_temperature": [
        {
          "type": "TemperatureInfo",
          "desired_temperature": 5000.0,
          "achieved_distance": { "type": "Iteration", "value": 10 },
          "save": false
        },
        {
          "type": "SimpleTemperatureInfoGenerator",
          "start_temperature": 1000.0,
          "end_temperature": 3000.0,
          "temperature_step": 500.0,
          "acceptance_distance": { "type": "Iteration", "value": 0 },
          "achieved_distance": { "type": "Iteration", "value": 500 },
          "save_step": 2
        },
        {
          "type": "TemperatureInfo",
          "desired_temperature": 1000.0,
          "achieved_distance": { "type": "Iteration", "value": 40 },
          "save": true
        }
      ]
    }"#;

    let file: IntegrationAlgorithmFile = serde_json::from_str(json).unwrap();
    let runtime = file.to_runtime().unwrap();

    let IntegrationAlgorithm::NoseHooverVerlet { desired_temperature, .. } = runtime else {
      panic!("expected NoseHooverVerlet");
    };

    let temps: Vec<f64> = desired_temperature
      .iter()
      .map(|info| info.desired_temperature)
      .collect();

    assert_eq!(
      temps,
      vec![5000.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0, 1000.0]
    );
    assert!(!desired_temperature[0].save);
    assert!(desired_temperature.last().unwrap().save);
    assert_eq!(
      desired_temperature[1].achieved_distance,
      TimeIterationDistance::Iteration { value: 500 }
    );
  }
}
