use std::fmt;
use crate::data::units::TIME_U;
use crate::data::units::TEMPERATURE_U;

const TEMP_THRESHOLD: f64 = 15.;
const TEMP_THRESHOLD_UNITLESS: f64 = TEMP_THRESHOLD / TEMPERATURE_U;
const ACCEPTANCE_TIME_UNITLESS: f64 = 80. * 1e-18 / TIME_U;

#[derive(Debug, Clone, Copy)]
pub enum TimeIterationDistance {
  Time (f64),
  Iteration (usize),
}

#[derive(Debug, Clone, Copy)]
pub struct TemperatureInfo {
  pub desired_temperature: f64,
  pub distance: TimeIterationDistance,
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

pub enum IntegrationAlgorithmState {
  SemiImplicitEuler,
  VelocityVerlet,
  NoseHooverVerlet {
    temperature_index: usize,
    consecutive: usize,
    achieved: bool
  }
}

pub fn new_integration_algorithm_state(integration_algorithm: &IntegrationAlgorithm) -> IntegrationAlgorithmState {
  match integration_algorithm {
    IntegrationAlgorithm::SemiImplicitEuler => IntegrationAlgorithmState::SemiImplicitEuler,
    IntegrationAlgorithm::VelocityVerlet => IntegrationAlgorithmState::VelocityVerlet,
    IntegrationAlgorithm::NoseHooverVerlet { .. } => IntegrationAlgorithmState::NoseHooverVerlet {
      temperature_index: 0,
      achieved: false,
      consecutive: 0,
    }
  }
}

#[derive(Debug, Clone, Copy)]
pub enum IntegrationStateUpdateResponse {
  SemiImplicitEuler,
  VelocityVerlet,
  NoseHooverVerlet {
    updated: bool,
    temperature: f64,
    end: bool,
  }
}

impl IntegrationAlgorithmState {
  pub fn update_state(&mut self, time_step: f64, integration_algorithm: &IntegrationAlgorithm,
                      simulation_temperature_unitless: f64) -> IntegrationStateUpdateResponse {
    match (self, integration_algorithm) {
      (IntegrationAlgorithmState::SemiImplicitEuler, IntegrationAlgorithm::SemiImplicitEuler) => {
        IntegrationStateUpdateResponse::SemiImplicitEuler
      }

      (IntegrationAlgorithmState::VelocityVerlet, IntegrationAlgorithm::VelocityVerlet) => {
        IntegrationStateUpdateResponse::VelocityVerlet
      }

      (IntegrationAlgorithmState::NoseHooverVerlet { temperature_index, consecutive, achieved },
       IntegrationAlgorithm::NoseHooverVerlet { desired_temperature, q_effective_mass: _q_effective_mass }) => {

        let temp_info = *desired_temperature.get(*temperature_index).unwrap();

        if *temperature_index == desired_temperature.len() - 1 {
          return IntegrationStateUpdateResponse::NoseHooverVerlet{
            temperature: temp_info.desired_temperature,
            updated: false,
            end: true,
          }
        }

        if !*achieved {
          if (simulation_temperature_unitless - temp_info.desired_temperature).abs() < TEMP_THRESHOLD_UNITLESS {
            *consecutive = *consecutive + 1;
            let current_time = time_step * *consecutive as f64;

            if current_time > ACCEPTANCE_TIME_UNITLESS {
              *achieved = true;
            }

            IntegrationStateUpdateResponse::NoseHooverVerlet{
              temperature: temp_info.desired_temperature,
              updated: false,
              end: false,
            }

          } else {
            *consecutive = 0;
            IntegrationStateUpdateResponse::NoseHooverVerlet{
              temperature: temp_info.desired_temperature,
              updated: false,
              end: false,
            }

          }

        } else {
          *consecutive = *consecutive + 1;

          let switch = match temp_info.distance {
            TimeIterationDistance::Time(time) => {
              let current_time = time_step * *consecutive as f64;
              current_time > time
            }
            TimeIterationDistance::Iteration(iteration) => *consecutive > iteration
          };

          if switch {
            *achieved = false;

            match desired_temperature.get(*temperature_index + 1) {
              Some(temp) => {
                *temperature_index = *temperature_index + 1;
                IntegrationStateUpdateResponse::NoseHooverVerlet{
                  temperature: temp.desired_temperature,
                  updated: true,
                  end: false,
                }
              }
              _ => {
                IntegrationStateUpdateResponse::NoseHooverVerlet{
                  temperature: temp_info.desired_temperature,
                  updated: true,
                  end: true,
                }
              }
            }

          } else {
            IntegrationStateUpdateResponse::NoseHooverVerlet{
              temperature: temp_info.desired_temperature,
              updated: false,
              end: false,
            }
          }
        }
      }

      _ => {
        panic!("IntegrationAlgorithmState does not match IntegrationAlgorithm variant");
      }
    }
  }
}

impl fmt::Display for IntegrationAlgorithm {

  fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
    match self {
      IntegrationAlgorithm::SemiImplicitEuler => write!(f, "SemiImplicitEuler"),
      IntegrationAlgorithm::VelocityVerlet => write!(f, "VelocityVerlet"),
      IntegrationAlgorithm::NoseHooverVerlet {
        desired_temperature: _t, q_effective_mass: _q
      } => {
        // let s = t.iter()
        //   .map(|v| v.to_string())
        //   .collect::<Vec<_>>()
        //   .join(", ");
        // write!(f, "NoseHooverVerlet temperature: {s}, q_effective_mass: {q}", )
        write!(f, "NoseHooverVerlet")
      },
    }
  }
}


