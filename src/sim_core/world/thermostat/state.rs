use std::io;
use std::path::Path;

use crate::data::units::TEMPERATURE_U;
use crate::data::{ParticleConfig, TimeIterationDistance, ValueUnits};
use crate::particle::Particle;
use crate::persistence::json::control_velocity_manager_file::ControlVelocityManagerFile;
use crate::persistence::json::particle_config::ParticleConfigFile;
use crate::persistence::json::velocity_manager_file::VelocityManagerFile;

use super::types::{
  IntegrationAlgorithm, IntegrationStateUpdateResponse, NoseHooverStage, TemperatureHistoryEntry,
  TemperatureIteration,
};

pub enum IntegrationAlgorithmState {
  SemiImplicitEuler,
  VelocityVerlet,
  NoseHooverVerlet {
    temperature_index: usize,
    stage: NoseHooverStage,
    acceptance_consecutive: usize,
    after_achieved_consecutive: usize,
    history: Vec<TemperatureHistoryEntry>,
  },
}

pub fn new_integration_algorithm_state(
  integration_algorithm: &IntegrationAlgorithm,
) -> IntegrationAlgorithmState {
  match integration_algorithm {
    IntegrationAlgorithm::SemiImplicitEuler => IntegrationAlgorithmState::SemiImplicitEuler,
    IntegrationAlgorithm::VelocityVerlet => IntegrationAlgorithmState::VelocityVerlet,
    IntegrationAlgorithm::NoseHooverVerlet {
      desired_temperature,
      ..
    } => {
      let history = vec![TemperatureHistoryEntry::default(); desired_temperature.len()];
      let stage = if desired_temperature.len() <= 1 {
        NoseHooverStage::LastEntry
      } else {
        NoseHooverStage::StabilizeTemperature
      };

      IntegrationAlgorithmState::NoseHooverVerlet {
        temperature_index: 0,
        stage,
        acceptance_consecutive: 0,
        after_achieved_consecutive: 0,
        history,
      }
    }
  }
}

impl IntegrationAlgorithmState {
  pub fn temperature_history(&self) -> Option<&[TemperatureHistoryEntry]> {
    match self {
      IntegrationAlgorithmState::NoseHooverVerlet { history, .. } => Some(history.as_slice()),
      _ => None,
    }
  }

  pub fn get_previous_history_entry(&self) -> Option<&TemperatureHistoryEntry> {
    match self {
      IntegrationAlgorithmState::SemiImplicitEuler => None,
      IntegrationAlgorithmState::VelocityVerlet => None,
      IntegrationAlgorithmState::NoseHooverVerlet {
        temperature_index,
        stage,
        acceptance_consecutive,
        after_achieved_consecutive,
        history,
      } => Some(history.get(*temperature_index - 1).unwrap()),
    }
  }

  fn distance_reached(distance: TimeIterationDistance, time_step: f64, consecutive: usize) -> bool {
    match distance {
      TimeIterationDistance::Time { value: time } => time_step * consecutive as f64 > time,
      TimeIterationDistance::Iteration { value: iteration } => consecutive > iteration,
    }
  }

  /// An `acceptance_distance` of zero means "don't wait on the actual temperature at all" -
  /// the stage should move straight to `TemperatureAchieved` and let `achieved_distance`
  /// govern the timing instead of the simulated temperature.
  fn skips_acceptance_check(distance: TimeIterationDistance) -> bool {
    match distance {
      TimeIterationDistance::Time { value } => value == 0.0,
      TimeIterationDistance::Iteration { value } => value == 0,
    }
  }

  fn response_nose_hoover(updated: bool, temperature: f64) -> IntegrationStateUpdateResponse {
    IntegrationStateUpdateResponse::NoseHooverVerlet {
      updated,
      temperature,
    }
  }

  fn record_started_if_missing(
    history: &mut [TemperatureHistoryEntry],
    index: usize,
    iteration: usize,
    simulation_temperature_unitless: f64,
  ) {
    if history[index].temperature_started.is_none() {
      history[index].temperature_started = Some(TemperatureIteration {
        iteration,
        temperature: simulation_temperature_unitless,
      });
    }
  }

  fn record_achieved(
    history: &mut [TemperatureHistoryEntry],
    index: usize,
    iteration: usize,
    simulation_temperature_unitless: f64,
  ) {
    history[index].temperature_achieved = Some(TemperatureIteration {
      iteration,
      temperature: simulation_temperature_unitless,
    });
  }

  fn record_switched<'a>(
    history: &mut [TemperatureHistoryEntry],
    index: usize,
    iteration: usize,
    simulation_temperature_unitless: f64,
    save: bool,
    particles: impl Iterator<Item = &'a Particle>,
  ) {
    history[index].temperature_switched = Some(TemperatureIteration {
      iteration,
      temperature: simulation_temperature_unitless,
    });
    if save {
      history[index].particles = Some(particles.cloned().collect());
    }
  }

  pub fn update_state<'a>(
    &mut self,
    current_iteration: usize,
    time_step: f64,
    integration_algorithm: &IntegrationAlgorithm,
    simulation_temperature_unitless: f64,
    particles: impl Iterator<Item = &'a Particle>,
  ) -> IntegrationStateUpdateResponse {

    match (self, integration_algorithm) {
      (IntegrationAlgorithmState::SemiImplicitEuler, IntegrationAlgorithm::SemiImplicitEuler) => {
        IntegrationStateUpdateResponse::SemiImplicitEuler
      }

      (IntegrationAlgorithmState::VelocityVerlet, IntegrationAlgorithm::VelocityVerlet) => {
        IntegrationStateUpdateResponse::VelocityVerlet
      }

      (
        IntegrationAlgorithmState::NoseHooverVerlet {
          temperature_index,
          stage,
          acceptance_consecutive,
          after_achieved_consecutive,
          history,
        },
        IntegrationAlgorithm::NoseHooverVerlet {
          desired_temperature,
          q_effective_mass: _q_effective_mass,
        },
      ) => {
        let temp_info = *desired_temperature.get(*temperature_index).unwrap();
        IntegrationAlgorithmState::record_started_if_missing(
          history,
          *temperature_index,
          current_iteration,
          simulation_temperature_unitless,
        );

        match *stage {
          NoseHooverStage::LastEntry => {
            IntegrationAlgorithmState::response_nose_hoover(false, temp_info.desired_temperature)
          }
          NoseHooverStage::StabilizeTemperature => {
            let within_threshold =
              IntegrationAlgorithmState::skips_acceptance_check(temp_info.acceptance_distance)
                || (simulation_temperature_unitless - temp_info.desired_temperature).abs()
                  < temp_info.threshold;

            if within_threshold {
              *acceptance_consecutive += 1;

              if IntegrationAlgorithmState::distance_reached(
                temp_info.acceptance_distance,
                time_step,
                *acceptance_consecutive,
              ) {
                *stage = NoseHooverStage::TemperatureAchieved;
                *after_achieved_consecutive = 0;
                IntegrationAlgorithmState::record_achieved(
                  history,
                  *temperature_index,
                  current_iteration,
                  simulation_temperature_unitless,
                );
              }
            } else {
              *acceptance_consecutive = 0;
            }

            IntegrationAlgorithmState::response_nose_hoover(false, temp_info.desired_temperature)
          }
          NoseHooverStage::TemperatureAchieved => {
            *after_achieved_consecutive += 1;

            let should_switch = IntegrationAlgorithmState::distance_reached(
              temp_info.achieved_distance,
              time_step,
              *after_achieved_consecutive,
            );

            if !should_switch {
              return IntegrationAlgorithmState::response_nose_hoover(
                false,
                temp_info.desired_temperature,
              );
            }

            IntegrationAlgorithmState::record_switched(
              history,
              *temperature_index,
              current_iteration,
              simulation_temperature_unitless,
              temp_info.save,
              particles,
            );
            *acceptance_consecutive = 0;
            *after_achieved_consecutive = 0;

            let next_temperature_index = *temperature_index + 1;
            debug_assert!(
              next_temperature_index < desired_temperature.len(),
              "Invalid NoseHoover state: switch requested from last temperature entry"
            );

            *temperature_index = next_temperature_index;
            let next_temp_info = *desired_temperature.get(*temperature_index).unwrap();

            *stage = if *temperature_index == desired_temperature.len() - 1 {
              NoseHooverStage::LastEntry
            } else {
              NoseHooverStage::StabilizeTemperature
            };

            IntegrationAlgorithmState::response_nose_hoover(
              true,
              next_temp_info.desired_temperature,
            )
          }
        }
      }

      _ => {
        panic!("IntegrationAlgorithmState does not match IntegrationAlgorithm variant");
      }
    }
  }

  pub fn save_temperature_particles(
    &self,
    integration_algorithm: &IntegrationAlgorithm,
    base_path: &Path,
    velocity_managers_file: Vec<VelocityManagerFile>,
    control_velocity_managers_file: Vec<ControlVelocityManagerFile>,
  ) -> io::Result<()> {
    let (history, desired_temperature) = match (self, integration_algorithm) {
      (
        IntegrationAlgorithmState::NoseHooverVerlet { history, .. },
        IntegrationAlgorithm::NoseHooverVerlet { desired_temperature, .. },
      ) => (history, desired_temperature),
      _ => return Ok(()),
    };

    let dir_path = base_path.join("particles").join("temperature");

    for (entry, temp_info) in history.iter().zip(desired_temperature.iter()) {
      let Some(particles) = &entry.particles else {
        continue;
      };

      let temperature_k = temp_info.desired_temperature * TEMPERATURE_U;
      let file_path = dir_path.join(format!("{}.json", temperature_k));

      std::fs::create_dir_all(&dir_path)?;

      let velocity_schedules = velocity_managers_file.iter().map(VelocityManagerFile::to_schedule).collect();
      let control_velocity_schedules = control_velocity_managers_file
        .iter()
        .map(ControlVelocityManagerFile::to_schedule)
        .collect();
      let config =
        ParticleConfig::new_with_all_schedules(particles.clone(), velocity_schedules, control_velocity_schedules);
      let config_file = ParticleConfigFile::from_runtime(&config, ValueUnits::Si);
      let json = serde_json::to_string_pretty(&config_file)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
      std::fs::write(file_path, json)?;
    }

    Ok(())
  }
}
