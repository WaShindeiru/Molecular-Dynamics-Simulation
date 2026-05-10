use std::path::Path;
use std::{fs, io};

use nalgebra::Vector3;
use serde_json::Value;

use crate::data::SimulationConfig;
use crate::data::units::{R_U, TIME_U, ValueUnits};
use crate::sim_core::world::integration::IntegrationAlgorithm;
use crate::sim_core::world::saver::FrameSamplingConfig;
use crate::utils::logging::get_save_path;

use super::save_options::SaveOptionsFile;
use super::types::{EdgeConditionFile, WorldTypeFile};

pub(crate) fn refresh_save_path_with_current_date(existing_path: &str) -> String {
  let path = Path::new(existing_path);
  let prefix_path = path
    .parent()
    .filter(|p| !p.as_os_str().is_empty())
    .unwrap_or(path);

  let mut prefix = prefix_path.to_string_lossy().into_owned();
  if !prefix.ends_with('/') {
    prefix.push('/');
  }

  get_save_path(prefix)
}

fn validate_sampling_if_present(
  config: &SimulationConfig,
  raw_json: &Value,
  frame_count_pointer: &str,
  sampling: &FrameSamplingConfig,
  save_all_iterations: bool,
  sampling_name: &str,
) -> io::Result<()> {
  if raw_json.pointer(frame_count_pointer).is_none() {
    return Ok(());
  }

  if config.time_step <= 0.0 {
    return Err(io::Error::new(
      io::ErrorKind::InvalidData,
      format!(
        "Cannot validate {}: time_step must be > 0, got {}",
        sampling_name, config.time_step
      ),
    ));
  }

  let expected =
    FrameSamplingConfig::iterations_from_duration(sampling.one_frame_duration, config.time_step);
  let actual = sampling.frame_iteration_count;

  let valid = if save_all_iterations {
    actual == 1
  } else {
    actual == expected
  };

  if !valid {
    let expected_message = if save_all_iterations {
      "1 (save-all mode)"
    } else {
      "nearest rounding of one_frame_duration/time_step"
    };
    return Err(io::Error::new(
      io::ErrorKind::InvalidData,
      format!(
        "Invalid {}.frame_iteration_count: got {}, expected {} [{}]. one_frame_duration={}, time_step={}, save_all_iterations={}",
        sampling_name,
        actual,
        expected,
        expected_message,
        sampling.one_frame_duration,
        config.time_step,
        save_all_iterations
      ),
    ));
  }

  Ok(())
}

#[derive(Clone, serde::Serialize, serde::Deserialize)]
pub struct SimulationConfigFile {
  #[serde(default)]
  pub value_units: ValueUnits,
  pub world_size: Vector3<f64>,
  pub potential_gravity_max: f64,
  pub time_step: f64,
  pub num_of_iterations: usize,
  pub max_iteration_till_reset: usize,
  pub save_options: SaveOptionsFile,
  pub integration_algorithm: IntegrationAlgorithm,
  pub world_type: WorldTypeFile,
  pub edge_condition: EdgeConditionFile,
}

impl SimulationConfigFile {
  pub fn from_runtime(config: &SimulationConfig, target_units: ValueUnits) -> Self {
    let unitless = Self {
      value_units: ValueUnits::Unitless,
      world_size: config.world_size,
      potential_gravity_max: config.potential_gravity_max,
      time_step: config.time_step,
      num_of_iterations: config.num_of_iterations,
      max_iteration_till_reset: config.max_iteration_till_reset,
      save_options: SaveOptionsFile::from_runtime(&config.save_options),
      integration_algorithm: config.integration_algorithm.clone(),
      world_type: WorldTypeFile::from_runtime(config.world_type),
      edge_condition: EdgeConditionFile::from_runtime(config.edge_condition),
    };

    unitless.to_value_units(target_units)
  }

  pub fn into_runtime_unitless(self) -> SimulationConfig {
    let unitless = self.to_value_units(ValueUnits::Unitless);
    debug_assert_eq!(unitless.value_units, ValueUnits::Unitless);

    SimulationConfig {
      world_size: unitless.world_size,
      potential_gravity_max: unitless.potential_gravity_max,
      time_step: unitless.time_step,
      num_of_iterations: unitless.num_of_iterations,
      max_iteration_till_reset: unitless.max_iteration_till_reset,
      save_options: unitless.save_options.to_runtime(),
      integration_algorithm: unitless.integration_algorithm,
      world_type: unitless.world_type.to_runtime(),
      edge_condition: unitless.edge_condition.to_runtime(),
    }
  }

  /// Unit conversion, [`SaveOptionsFile`] -> [`SaveOptions`] (including `keep_path` / path refresh), and sampling validation when `frame_iteration_count` is present in `raw_json`.
  pub fn try_into_simulation_config(self, raw_json: &Value) -> io::Result<SimulationConfig> {
    let config = self.into_runtime_unitless();

    validate_sampling_if_present(
      &config,
      raw_json,
      "/save_options/laamps_sampling/frame_iteration_count",
      &config.save_options.laamps_sampling,
      config.save_options.save_all_iterations_laamps,
      "laamps_sampling",
    )?;
    validate_sampling_if_present(
      &config,
      raw_json,
      "/save_options/energy_sampling/frame_iteration_count",
      &config.save_options.energy_sampling,
      config.save_options.save_all_iterations_energy,
      "energy_sampling",
    )?;

    Ok(config)
  }

  pub fn to_value_units(&self, target: ValueUnits) -> Self {
    let source = self.value_units;

    Self {
      value_units: target,
      world_size: self.world_size * ValueUnits::scale_between(source, target, R_U),
      potential_gravity_max: self.potential_gravity_max,
      time_step: self.time_step * ValueUnits::scale_between(source, target, TIME_U),
      num_of_iterations: self.num_of_iterations,
      max_iteration_till_reset: self.max_iteration_till_reset,
      save_options: self.save_options.to_value_units(source, target),
      integration_algorithm: self.integration_algorithm.to_value_units(source, target),
      world_type: self.world_type,
      edge_condition: self.edge_condition,
    }
  }

  pub fn to_json_string(&self) -> Result<String, serde_json::Error> {
    serde_json::to_string_pretty(self)
  }

  pub fn from_json_str(s: &str) -> Result<Self, serde_json::Error> {
    serde_json::from_str(s)
  }

  pub fn to_json_file(&self, path: &str) -> io::Result<()> {
    let json = self
      .to_json_string()
      .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    fs::write(path, json)
  }

  pub fn from_json_file(path: &str) -> io::Result<Self> {
    let content = fs::read_to_string(path)?;
    Self::from_json_str(&content).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
  }
}

pub fn read_simulation_config_from_json_str(content: &str) -> io::Result<SimulationConfig> {
  let raw: Value =
    serde_json::from_str(content).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
  let config_file: SimulationConfigFile = serde_json::from_value(raw.clone())
    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
  config_file.try_into_simulation_config(&raw)
}

pub fn read_simulation_config_from_json_file(path: &str) -> io::Result<SimulationConfig> {
  let content = fs::read_to_string(path)?;
  read_simulation_config_from_json_str(&content)
}
