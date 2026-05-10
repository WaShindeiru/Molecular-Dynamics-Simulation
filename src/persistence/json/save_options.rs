use crate::data::units::ValueUnits;
use crate::sim_core::world::saver::SaveOptions;

use super::frame_sampling::FrameSamplingConfigFile;
use super::simulation_config::refresh_save_path_with_current_date;

#[derive(Clone, serde::Serialize, serde::Deserialize)]
pub struct SaveOptionsFile {
  pub save: bool,
  pub save_path: String,
  #[serde(default)]
  pub keep_path: bool,
  pub save_laamps: bool,
  pub save_verbose: bool,
  pub save_all_iterations_laamps: bool,
  pub save_all_iterations_energy: bool,
  pub laamps_sampling: FrameSamplingConfigFile,
  pub energy_sampling: FrameSamplingConfigFile,
}

impl SaveOptionsFile {
  pub fn from_runtime(save_options: &SaveOptions) -> Self {
    Self {
      save: save_options.save,
      save_path: save_options.save_path.clone(),
      keep_path: false,
      save_laamps: save_options.save_laamps,
      save_verbose: save_options.save_verbose,
      save_all_iterations_laamps: save_options.save_all_iterations_laamps,
      save_all_iterations_energy: save_options.save_all_iterations_energy,
      laamps_sampling: FrameSamplingConfigFile::from_runtime(&save_options.laamps_sampling),
      energy_sampling: FrameSamplingConfigFile::from_runtime(&save_options.energy_sampling),
    }
  }

  pub fn to_runtime(&self) -> SaveOptions {
    let save_path = if self.keep_path {
      self.save_path.clone()
    } else {
      refresh_save_path_with_current_date(&self.save_path)
    };

    SaveOptions {
      save: self.save,
      save_path,
      save_laamps: self.save_laamps,
      save_verbose: self.save_verbose,
      save_all_iterations_laamps: self.save_all_iterations_laamps,
      save_all_iterations_energy: self.save_all_iterations_energy,
      laamps_sampling: self.laamps_sampling.to_runtime(),
      energy_sampling: self.energy_sampling.to_runtime(),
    }
  }

  pub fn to_value_units(&self, source: ValueUnits, target: ValueUnits) -> Self {
    Self {
      save: self.save,
      save_path: self.save_path.clone(),
      keep_path: self.keep_path,
      save_laamps: self.save_laamps,
      save_verbose: self.save_verbose,
      save_all_iterations_laamps: self.save_all_iterations_laamps,
      save_all_iterations_energy: self.save_all_iterations_energy,
      laamps_sampling: self.laamps_sampling.to_value_units(source, target),
      energy_sampling: self.energy_sampling.to_value_units(source, target),
    }
  }
}
