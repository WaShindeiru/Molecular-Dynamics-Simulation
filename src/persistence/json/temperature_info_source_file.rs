use crate::sim_core::world::thermostat::{
  collect_temperature_infos, TemperatureInfo, TemperatureInfoSource,
};

use super::simple_temperature_info_generator_file::SimpleTemperatureInfoGeneratorFile;
use super::temperature_info_file::TemperatureInfoFile;

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
#[serde(tag = "type")]
pub enum TemperatureInfoSourceFile {
  #[serde(rename = "TemperatureInfo")]
  Direct(TemperatureInfoFile),
  #[serde(rename = "SimpleTemperatureInfoGenerator")]
  Simple(SimpleTemperatureInfoGeneratorFile),
}

impl TemperatureInfoSourceFile {
  pub fn from_runtime(info: &TemperatureInfo) -> Self {
    TemperatureInfoSourceFile::Direct(TemperatureInfoFile::from_runtime(info))
  }

  pub fn to_runtime_source(&self) -> Result<TemperatureInfoSource, String> {
    match self {
      TemperatureInfoSourceFile::Direct(file) => {
        Ok(TemperatureInfoSource::Direct(file.to_runtime()))
      }
      TemperatureInfoSourceFile::Simple(file) => {
        Ok(TemperatureInfoSource::Generated(Box::new(file.to_generator())))
      }
    }
  }
}

pub fn collect_temperature_infos_from_file(
  sources: &[TemperatureInfoSourceFile],
) -> Result<Vec<TemperatureInfo>, String> {
  let runtime_sources: Vec<TemperatureInfoSource> = sources
    .iter()
    .map(|source| source.to_runtime_source())
    .collect::<Result<_, _>>()?;

  collect_temperature_infos(&runtime_sources)
}
