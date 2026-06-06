use crate::data::config::correction_param::{CorrectionParam, SmallDistanceCorrection};
use crate::data::units::{R_U, ValueUnits};

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize, Default)]
pub struct CorrectionParamFile {
  #[serde(default)]
  pub small_distance: Option<SmallDistanceCorrectionFile>,
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct SmallDistanceCorrectionFile {
  pub enabled: bool,
  pub substep_count: usize,
  pub distance_threshold: f64,
}

impl CorrectionParamFile {
  pub fn from_runtime(value: CorrectionParam) -> Self {
    Self {
      small_distance: Some(SmallDistanceCorrectionFile::from_runtime(value.small_distance)),
    }
  }

  pub fn to_runtime(&self) -> CorrectionParam {
    CorrectionParam {
      small_distance: self
        .small_distance
        .as_ref()
        .map(|sd| sd.to_runtime())
        .unwrap_or_default(),
    }
  }

  pub fn to_value_units(&self, source: ValueUnits, target: ValueUnits) -> Self {
    Self {
      small_distance: self.small_distance.as_ref().map(|sd| sd.to_value_units(source, target)),
    }
  }
}

impl SmallDistanceCorrectionFile {
  pub fn from_runtime(value: SmallDistanceCorrection) -> Self {
    Self {
      enabled: value.enabled,
      substep_count: value.substep_count,
      distance_threshold: value.distance_threshold,
    }
  }

  pub fn to_runtime(&self) -> SmallDistanceCorrection {
    SmallDistanceCorrection {
      enabled: self.enabled,
      substep_count: self.substep_count,
      distance_threshold: self.distance_threshold,
    }
  }

  pub fn to_value_units(&self, source: ValueUnits, target: ValueUnits) -> Self {
    Self {
      enabled: self.enabled,
      substep_count: self.substep_count,
      distance_threshold: self.distance_threshold
        * ValueUnits::scale_between(source, target, R_U),
    }
  }
}
