use crate::data::TimeIterationDistance;
use crate::data::units::{TEMPERATURE_U, ValueUnits};
use crate::sim_core::world::thermostat::SimpleTemperatureInfoGenerator;

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct SimpleTemperatureInfoGeneratorFile {
  pub start_temperature: f64,
  pub end_temperature: f64,
  pub acceptance_distance: TimeIterationDistance,
  pub achieved_distance: TimeIterationDistance,
  #[serde(default)]
  pub temperature_step: Option<f64>,
  #[serde(default)]
  pub threshold: Option<f64>,
  #[serde(default)]
  pub save_step: Option<usize>,
}

impl SimpleTemperatureInfoGeneratorFile {
  pub fn to_generator(&self) -> SimpleTemperatureInfoGenerator {
    let mut generator = SimpleTemperatureInfoGenerator::new()
      .start_temperature(self.start_temperature)
      .end_temperature(self.end_temperature)
      .acceptance_distance(self.acceptance_distance)
      .achieved_distance(self.achieved_distance);

    if let Some(temperature_step) = self.temperature_step {
      generator = generator.temperature_step(temperature_step);
    }
    if let Some(threshold) = self.threshold {
      generator = generator.threshold(threshold);
    }
    if let Some(save_step) = self.save_step {
      generator = generator.save_step(save_step);
    }

    generator
  }

  pub fn to_value_units(&self, source: ValueUnits, target: ValueUnits) -> Self {
    let temperature_scale = ValueUnits::scale_between(source, target, TEMPERATURE_U);

    Self {
      start_temperature: self.start_temperature * temperature_scale,
      end_temperature: self.end_temperature * temperature_scale,
      acceptance_distance: self.acceptance_distance.to_value_units(source, target),
      achieved_distance: self.achieved_distance.to_value_units(source, target),
      temperature_step: self.temperature_step.map(|step| step * temperature_scale),
      threshold: self.threshold.map(|threshold| threshold * temperature_scale),
      save_step: self.save_step,
    }
  }
}
