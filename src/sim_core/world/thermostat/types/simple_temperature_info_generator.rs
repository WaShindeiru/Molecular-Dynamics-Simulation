use crate::data::TimeIterationDistance;

use super::{TemperatureInfo, TemperatureInfoGenerator};

const DEFAULT_TEMPERATURE_STEP: f64 = 70.;
const DEFAULT_THRESHOLD: f64 = 30.;
const DEFAULT_SAVE_STEP: usize = 3;

pub struct SimpleTemperatureInfoGenerator {
  start_temperature: Option<f64>,
  end_temperature: Option<f64>,
  acceptance_distance: Option<TimeIterationDistance>,
  achieved_distance: Option<TimeIterationDistance>,
  temperature_step: Option<f64>,
  threshold: Option<f64>,
  save_step: Option<usize>,
}

impl SimpleTemperatureInfoGenerator {
  pub fn new() -> Self {
    Self {
      start_temperature: None,
      end_temperature: None,
      acceptance_distance: None,
      achieved_distance: None,
      temperature_step: None,
      threshold: None,
      save_step: None,
    }
  }

  pub fn start_temperature(mut self, value: f64) -> Self {
    self.start_temperature = Some(value);
    self
  }

  pub fn end_temperature(mut self, value: f64) -> Self {
    self.end_temperature = Some(value);
    self
  }

  pub fn acceptance_distance(mut self, value: TimeIterationDistance) -> Self {
    self.acceptance_distance = Some(value);
    self
  }

  pub fn achieved_distance(mut self, value: TimeIterationDistance) -> Self {
    self.achieved_distance = Some(value);
    self
  }

  pub fn temperature_step(mut self, value: f64) -> Self {
    self.temperature_step = Some(value);
    self
  }

  pub fn threshold(mut self, value: f64) -> Self {
    self.threshold = Some(value);
    self
  }

  pub fn save_step(mut self, value: usize) -> Self {
    self.save_step = Some(value);
    self
  }
}

impl TemperatureInfoGenerator for SimpleTemperatureInfoGenerator {
  fn generate(&self) -> Result<Vec<TemperatureInfo>, String> {
    let start_temperature = self
      .start_temperature
      .ok_or("start_temperature must be set")?;
    let end_temperature = self.end_temperature.ok_or("end_temperature must be set")?;
    let acceptance_distance = self
      .acceptance_distance
      .ok_or("acceptance_distance must be set")?;
    let achieved_distance = self
      .achieved_distance
      .ok_or("achieved_distance must be set")?;

    let temperature_step = self.temperature_step.unwrap_or(DEFAULT_TEMPERATURE_STEP);
    let threshold = self.threshold.unwrap_or(DEFAULT_THRESHOLD);
    let save_step = self.save_step.unwrap_or(DEFAULT_SAVE_STEP);

    if save_step == 0 {
      return Err("save_step must be greater than 0".to_string());
    }

    let temperatures = build_temperature_series(
      start_temperature,
      end_temperature,
      temperature_step,
    )?;

    let len = temperatures.len();
    Ok(temperatures
      .into_iter()
      .enumerate()
      .map(|(index, desired_temperature)| TemperatureInfo {
        desired_temperature,
        acceptance_distance,
        achieved_distance,
        threshold,
        save: should_save(index, len, save_step),
      })
      .collect())
  }
}

fn build_temperature_series(
  start_temperature: f64,
  end_temperature: f64,
  temperature_step: f64,
) -> Result<Vec<f64>, String> {
  let delta = end_temperature - start_temperature;

  if delta == 0.0 {
    return Ok(vec![start_temperature]);
  }

  if temperature_step == 0.0 {
    return Err(
      "temperature_step must be non-zero when start_temperature and end_temperature differ"
        .to_string(),
    );
  }

  if delta.signum() != temperature_step.signum() {
    return Err(
      "temperature_step must have the same sign as (end_temperature - start_temperature)"
        .to_string(),
    );
  }

  let mut temperatures = vec![start_temperature];
  let mut current = start_temperature;

  if temperature_step > 0.0 {
    while current + temperature_step < end_temperature {
      current += temperature_step;
      temperatures.push(current);
    }
  } else {
    while current + temperature_step > end_temperature {
      current += temperature_step;
      temperatures.push(current);
    }
  }

  if !((temperatures.last().copied().unwrap() - Some(end_temperature).unwrap()).abs() < 2.0) {
    temperatures.push(end_temperature);
  }

  Ok(temperatures)
}

fn should_save(index: usize, len: usize, save_step: usize) -> bool {
  if index == 0 || index == len - 1 {
    return true;
  }

  index % save_step == 0
}

#[cfg(test)]
mod tests {
  use super::*;

  fn base_generator() -> SimpleTemperatureInfoGenerator {
    SimpleTemperatureInfoGenerator::new()
      .start_temperature(2000.)
      .end_temperature(200.)
      .acceptance_distance(TimeIterationDistance::Time { value: 1.0 })
      .achieved_distance(TimeIterationDistance::Iteration { value: 100 })
  }

  #[test]
  fn generates_cooling_series_ending_exactly_at_end_temperature() {
    let infos = base_generator()
      .temperature_step(-70.)
      .generate()
      .unwrap();

    assert_eq!(infos.first().unwrap().desired_temperature, 2000.);
    assert_eq!(infos.last().unwrap().desired_temperature, 200.);

    for info in &infos {
      assert_eq!(info.acceptance_distance, TimeIterationDistance::Time { value: 1.0 });
      assert_eq!(
        info.achieved_distance,
        TimeIterationDistance::Iteration { value: 100 }
      );
      assert_eq!(info.threshold, DEFAULT_THRESHOLD);
    }

    let temps: Vec<f64> = infos.iter().map(|i| i.desired_temperature).collect();
    assert!(temps.windows(2).all(|w| (w[0] - w[1]).abs() <= 70.));
    assert!(temps.contains(&200.));
  }

  #[test]
  fn generates_heating_series() {
    let infos = SimpleTemperatureInfoGenerator::new()
      .start_temperature(200.)
      .end_temperature(500.)
      .temperature_step(100.)
      .acceptance_distance(TimeIterationDistance::Iteration { value: 0 })
      .achieved_distance(TimeIterationDistance::Iteration { value: 0 })
      .generate()
      .unwrap();

    let temps: Vec<f64> = infos.iter().map(|i| i.desired_temperature).collect();
    assert_eq!(temps, vec![200., 300., 400., 500.]);
  }

  #[test]
  fn single_entry_when_start_equals_end() {
    let infos = SimpleTemperatureInfoGenerator::new()
      .start_temperature(300.)
      .end_temperature(300.)
      .acceptance_distance(TimeIterationDistance::Iteration { value: 0 })
      .achieved_distance(TimeIterationDistance::Iteration { value: 0 })
      .generate()
      .unwrap();

    assert_eq!(infos.len(), 1);
    assert_eq!(infos[0].desired_temperature, 300.);
    assert!(infos[0].save);
  }

  #[test]
  fn save_step_marks_every_nth_entry_and_always_start_and_end() {
    let infos = SimpleTemperatureInfoGenerator::new()
      .start_temperature(0.)
      .end_temperature(500.)
      .temperature_step(100.)
      .save_step(2)
      .acceptance_distance(TimeIterationDistance::Iteration { value: 0 })
      .achieved_distance(TimeIterationDistance::Iteration { value: 0 })
      .generate()
      .unwrap();

    let save_flags: Vec<bool> = infos.iter().map(|i| i.save).collect();
    assert_eq!(save_flags, vec![true, false, true, false, true, true]);
  }

  #[test]
  fn default_optional_fields_are_applied() {
    let infos = SimpleTemperatureInfoGenerator::new()
      .start_temperature(1000.)
      .end_temperature(1000.)
      .acceptance_distance(TimeIterationDistance::Iteration { value: 0 })
      .achieved_distance(TimeIterationDistance::Iteration { value: 0 })
      .generate()
      .unwrap();

    assert_eq!(infos[0].threshold, DEFAULT_THRESHOLD);
    assert!(infos[0].save);
  }

  #[test]
  fn rejects_missing_required_fields() {
    let err = SimpleTemperatureInfoGenerator::new()
      .end_temperature(200.)
      .acceptance_distance(TimeIterationDistance::Iteration { value: 0 })
      .achieved_distance(TimeIterationDistance::Iteration { value: 0 })
      .generate()
      .unwrap_err();
    assert_eq!(err, "start_temperature must be set");
  }

  #[test]
  fn rejects_opposite_sign_step() {
    let err = base_generator()
      .temperature_step(70.)
      .generate()
      .unwrap_err();
    assert_eq!(
      err,
      "temperature_step must have the same sign as (end_temperature - start_temperature)"
    );
  }

  #[test]
  fn rejects_zero_step_when_range_is_non_zero() {
    let err = base_generator()
      .temperature_step(0.)
      .generate()
      .unwrap_err();
    assert_eq!(
      err,
      "temperature_step must be non-zero when start_temperature and end_temperature differ"
    );
  }
}
