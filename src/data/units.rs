pub const ENERGY_U: f64 = 1.602e-19;
pub const TIME_U: f64 = 10.18e-15; // seconds
pub const R_U: f64 = 1.0e-10; // m
pub const K_B: f64 = 1.;
pub const TEMPERATURE_U: f64 = 11608.7; // K
pub const MASS_U: f64 = 1.66e-27; // kg
pub const VELOCITY_U: f64 = 9823.75; // m/s

#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ValueUnits {
  Unitless,
  Si,
}

impl Default for ValueUnits {
  fn default() -> Self {
    ValueUnits::Si
  }
}

impl ValueUnits {
  pub fn scale_between(source: ValueUnits, target: ValueUnits, base_si_unit: f64) -> f64 {
    match (source, target) {
      (ValueUnits::Unitless, ValueUnits::Si) => base_si_unit,
      (ValueUnits::Si, ValueUnits::Unitless) => 1.0 / base_si_unit,
      _ => 1.0,
    }
  }
}
