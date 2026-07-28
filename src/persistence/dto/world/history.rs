use crate::persistence::dto::world::boxed::box_container::BoxContainerDTO;

pub struct HistoryDTO {
  pub box_container: Vec<BoxContainerDTO>,
  pub thermostat_epsilon: Vec<f64>,
  pub temperature: Vec<f64>,
  /// `Some` only for `OptimizedWorld` when the optional nanotube thermostat is configured.
  pub nanotube_thermostat_epsilon: Option<Vec<f64>>,
  pub nanotube_temperature: Option<Vec<f64>>,
}
