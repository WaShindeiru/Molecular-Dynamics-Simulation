use crate::persistence::dto::world::boxed::box_container::BoxContainerDTO;

pub struct HistoryDTO {
  pub box_container: Vec<BoxContainerDTO>,
  pub thermostat_epsilon: Vec<f64>,
}
