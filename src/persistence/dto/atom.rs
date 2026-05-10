#[derive(Clone, serde::Serialize, serde::Deserialize)]
pub struct AtomDTO {
  pub id: usize,
  pub iteration: usize,
  pub atom_type: u64,
  pub x: f64,
  pub y: f64,
  pub z: f64,

  pub kinetic_energy: f64,
  pub potential_energy: f64,
  pub thermostat_work: f64,
  pub potential_gravity_energy: f64,

  pub force_x: f64,
  pub force_y: f64,
  pub force_z: f64,
}
