use crate::sim_core::world::utils::SparseManager;

pub struct GravityManager {
  inner: SparseManager<f64>,
}

impl GravityManager {
  pub fn constant(gravity: f64) -> Self {
    Self::new(vec![(0, gravity)])
  }

  pub fn new(gravities: Vec<(usize, f64)>) -> Self {
    GravityManager {
      inner: SparseManager::new(gravities, "GravityManager"),
    }
  }

  pub fn get_gravity(&mut self, iteration: usize) -> f64 {
    self.inner.get(iteration)
  }

  pub fn initial_gravity(&self) -> f64 {
    self.inner.initial_value()
  }

  pub fn scheduled_changes(&self) -> Vec<(usize, f64)> {
    self.inner.scheduled_changes()
  }
}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn constant_gravity_returns_same_value() {
    let mut manager = GravityManager::constant(1.5);
    assert_eq!(manager.get_gravity(0), 1.5);
    assert_eq!(manager.get_gravity(100), 1.5);
    assert_eq!(manager.get_gravity(1000), 1.5);
  }

  #[test]
  fn sparse_schedule_advances_in_order() {
    let mut manager = GravityManager::new(vec![(0, 1.0), (5, 0.5), (10, 0.0)]);

    assert_eq!(manager.get_gravity(0), 1.0);
    assert_eq!(manager.get_gravity(4), 1.0);
    assert_eq!(manager.get_gravity(5), 0.5);
    assert_eq!(manager.get_gravity(9), 0.5);
    assert_eq!(manager.get_gravity(10), 0.0);
    assert_eq!(manager.get_gravity(100), 0.0);
  }

  #[test]
  #[should_panic(expected = "iteration 10 was skipped")]
  fn skipping_iteration_panics() {
    let mut manager = GravityManager::new(vec![(0, 1.0), (5, 0.5)]);
    manager.get_gravity(10);
  }
}
