use crate::sim_core::world::utils::SparseManager;

pub struct TimestepManager {
  inner: SparseManager<f64>,
}

impl TimestepManager {
  pub fn constant(timestep: f64) -> Self {
    Self::new(vec![(0, timestep)])
  }

  pub fn new(timesteps: Vec<(usize, f64)>) -> Self {
    TimestepManager {
      inner: SparseManager::new(timesteps, "TimestepManager"),
    }
  }

  pub fn get_timestep(&mut self, iteration: usize) -> f64 {
    self.inner.get(iteration)
  }

  pub fn initial_timestep(&self) -> f64 {
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
  fn constant_timestep_returns_same_value() {
    let mut manager = TimestepManager::constant(1e-15);
    assert_eq!(manager.get_timestep(0), 1e-15);
    assert_eq!(manager.get_timestep(100), 1e-15);
    assert_eq!(manager.get_timestep(1000), 1e-15);
  }

  #[test]
  fn sparse_schedule_advances_in_order() {
    let mut manager =
      TimestepManager::new(vec![(0, 1e-15), (100, 5e-16), (200, 1e-16)]);

    assert_eq!(manager.get_timestep(0), 1e-15);
    assert_eq!(manager.get_timestep(99), 1e-15);
    assert_eq!(manager.get_timestep(100), 5e-16);
    assert_eq!(manager.get_timestep(199), 5e-16);
    assert_eq!(manager.get_timestep(200), 1e-16);
    assert_eq!(manager.get_timestep(500), 1e-16);
  }

  #[test]
  #[should_panic(expected = "iteration 200 was skipped")]
  fn skipping_iteration_panics() {
    let mut manager = TimestepManager::new(vec![(0, 1e-15), (100, 5e-16)]);
    manager.get_timestep(200);
  }
}
