use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap};

use log::info;

pub struct GravityManager {
  current_gravity: f64,
  gravity_list: HashMap<usize, f64>,
  gravity_heap: BinaryHeap<Reverse<usize>>,
}

impl GravityManager {
  pub fn constant(gravity: f64) -> Self {
    Self::new(vec![(0, gravity)])
  }

  pub fn new(gravities: Vec<(usize, f64)>) -> Self {
    assert!(
      !gravities.is_empty(),
      "GravityManager requires at least one gravity entry"
    );
    assert!(
      gravities.iter().any(|(i, _)| *i == 0),
      "GravityManager requires an entry for iteration 0 (initial gravity)"
    );

    let gravity_list: HashMap<usize, f64> = gravities.into_iter().collect();
    let initial_gravity = *gravity_list.get(&0).unwrap();

    // Iteration 0 is consumed immediately as the initial state — it must not sit in
    // the heap, otherwise the first real call (e.g. get_gravity(1)) would see
    // heap-top 0 < 1 and panic with "iteration skipped".
    let gravity_heap: BinaryHeap<Reverse<usize>> = gravity_list
      .keys()
      .filter(|&&k| k != 0)
      .map(|k| Reverse(*k))
      .collect();

    GravityManager {
      current_gravity: initial_gravity,
      gravity_list,
      gravity_heap,
    }
  }

  /// Returns the gravity for the given iteration, advancing the internal state if needed.
  /// Panics if `iteration` skips over a scheduled change (heap top < iteration).
  pub fn get_gravity(&mut self, iteration: usize) -> f64 {
    loop {
      match self.gravity_heap.peek() {
        None => return self.current_gravity,
        Some(Reverse(top)) if *top > iteration => return self.current_gravity,
        Some(Reverse(top)) if *top == iteration => {
          let top = *top;
          self.gravity_heap.pop();
          self.current_gravity = *self.gravity_list.get(&top).unwrap();

          info!("Gravity changed to: {}", self.current_gravity);

          return self.current_gravity;
        }
        Some(Reverse(top)) => panic!(
          "GravityManager: iteration {} was skipped (heap top: {}). \
           Iterations must be consumed in order.",
          iteration, top
        ),
      }
    }
  }

  pub fn initial_gravity(&self) -> f64 {
    *self
      .gravity_list
      .get(&0)
      .expect("GravityManager requires an entry for iteration 0")
  }

  pub fn scheduled_changes(&self) -> Vec<(usize, f64)> {
    let mut entries: Vec<_> = self
      .gravity_list
      .iter()
      .map(|(&iteration, &gravity)| (iteration, gravity))
      .collect();
    entries.sort_unstable_by_key(|(iteration, _)| *iteration);
    entries
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
