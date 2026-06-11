use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap};

use log::info;

pub struct SparseManager<T: Copy + std::fmt::Debug> {
  type_name: &'static str,
  current_value: T,
  value_list: HashMap<usize, T>,
  value_heap: BinaryHeap<Reverse<usize>>,
}

impl<T: Copy + std::fmt::Debug> SparseManager<T> {
  pub fn new(entries: Vec<(usize, T)>, type_name: &'static str) -> Self {
    assert!(
      !entries.is_empty(),
      "{} requires at least one entry",
      type_name
    );
    assert!(
      entries.iter().any(|(i, _)| *i == 0),
      "{} requires an entry for iteration 0",
      type_name
    );

    let value_list: HashMap<usize, T> = entries.into_iter().collect();
    let initial_value = *value_list.get(&0).unwrap();

    let value_heap: BinaryHeap<Reverse<usize>> = value_list
      .keys()
      .filter(|&&k| k != 0)
      .map(|k| Reverse(*k))
      .collect();

    SparseManager {
      type_name,
      current_value: initial_value,
      value_list,
      value_heap,
    }
  }

  pub fn constant(value: T, type_name: &'static str) -> Self {
    Self::new(vec![(0, value)], type_name)
  }

  pub fn get(&mut self, iteration: usize) -> T {
    loop {
      match self.value_heap.peek() {
        None => return self.current_value,
        Some(Reverse(top)) if *top > iteration => return self.current_value,
        Some(Reverse(top)) if *top == iteration => {
          let top = *top;
          self.value_heap.pop();
          self.current_value = *self.value_list.get(&top).unwrap();
          info!("{} changed to: {:?}", self.type_name, self.current_value);
          return self.current_value;
        }
        Some(Reverse(top)) => panic!(
          "{}: iteration {} was skipped (heap top: {}). \
           Iterations must be consumed in order.",
          self.type_name, iteration, top
        ),
      }
    }
  }

  pub fn initial_value(&self) -> T {
    *self
      .value_list
      .get(&0)
      .expect("SparseManager requires an entry for iteration 0")
  }

  pub fn scheduled_changes(&self) -> Vec<(usize, T)> {
    let mut entries: Vec<_> = self
      .value_list
      .iter()
      .map(|(&iteration, &value)| (iteration, value))
      .collect();
    entries.sort_unstable_by_key(|(iteration, _)| *iteration);
    entries
  }
}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn constant_returns_same_value() {
    let mut m: SparseManager<f64> = SparseManager::constant(1.5, "Test");
    assert_eq!(m.get(0), 1.5);
    assert_eq!(m.get(100), 1.5);
    assert_eq!(m.get(1000), 1.5);
  }

  #[test]
  fn sparse_schedule_advances_in_order() {
    let mut m: SparseManager<f64> =
      SparseManager::new(vec![(0, 1.0), (5, 0.5), (10, 0.0)], "Test");

    assert_eq!(m.get(0), 1.0);
    assert_eq!(m.get(4), 1.0);
    assert_eq!(m.get(5), 0.5);
    assert_eq!(m.get(9), 0.5);
    assert_eq!(m.get(10), 0.0);
    assert_eq!(m.get(100), 0.0);
  }

  #[test]
  #[should_panic(expected = "iteration 10 was skipped")]
  fn skipping_iteration_panics() {
    let mut m: SparseManager<f64> = SparseManager::new(vec![(0, 1.0), (5, 0.5)], "Test");
    m.get(10);
  }

  #[test]
  fn initial_value_returns_iteration_zero_entry() {
    let m: SparseManager<f64> = SparseManager::new(vec![(0, 3.14), (10, 2.71)], "Test");
    assert_eq!(m.initial_value(), 3.14);
  }

  #[test]
  fn scheduled_changes_sorted() {
    let m: SparseManager<f64> = SparseManager::new(vec![(10, 2.0), (0, 1.0), (5, 1.5)], "Test");
    let changes = m.scheduled_changes();
    assert_eq!(changes, vec![(0, 1.0), (5, 1.5), (10, 2.0)]);
  }
}
