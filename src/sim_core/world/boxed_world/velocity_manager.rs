use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap};

use nalgebra::Vector3;

use crate::data::config::particle_config::VelocityScheduleConfig;

pub struct ParticleVelocityManager {
  current_velocity: Vector3<f64>,
  velocity_list: HashMap<usize, Vector3<f64>>,
  velocity_heap: BinaryHeap<Reverse<usize>>,
}

impl ParticleVelocityManager {
  pub fn new(velocities: Vec<(usize, Vector3<f64>)>) -> Self {
    assert!(
      !velocities.is_empty(),
      "ParticleVelocityManager requires at least one velocity entry"
    );
    assert!(
      velocities.iter().any(|(i, _)| *i == 0),
      "ParticleVelocityManager requires an entry for iteration 0 (initial velocity)"
    );

    let velocity_list: HashMap<usize, Vector3<f64>> = velocities.into_iter().collect();
    let initial_velocity = *velocity_list.get(&0).unwrap();

    // Iteration 0 is consumed immediately as the initial state — it must not sit in
    // the heap, otherwise the first real call (e.g. get_velocity(1)) would see
    // heap-top 0 < 1 and panic with "iteration skipped".
    let velocity_heap: BinaryHeap<Reverse<usize>> = velocity_list
      .keys()
      .filter(|&&k| k != 0)
      .map(|k| Reverse(*k))
      .collect();

    ParticleVelocityManager {
      current_velocity: initial_velocity,
      velocity_list,
      velocity_heap,
    }
  }

  /// Returns the velocity for the given iteration, advancing the internal state if needed.
  /// Panics if `iteration` skips over a scheduled change (heap top < iteration).
  pub fn get_velocity(&mut self, iteration: usize) -> Vector3<f64> {
    loop {
      match self.velocity_heap.peek() {
        None => return self.current_velocity,
        Some(Reverse(top)) if *top > iteration => return self.current_velocity,
        Some(Reverse(top)) if *top == iteration => {
          let top = *top;
          self.velocity_heap.pop();
          self.current_velocity = *self.velocity_list.get(&top).unwrap();
          return self.current_velocity;
        }
        Some(Reverse(top)) => panic!(
          "VelocityManager: iteration {} was skipped (heap top: {}). \
           Iterations must be consumed in order.",
          iteration, top
        ),
      }
    }
  }
}

pub struct VelocityManager {
  managers: HashMap<usize, ParticleVelocityManager>,
}

impl VelocityManager {
  pub fn new(managers: HashMap<usize, ParticleVelocityManager>) -> Self {
    VelocityManager { managers }
  }

  pub fn empty() -> Self {
    VelocityManager {
      managers: HashMap::new(),
    }
  }

  pub fn from_schedules(schedules: &[VelocityScheduleConfig]) -> Self {
    let managers = schedules
      .iter()
      .map(|s| {
        let manager = ParticleVelocityManager::new(s.velocities.clone());
        (s.particle_id, manager)
      })
      .collect();
    VelocityManager { managers }
  }

  /// Returns the velocity for each managed particle at the given iteration,
  /// advancing the internal state of each ParticleVelocityManager as needed.
  pub fn compute_velocities_for_iteration(
    &mut self,
    iteration: usize,
  ) -> HashMap<usize, Vector3<f64>> {
    self
      .managers
      .iter_mut()
      .map(|(id, manager)| (*id, manager.get_velocity(iteration)))
      .collect()
  }

  pub fn get_velocity_for_iteration_particle(
    &mut self,
    iteration: usize,
    particle_id: usize,
  ) -> Vector3<f64> {
    self
      .managers
      .get_mut(&particle_id)
      .unwrap()
      .get_velocity(iteration);

    Vector3::zeros()
  }
}
