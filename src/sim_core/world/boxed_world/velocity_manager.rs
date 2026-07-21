use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap};

use nalgebra::Vector3;

use crate::data::ParticleConfig;
use crate::particle::Particle;

pub struct GenericParticleManager<T> {
  current_value: T,
  value_list: HashMap<usize, T>,
  value_heap: BinaryHeap<Reverse<usize>>,
}

impl<T: Copy> GenericParticleManager<T> {
  pub fn new(values: Vec<(usize, T)>) -> Self {
    assert!(
      !values.is_empty(),
      "ParticleVelocityManager requires at least one velocity entry"
    );
    assert!(
      values.iter().any(|(i, _)| *i == 0),
      "ParticleVelocityManager requires an entry for iteration 0 (initial velocity)"
    );

    let value_list: HashMap<usize, T> = values.into_iter().collect();
    let initial_value = *value_list.get(&0).unwrap();

    // Iteration 0 is consumed immediately as the initial state — it must not sit in
    // the heap, otherwise the first real call (e.g. get_value(1)) would see
    // heap-top 0 < 1 and panic with "iteration skipped".
    let value_heap: BinaryHeap<Reverse<usize>> = value_list
      .keys()
      .filter(|&&k| k != 0)
      .map(|k| Reverse(*k))
      .collect();

    GenericParticleManager {
      current_value: initial_value,
      value_list,
      value_heap,
    }
  }

  /// Returns the value for the given iteration, advancing the internal state if needed.
  /// Panics if `iteration` skips over a scheduled change (heap top < iteration).
  pub fn get_value(&mut self, iteration: usize) -> T {
    loop {
      match self.value_heap.peek() {
        None => return self.current_value,
        Some(Reverse(top)) if *top > iteration => return self.current_value,
        Some(Reverse(top)) if *top == iteration => {
          let top = *top;
          self.value_heap.pop();
          self.current_value = *self.value_list.get(&top).unwrap();
          return self.current_value;
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

pub struct GenericVelocityManager<T> {
  // Indexed by particle_velocity_manager_id. Particles that share a velocity schedule
  // share a single entry here instead of duplicating one GenericParticleManager each.
  managers: Vec<GenericParticleManager<T>>,
  // Sparse (particle_id, manager_id) pairs, one per driven particle only —
  // most particles have no manager and would just waste a slot in a dense mapping.
  particle_id_to_manager_id: Vec<(usize, usize)>,
  // Manager outputs from the last compute_for_iteration call, used to detect
  // whether the per-particle result actually needs to be rebuilt (schedules
  // change rarely, so most iterations hit this cache).
  cached_manager_values: Vec<T>,
  cached_particle_values: Vec<(usize, T)>,
}

impl<T: Copy + PartialEq> GenericVelocityManager<T> {
  pub fn new(
    managers: Vec<GenericParticleManager<T>>,
    particle_id_to_manager_id: Vec<(usize, usize)>,
  ) -> Self {
    GenericVelocityManager {
      managers,
      particle_id_to_manager_id,
      cached_manager_values: Vec::new(),
      cached_particle_values: Vec::new(),
    }
  }

  pub fn empty() -> Self {
    GenericVelocityManager {
      managers: Vec::new(),
      particle_id_to_manager_id: Vec::new(),
      cached_manager_values: Vec::new(),
      cached_particle_values: Vec::new(),
    }
  }

  /// Returns (particle_id, value) pairs for every particle driven by a manager, advancing
  /// the internal state of each GenericParticleManager as needed. Callers should iterate the
  /// slice directly rather than looking entries up by id.
  ///
  /// Schedules change sparsely, so most calls leave every manager's value unchanged; in that
  /// case the previously built list is reused instead of rescanning particle_id_to_manager_id.
  pub(crate) fn compute_for_iteration(&mut self, iteration: usize) -> &[(usize, T)] {
    let manager_values: Vec<T> = self.managers.iter_mut().map(|manager| manager.get_value(iteration)).collect();

    if manager_values != self.cached_manager_values {
      self.cached_particle_values.clear();
      self
        .cached_particle_values
        .extend(self.particle_id_to_manager_id.iter().map(
          |&(particle_id, manager_id)| (particle_id, manager_values[manager_id]),
        ));
      self.cached_manager_values = manager_values;
    }

    &self.cached_particle_values
  }
}

pub type ParticleVelocityManager = GenericParticleManager<Vector3<f64>>;
pub type VelocityManager = GenericVelocityManager<Vector3<f64>>;

impl VelocityManager {
  pub fn from_config(particle_config: &ParticleConfig) -> Self {
    let manager_count = particle_config
      .velocity_schedules
      .iter()
      .map(|s| s.particle_velocity_manager_id)
      .max()
      .map_or(0, |id| id + 1);

    let mut manager_slots: Vec<Option<ParticleVelocityManager>> =
      (0..manager_count).map(|_| None).collect();
    for schedule in &particle_config.velocity_schedules {
      manager_slots[schedule.particle_velocity_manager_id] =
        Some(ParticleVelocityManager::new(schedule.velocities.clone()));
    }
    let managers: Vec<ParticleVelocityManager> = manager_slots
      .into_iter()
      .enumerate()
      .map(|(id, manager)| {
        manager.unwrap_or_else(|| {
          panic!(
            "VelocityManager: no velocity schedule provided for particle_velocity_manager_id {}",
            id
          )
        })
      })
      .collect();

    let mut particle_id_to_manager_id: Vec<(usize, usize)> = Vec::new();
    for atom in &particle_config.atoms {
      if let Particle::CustomVelocityAtom(p) = atom {
        particle_id_to_manager_id.push((p.get_id(), p.get_particle_velocity_manager_id()));
      }
    }

    VelocityManager {
      managers,
      particle_id_to_manager_id,
      cached_manager_values: Vec::new(),
      cached_particle_values: Vec::new(),
    }
  }

  /// Returns (particle_id, velocity) pairs for every particle driven by a velocity manager,
  /// advancing the internal state of each ParticleVelocityManager as needed. Callers should
  /// iterate the slice directly rather than looking entries up by id.
  pub fn compute_velocities_for_iteration(&mut self, iteration: usize) -> &[(usize, Vector3<f64>)] {
    self.compute_for_iteration(iteration)
  }
}
