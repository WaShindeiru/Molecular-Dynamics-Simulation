use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap};

use nalgebra::Vector3;

use crate::data::ParticleConfig;
use crate::particle::Particle;

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
  // Indexed by particle_velocity_manager_id. Particles that share a velocity schedule
  // share a single entry here instead of duplicating one ParticleVelocityManager each.
  managers: Vec<ParticleVelocityManager>,
  // Indexed by particle id (dense over all particles, not just CustomVelocityAtoms),
  // for O(1) lookup of which manager (if any) drives a given particle's velocity.
  particle_id_to_manager_id: Vec<Option<usize>>,
}

impl VelocityManager {
  pub fn new(
    managers: Vec<ParticleVelocityManager>,
    particle_id_to_manager_id: Vec<Option<usize>>,
  ) -> Self {
    VelocityManager {
      managers,
      particle_id_to_manager_id,
    }
  }

  pub fn empty() -> Self {
    VelocityManager {
      managers: Vec::new(),
      particle_id_to_manager_id: Vec::new(),
    }
  }

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

    let mut particle_id_to_manager_id: Vec<Option<usize>> =
      vec![None; particle_config.atoms.len()];
    for atom in &particle_config.atoms {
      if let Particle::CustomVelocityAtom(p) = atom {
        particle_id_to_manager_id[p.get_id()] = Some(p.get_particle_velocity_manager_id());
      }
    }

    VelocityManager {
      managers,
      particle_id_to_manager_id,
    }
  }

  /// Returns the velocity for each particle at the given iteration, indexed by particle id
  /// (particles with no velocity manager are absent from the map), advancing the internal
  /// state of each ParticleVelocityManager as needed.
  pub fn compute_velocities_for_iteration(&mut self, iteration: usize) -> HashMap<usize, Vector3<f64>> {
    let manager_velocities: Vec<Vector3<f64>> = self
      .managers
      .iter_mut()
      .map(|manager| manager.get_velocity(iteration))
      .collect();

    self
      .particle_id_to_manager_id
      .iter()
      .enumerate()
      .filter_map(|(particle_id, manager_id)| {
        manager_id.map(|id| (particle_id, manager_velocities[id]))
      })
      .collect()
  }
}
