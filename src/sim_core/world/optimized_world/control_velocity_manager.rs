use nalgebra::Vector3;

use crate::data::ParticleConfig;
use crate::particle::Particle;
use crate::sim_core::world::boxed_world::velocity_manager::{GenericParticleManager, GenericVelocityManager};

/// Per-iteration output is (component_velocity, desired_velocity) for a VelocityControlledParticle.
pub type ControlVelocityManager = GenericVelocityManager<(Vector3<f64>, Vector3<f64>)>;

impl ControlVelocityManager {
  pub fn from_config(particle_config: &ParticleConfig) -> Self {
    let manager_count = particle_config
      .control_velocity_schedules
      .iter()
      .map(|s| s.control_velocity_manager_id)
      .max()
      .map_or(0, |id| id + 1);

    let mut manager_slots: Vec<Option<GenericParticleManager<(Vector3<f64>, Vector3<f64>)>>> =
      (0..manager_count).map(|_| None).collect();
    for schedule in &particle_config.control_velocity_schedules {
      let entries = schedule
        .entries
        .iter()
        .map(|&(iteration, component_velocity, desired_velocity)| {
          (iteration, (component_velocity, desired_velocity))
        })
        .collect();
      manager_slots[schedule.control_velocity_manager_id] =
        Some(GenericParticleManager::new(entries));
    }
    let managers: Vec<GenericParticleManager<(Vector3<f64>, Vector3<f64>)>> = manager_slots
      .into_iter()
      .enumerate()
      .map(|(id, manager)| {
        manager.unwrap_or_else(|| {
          panic!(
            "ControlVelocityManager: no control velocity schedule provided for control_velocity_manager_id {}",
            id
          )
        })
      })
      .collect();

    let mut particle_id_to_manager_id: Vec<(usize, usize)> = Vec::new();
    for atom in &particle_config.atoms {
      if let Particle::VelocityControlledParticle(p) = atom {
        particle_id_to_manager_id.push((p.get_id(), p.get_control_velocity_manager_id()));
      }
    }

    ControlVelocityManager::new(managers, particle_id_to_manager_id)
  }

  /// Returns (particle_id, (component_velocity, desired_velocity)) pairs for every
  /// VelocityControlledParticle, advancing the internal state of each manager as needed.
  pub fn compute_controlled_velocities_for_iteration(
    &mut self,
    iteration: usize,
  ) -> &[(usize, (Vector3<f64>, Vector3<f64>))] {
    self.compute_for_iteration(iteration)
  }
}
