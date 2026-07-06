use nalgebra::Vector3;

use crate::data::config::particle_config::VelocityScheduleConfig;
use crate::persistence::json::particle_config::Vector3Record;

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct VelocityChangeEntry {
  pub iteration: usize,
  pub velocity: Vector3Record,
}

impl VelocityChangeEntry {
  pub fn from_runtime(iteration: usize, velocity: Vector3<f64>) -> Self {
    VelocityChangeEntry {
      iteration,
      velocity: Vector3Record::from(velocity),
    }
  }

  pub fn to_runtime(&self) -> (usize, Vector3<f64>) {
    (self.iteration, self.velocity.to_runtime())
  }

  pub fn scale(&self, velocity_scale: f64) -> Self {
    VelocityChangeEntry {
      iteration: self.iteration,
      velocity: self.velocity.scale(velocity_scale),
    }
  }
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct VelocityManagerFile {
  pub id: usize,
  pub velocities: Vec<VelocityChangeEntry>,
}

impl VelocityManagerFile {
  pub fn from_schedule(schedule: &VelocityScheduleConfig) -> Self {
    VelocityManagerFile {
      id: schedule.particle_velocity_manager_id,
      velocities: schedule
        .velocities
        .iter()
        .map(|(iteration, velocity)| VelocityChangeEntry::from_runtime(*iteration, *velocity))
        .collect(),
    }
  }

  pub fn scale_velocities(&self, velocity_scale: f64) -> Self {
    VelocityManagerFile {
      id: self.id,
      velocities: self.velocities.iter().map(|e| e.scale(velocity_scale)).collect(),
    }
  }

  pub fn to_schedule(&self) -> VelocityScheduleConfig {
    VelocityScheduleConfig {
      particle_velocity_manager_id: self.id,
      velocities: self.velocities.iter().map(|e| e.to_runtime()).collect(),
    }
  }
}
