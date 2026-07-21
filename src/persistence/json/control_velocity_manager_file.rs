use nalgebra::Vector3;

use crate::data::config::particle_config::ControlVelocityScheduleConfig;
use crate::persistence::json::particle_config::Vector3Record;

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct ControlVelocityChangeEntry {
  pub iteration: usize,
  pub component_velocity: Vector3Record,
  pub desired_velocity: Vector3Record,
}

impl ControlVelocityChangeEntry {
  pub fn from_runtime(iteration: usize, component_velocity: Vector3<f64>, desired_velocity: Vector3<f64>) -> Self {
    ControlVelocityChangeEntry {
      iteration,
      component_velocity: Vector3Record::from(component_velocity),
      desired_velocity: Vector3Record::from(desired_velocity),
    }
  }

  pub fn to_runtime(&self) -> (usize, Vector3<f64>, Vector3<f64>) {
    (self.iteration, self.component_velocity.to_runtime(), self.desired_velocity.to_runtime())
  }

  pub fn scale(&self, velocity_scale: f64) -> Self {
    ControlVelocityChangeEntry {
      iteration: self.iteration,
      component_velocity: self.component_velocity.scale(velocity_scale),
      desired_velocity: self.desired_velocity.scale(velocity_scale),
    }
  }
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct ControlVelocityManagerFile {
  pub id: usize,
  pub changes: Vec<ControlVelocityChangeEntry>,
}

impl ControlVelocityManagerFile {
  pub fn from_schedule(schedule: &ControlVelocityScheduleConfig) -> Self {
    ControlVelocityManagerFile {
      id: schedule.control_velocity_manager_id,
      changes: schedule
        .entries
        .iter()
        .map(|(iteration, component_velocity, desired_velocity)| {
          ControlVelocityChangeEntry::from_runtime(*iteration, *component_velocity, *desired_velocity)
        })
        .collect(),
    }
  }

  pub fn scale_velocities(&self, velocity_scale: f64) -> Self {
    ControlVelocityManagerFile {
      id: self.id,
      changes: self.changes.iter().map(|e| e.scale(velocity_scale)).collect(),
    }
  }

  pub fn to_schedule(&self) -> ControlVelocityScheduleConfig {
    ControlVelocityScheduleConfig {
      control_velocity_manager_id: self.id,
      entries: self.changes.iter().map(|e| e.to_runtime()).collect(),
    }
  }
}
