use std::collections::HashMap;
use std::sync::{Arc, RwLock};
use nalgebra::Vector3;
use crate::particle::Particle;
use crate::sim_core::world::boxed_world::box_container::BoxContainer;

mod handle_task;
pub mod threads;

pub enum BoxTask {
  VelocityTask {
    box_container: Arc<RwLock<BoxContainer>>,
    box_id: usize,
    time_step: f64,
    previous_thermostat_epsilon: f64,
    current_iteration: usize,
  },
  ForceTask,
}

pub enum BoxResult {
  VelocityResult {
    half_velocity_cache: HashMap<usize, Vector3<f64>>,
    new_position_atoms: HashMap<usize, Particle>,
  },
  ForceResult,
}