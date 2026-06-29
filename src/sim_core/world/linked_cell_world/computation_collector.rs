use std::collections::HashMap;

use nalgebra::Vector3;

use crate::data::SimulationConfig;
use crate::data::units::K_B;
use crate::particle::Particle;
use crate::sim_core::world::boundary_constraint::{EdgeCondition, ParticleCompliance};
use crate::sim_core::world::boxed_world::box_task::ForceTaskParticleData;
use crate::sim_core::world::computation::{compute_new_velocity};
use crate::sim_core::world::linked_cell_world::LinkedCellContainer;

pub struct ComputationCollector {
  config: SimulationConfig,
  half_velocity_cache: Vec<Vector3<f64>>,
  particle_compliance: Vec<ParticleCompliance>,
  local_container: LinkedCellContainer,
  force_visited: Vec<usize>,
}

impl ComputationCollector {
  pub fn from_integration_cache(
    config: SimulationConfig,
    num_of_particles: usize,
    half_velocity_cache: Vec<Vector3<f64>>,
    particle_compliance: Vec<ParticleCompliance>,
    local_container: LinkedCellContainer,
  ) -> Self {
    ComputationCollector {
      config,
      half_velocity_cache,
      particle_compliance,
      local_container,
      force_visited: vec![0; num_of_particles],
    }
  }

  pub fn apply_force_results<'a>(
    &mut self,
    force_results: impl IntoIterator<Item = (&'a usize, &'a ForceTaskParticleData)>,
  ) {
    for (particle_id, force_data) in force_results {
      let particle = self.local_container.particle_mut(*particle_id);
      let acceleration = force_data.force / particle.get_mass();
      particle.set_force(particle.get_force() + force_data.force);
      particle.set_potential_energy(particle.get_potential_energy() + force_data.potential_energy);
      particle.set_acceleration(particle.get_acceleration() + acceleration);
      self.force_visited[*particle_id] += 1;
    }
  }

  pub fn apply_gravity(&mut self, iteration: usize) {
    let potential_gravity_max = self
      .config
      .gravity_manager
      .lock()
      .expect("gravity manager lock poisoned")
      .get_gravity(iteration);
    let z_max = self.local_container.config().world_size.z;

    for id in 0..self.local_container.particles().len() {
      if self.local_container.particles()[id].is_none() {
        continue;
      }
      let particle = self.local_container.particle_mut(id);
      let new_force = particle.get_force()
        - Vector3::new(0., 0., 1.) * potential_gravity_max * particle.get_mass() / z_max;
      particle.set_force(new_force);
      particle.set_acceleration(new_force / particle.get_mass());
      particle.set_potential_gravity_energy(
        potential_gravity_max * particle.get_mass() * particle.get_position().z / z_max,
      );
    }
  }

  pub fn set_velocity(&mut self, thermostat_epsilon: f64) {
    let time_step = self.config.time_step;
    let edge_condition = self.config.edge_condition;
    let subtask_size = match edge_condition {
      EdgeCondition::Periodic { trigger_small_subtask_size, .. }
      | EdgeCondition::Simple { trigger_small_subtask_size, .. } => trigger_small_subtask_size,
      EdgeCondition::PeriodicAll => 1,
    };
    let collision_split = edge_condition.collision_split_enabled();

    for id in 0..self.local_container.particles().len() {
      if self.local_container.particles()[id].is_none() {
        continue;
      }
      if self.local_container.particles()[id].as_ref().unwrap().is_custom_velocity_atom() {
        continue;
      }
      let half_velocity = self.half_velocity_cache[id];
      let compliance = &self.particle_compliance[id];
      let effective_time_step = if !compliance.compliant && collision_split && subtask_size > 1 {
        time_step / subtask_size as f64
      } else {
        time_step
      };
      let particle = self.local_container.particle_mut(id);
      let new_velocity = compute_new_velocity(
        half_velocity,
        *particle.get_acceleration(),
        thermostat_epsilon,
        effective_time_step,
      );
      particle.set_velocity(new_velocity);
    }
  }

  pub fn apply_custom_velocities(&mut self, custom_velocities: &HashMap<usize, Vector3<f64>>) {
    for id in 0..self.local_container.particles().len() {
      if self.local_container.particles()[id].is_none() {
        continue;
      }
      if let Some(&vel) = custom_velocities.get(&id) {
        if let Particle::CustomVelocityAtom(p) = self.local_container.particle_mut(id) {
          p.set_velocity(vel);
        }
      }
    }
  }

  pub fn get_mean_temperature(&self) -> f64 {
    let relevant: Vec<_> = self
      .local_container
      .particles()
      .iter()
      .filter_map(|opt| opt.as_ref())
      .filter(|p| !p.is_custom_velocity_atom())
      .collect();

    let count = relevant.len();
    if count == 0 {
      return 0.0;
    }

    let temperature: f64 = relevant
      .iter()
      .map(|p| p.get_mass() * p.get_velocity().magnitude().powi(2) / (3. * K_B))
      .sum();
    temperature / count as f64
  }

  pub fn particles(&self) -> impl Iterator<Item = &Particle> {
    self.local_container.particles()
      .iter()
      .filter_map(|opt| opt.as_ref())
      .map(|arc| arc.as_ref())
  }

  #[cfg(debug_assertions)]
  pub fn assert_zero_net_force(&self) {
    let sum = self.local_container.particles().iter()
      .filter_map(|opt| opt.as_ref())
      .fold(nalgebra::Vector3::<f64>::zeros(), |acc, p| acc + p.get_force());
    let magnitude = sum.magnitude();
    debug_assert!(
      magnitude < 1e-10,
      "Newton's 3rd law violated: net force magnitude = {magnitude:.3e} (expected < 1e-10); sum = {sum:?}"
    );
  }

  pub fn build(self) -> LinkedCellContainer {
    self.local_container
  }
}
