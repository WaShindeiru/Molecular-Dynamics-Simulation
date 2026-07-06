use nalgebra::Vector3;
use std::collections::HashMap;
use std::sync::Arc;

use crate::data::SimulationConfig;
use crate::data::units::K_B;
use crate::sim_core::world::boundary_constraint::EdgeCondition;

use crate::particle::Particle;

use crate::sim_core::world::computation::compute_new_velocity;

use crate::sim_core::world::boxed_world::box_container::BoxContainer;
use crate::sim_core::world::boxed_world::box_container::sim_box::SimulationBox;
use crate::sim_core::world::boxed_world::box_task::ForceTaskParticleData;
use crate::sim_core::world::boxed_world::integration_cache::IntegrationCache;

pub struct ComputationCollector {
  config: SimulationConfig,
  integration_cache: Arc<IntegrationCache>,
  particles_modified: HashMap<usize, Particle>,
  force_visited: HashMap<usize, usize>,
}

impl ComputationCollector {
  pub fn from_integration_cache(config: SimulationConfig, cache: Arc<IntegrationCache>) -> Self {
    let mut particles_modified = cache.box_cache().all_particles_cloned();
    for particle in particles_modified.values_mut() {
      particle.set_iteration(particle.get_iteration() + 1);
    }

    ComputationCollector {
      config,
      integration_cache: cache,
      particles_modified,
      force_visited: HashMap::new(),
    }
  }

  pub fn apply_force_results<'a>(
    &mut self,
    force_results: impl IntoIterator<Item = (&'a usize, &'a ForceTaskParticleData)>,
  ) {

    for (particle_id, force_data) in force_results {
      if let Some(particle) = self.particles_modified.get_mut(particle_id) {

        let acceleration = force_data.force / particle.get_mass();
        particle.set_force(particle.get_force() + force_data.force);
        particle
          .set_potential_energy(particle.get_potential_energy() + force_data.potential_energy);
        particle.set_acceleration(particle.get_acceleration() + acceleration);
        self
          .force_visited
          .entry(*particle_id)
          .and_modify(|count| *count += 1)
          .or_insert(1);
      } else {
        panic!("sth wrong")
      }
    }
  }

  pub fn apply_gravity(&mut self, iteration: usize) {
    let potential_gravity_max = self
      .config
      .gravity_manager
      .lock()
      .expect("gravity manager lock poisoned")
      .get_gravity(iteration);
    let z_max = self.integration_cache.box_cache().config().world_size.z;

    for particle in self.particles_modified.values_mut() {
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
    let half_velocity_cache = self.integration_cache.half_velocity_cache();
    let particle_compliance = self.integration_cache.particle_compliance();
    let box_cache = self.integration_cache.box_cache();

    for (id, particle) in self.particles_modified.iter_mut() {
      if particle.is_custom_velocity_atom() {
        continue;
      }

      let half_velocity = half_velocity_cache.get(id).unwrap();
      let compliance = particle_compliance.get(id).unwrap();

      // Non-compliant particles that went through PartialVelocityStep use sub-steps of
      // time_step/subtask_size; their half_velocity must be combined with the same step here.
      let effective_time_step = if !compliance.compliant && collision_split && subtask_size > 1 {
        time_step / subtask_size as f64
      } else {
        time_step
      };

      let new_velocity = compute_new_velocity(*half_velocity, *particle.get_acceleration(), thermostat_epsilon, effective_time_step);

      #[cfg(debug_assertions)]
      log::debug!(
        "Updated particle {} in box {}: velocity={:?}, force={:?}, position={:?}",
        id,
        box_cache.particle_box_id(*id),
        new_velocity,
        particle.get_force(),
        particle.get_position()
      );

      particle.set_velocity(new_velocity);
    }
  }

  /// Sets the velocity on all `CustomVelocityAtom` particles to the prescribed value for
  /// the *next* iteration.  Workers will read this velocity directly from history to
  /// advance the particle's position in the next step.
  pub fn apply_custom_velocities(&mut self, custom_velocities: &HashMap<usize, Vector3<f64>>) {
    for (id, particle) in self.particles_modified.iter_mut() {
      if let Particle::CustomVelocityAtom(p) = particle {
        if let Some(&vel) = custom_velocities.get(id) {
          p.set_velocity(vel);
        }
      }
    }
  }

  pub fn get_mean_temperature(&self) -> f64 {
    let relevant: Vec<_> = self
      .particles_modified
      .values()
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
    self.particles_modified.values()
  }

  #[cfg(debug_assertions)]
  pub fn assert_zero_net_force(&self) {
    let sum = self.particles_modified.values().fold(
      nalgebra::Vector3::<f64>::zeros(),
      |acc, p| acc + p.get_force(),
    );
    let magnitude = sum.magnitude();
    debug_assert!(
      magnitude < 1e-10,
      "Newton's 3rd law violated: net force magnitude = {magnitude:.3e} (expected < 1e-10); sum = {sum:?}"
    );
  }

  pub fn into_box_container(self) -> BoxContainer<Arc<SimulationBox>> {
    let box_container_config = *self.integration_cache.box_cache().config();
    let particle_box_mapping = self.integration_cache.box_cache().box_id_cache().clone();
    BoxContainer::from_particles(
      box_container_config,
      &self.particles_modified,
      &particle_box_mapping,
    )
  }
}
