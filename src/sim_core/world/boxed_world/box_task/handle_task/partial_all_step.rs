use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use nalgebra::Vector3;

use crate::data::SimulationConfig;
use crate::particle::Particle;
use crate::sim_core::world::boundary_constraint::ParticleCompliance;
use crate::sim_core::world::boxed_world::box_container::BoxContainer;
use crate::sim_core::world::boxed_world::box_container::box_container_config::BoxContainerConfig;
use crate::sim_core::world::boxed_world::box_container::sim_box::SimulationBox;
use crate::sim_core::world::boxed_world::box_task::VelocityTaskParticleData;
use crate::sim_core::world::boxed_world::box_task::force_task_box_container::{
  ForceTaskBoxContainer, get_needed_box_id_periodic,
};
use crate::sim_core::world::boxed_world::box_task::handle_task::compute_forces::compute_forces_for_boxes_primary_only;
use crate::sim_core::world::boxed_world::box_task::handle_task::handle_partial_velocity_step::apply_velocity_constraint;
use crate::sim_core::world::boxed_world::integration::verlet_nose_hoover::computation::{
  HalfVelocityResult, verlet_noose_hoover_half_velocity_position,
};
use crate::sim_core::world::computation::compute_new_velocity;

pub struct PartialAllStepResult {
  pub particles: HashMap<usize, VelocityTaskParticleData>,
}

pub struct PartialAllStep {
  config: SimulationConfig,
  box_config: BoxContainerConfig,
  history: Arc<BoxContainer<Arc<SimulationBox>>>,

  current_iteration: usize,
  thermostat_epsilon: f64,
  dt_sub: f64,
  substep_count: usize,

  primary_box_ids: Vec<usize>,
  primary_ids: HashSet<usize>,

  working_box_ids: Vec<usize>,

  working: HashMap<usize, Particle>,
  half_velocity: HashMap<usize, Vector3<f64>>,
  compliance: HashMap<usize, ParticleCompliance>,

  current_substep: usize,
}

impl PartialAllStep {
  pub fn new(
    primary_box_ids: &[usize],
    history: Arc<BoxContainer<Arc<SimulationBox>>>,
    config: SimulationConfig,
    thermostat_epsilon: f64,
    current_iteration: usize,
  ) -> Self {
    let box_config = *history.config();
    let substep_count = config.correction.small_distance.substep_count.max(1);
    let dt_sub = config.time_step / substep_count as f64;

    let primary_ids = Self::collect_primary_ids(primary_box_ids, &history);
    let working = Self::initialize_working(&primary_ids, &history);

    let working_box_ids = primary_box_ids.to_vec();

    PartialAllStep {
      primary_box_ids: primary_box_ids.to_vec(),
      primary_ids,
      history,
      config,
      thermostat_epsilon,
      current_iteration,
      box_config,
      working,
      half_velocity: HashMap::new(),
      compliance: HashMap::new(),
      current_substep: 0,
      dt_sub,
      substep_count,
      working_box_ids: working_box_ids,
    }
  }

  fn collect_primary_ids(
    primary_box_ids: &[usize],
    history: &BoxContainer<Arc<SimulationBox>>,
  ) -> HashSet<usize> {
    let mut primary_ids = HashSet::new();
    for &box_id in primary_box_ids {
      let sim_box = history.get_box(box_id);
      for particle in sim_box.particles().values() {
        if !particle.is_custom_velocity_atom() {
          primary_ids.insert(particle.get_id());
        }
      }
    }
    primary_ids
  }

  fn initialize_working(
    primary_ids: &HashSet<usize>,
    history: &BoxContainer<Arc<SimulationBox>>,
  ) -> HashMap<usize, Particle> {
    primary_ids
      .iter()
      .map(|&id| {
        let mut particle = history.get_particle(id).as_ref().clone();
        particle.set_thermostat_work(0.0);
        (id, particle)
      })
      .collect()
  }

  fn apply_half_velocity(&mut self, half_result: &HalfVelocityResult) {
    for id in &self.primary_ids {
      self.half_velocity.insert(*id, half_result.half_velocity[id]);
      self.compliance.insert(*id, half_result.compliance[id]);

      let particle = self.working.get_mut(id).unwrap();
      particle.set_thermostat_work(
        particle.get_thermostat_work() + half_result.thermostat_work[id],
      );
      particle.update_position(half_result.new_position[id]);
    }
  }

  fn apply_edge_constraints(&mut self) {
    let constrained: HashMap<usize, VelocityTaskParticleData> = self
      .primary_ids
      .iter()
      .map(|id| {
        let particle = self.working.get(id).unwrap();
        (
          *id,
          VelocityTaskParticleData {
            half_velocity: *self.half_velocity.get(id).unwrap(),
            new_position: *particle.get_position(),
            thermostat_work: particle.get_thermostat_work(),
            compliance: *self.compliance.get(id).unwrap(),
          },
        )
      })
      .collect();

    let edge_adjusted = apply_velocity_constraint(constrained, self.config.edge_condition);
    for (id, data) in edge_adjusted {
      self.half_velocity.insert(id, data.half_velocity);
    }
  }

  fn build_force_container(&self) -> (ForceTaskBoxContainer, Vec<usize>) {
    let working_particles: Vec<Particle> = self.working.values().cloned().collect();
    let (working_boxes, working_box_ids) =
      BoxContainer::<SimulationBox>::new_local_with_particles_with_box_ids(
        self.box_config,
        &working_particles,
      );
    let working_box_ids: Vec<usize> = working_box_ids.iter().cloned().collect();

    let needed_ids = get_needed_box_id_periodic(&working_box_ids, &self.box_config);
    let static_ids: Vec<usize> = needed_ids.clone()
      .into_iter()
      .filter(|id| !self.primary_box_ids.contains(id))
      .collect();
    let static_view = self.history.view_select_boxes(&static_ids);

    let mut force_view = BoxContainer::merge_force_views(working_boxes, &static_view)
      .expect("working primaries and static halo must not share particle ids");

    // neighbour_atoms_periodic requires every periodic halo cell to be Some, even when empty.
    // merge only materializes boxes that contain particles; ensure the full halo exists.
    let mut boxes_to_ensure: HashSet<usize> = needed_ids.into_iter().collect();
    boxes_to_ensure.extend(self.primary_box_ids.iter().copied());
    let boxes_to_ensure: Vec<usize> = boxes_to_ensure.into_iter().collect();
    force_view.ensure_boxes_exist(&boxes_to_ensure);

    (
      ForceTaskBoxContainer::new(force_view, self.config.edge_condition),
      working_box_ids,
    )
  }

  fn apply_forces_and_velocity(&mut self, force_container: &ForceTaskBoxContainer) {
    let forces = compute_forces_for_boxes_primary_only(
      force_container,
      &self.working_box_ids,
      &self.primary_ids,
    );

    let z_max = self.config.world_size.z;
    let potential_gravity_max = self.config.potential_gravity_max;

    for id in &self.primary_ids {
      let particle = self.working.get_mut(id).unwrap();
      let force = *forces.get(id).unwrap_or(&Vector3::zeros());
      let gravity_force =
        -Vector3::new(0., 0., 1.) * potential_gravity_max * particle.get_mass() / z_max;
      let total_force = force + gravity_force;
      particle.set_force(total_force);
      particle.set_acceleration(total_force / particle.get_mass());

      let hv = *self.half_velocity.get(id).unwrap();
      let new_v = compute_new_velocity(
        hv,
        *particle.get_acceleration(),
        self.thermostat_epsilon,
        self.dt_sub,
      );
      particle.set_velocity(new_v);
    }
  }

  fn half_position_and_edge_constraints(&mut self) {
    let half_result = verlet_noose_hoover_half_velocity_position(
      self.working.values(),
      self.dt_sub,
      self.thermostat_epsilon,
      self.primary_ids.len(),
      self.current_iteration,
      &self.config.world_size,
      self.config.edge_condition,
    );

    self.apply_half_velocity(&half_result);
    self.apply_edge_constraints();
  }

  fn do_one_full_substep(&mut self) {
    self.half_position_and_edge_constraints();

    let (force_container, working_box_ids) = self.build_force_container();
    self.working_box_ids = working_box_ids;

    self.apply_forces_and_velocity(&force_container);
    self.current_substep += 1;
  }

  fn do_one_position_substep(&mut self) {
    self.half_position_and_edge_constraints();
    self.current_substep += 1;
  }

  fn build_result(&self) -> PartialAllStepResult {
    let particles: HashMap<usize, VelocityTaskParticleData> = self
      .primary_ids
      .iter()
      .map(|id| {
        let particle = self.working.get(id).unwrap();
        (
          *id,
          VelocityTaskParticleData {
            half_velocity: *self.half_velocity.get(id).unwrap(),
            new_position: *particle.get_position(),
            thermostat_work: particle.get_thermostat_work(),
            compliance: *self.compliance.get(id).unwrap(),
          },
        )
      })
      .collect();

    PartialAllStepResult { particles }
  }

  pub fn run(mut self) -> PartialAllStepResult {
    for _ in 0..self.substep_count.saturating_sub(1) {
      self.do_one_full_substep();
    }
    self.do_one_position_substep();

    debug_assert!(self.current_substep == self.substep_count);

    self.build_result()
  }
}

#[cfg(test)]
impl PartialAllStep {
  pub(crate) fn build_force_container_for_test(&self) -> ForceTaskBoxContainer {
    self.build_force_container().0
  }

  pub(crate) fn set_working_position_for_test(&mut self, id: usize, position: Vector3<f64>) {
    self.working.get_mut(&id).unwrap().update_position(position);
  }
}

#[cfg(test)]
mod build_force_container_test;
