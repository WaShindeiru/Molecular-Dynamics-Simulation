use crate::data::{ParticleConfig, SimulationConfig, config::ConfigAll};
use crate::particle::Particle;
use crate::sim_core::world::WorldType;
use crate::sim_core::world::boundary_constraint::EdgeCondition;
use crate::sim_core::world::boxed_world::box_task::task_manager::TaskManagerConfig;
use crate::sim_core::world::cell::TaskSplitVariant;
use crate::sim_core::world::thermostat::{
  IntegrationAlgorithm, TemperatureInfo, TimeIterationDistance,
};
use crate::sim_core::world::saver::{FrameSamplingConfig, SaveOptions};
use nalgebra::Vector3;

pub struct SimulationConfigBuilder {
  atoms: Option<Vec<Particle>>,
  world_size: Option<Vector3<f64>>,
  gravity_schedule: Option<Vec<(usize, f64)>>,
  time_step: Option<f64>,
  num_of_iterations: Option<usize>,
  max_iteration_till_reset: Option<usize>,
  save_all_iterations: Option<bool>,
  one_frame_duration: Option<f64>,
  save_all_iterations_laamps: Option<bool>,
  one_frame_duration_laamps: Option<f64>,
  save_all_iterations_energy: Option<bool>,
  one_frame_duration_energy: Option<f64>,
  save_options: Option<SaveOptions>,
  integration_algorithm: Option<IntegrationAlgorithm>,
  world_type: Option<WorldType>,
  edge_condition: Option<EdgeCondition>,
}

impl SimulationConfigBuilder {
  pub fn new() -> Self {
    SimulationConfigBuilder {
      atoms: None,
      world_size: None,
      gravity_schedule: None,
      time_step: None,
      num_of_iterations: None,
      max_iteration_till_reset: None,
      save_all_iterations: None,
      one_frame_duration: None,
      save_all_iterations_laamps: None,
      one_frame_duration_laamps: None,
      save_all_iterations_energy: None,
      one_frame_duration_energy: None,
      save_options: None,
      integration_algorithm: None,
      world_type: None,
      edge_condition: None,
    }
  }

  pub fn atoms(mut self, atoms: Vec<Particle>) -> Self {
    self.atoms = Some(atoms);
    self
  }

  pub fn world_size(mut self, world_size: Vector3<f64>) -> Self {
    self.world_size = Some(world_size);
    self
  }

  pub fn gravity_schedule(mut self, gravity_schedule: Vec<(usize, f64)>) -> Self {
    self.gravity_schedule = Some(gravity_schedule);
    self
  }

  pub fn constant_gravity(mut self, gravity: f64) -> Self {
    self.gravity_schedule = Some(vec![(0, gravity)]);
    self
  }

  pub fn time_step(mut self, time_step: f64) -> Self {
    self.time_step = Some(time_step);
    self
  }

  pub fn num_of_iterations(mut self, num_of_iterations: usize) -> Self {
    self.num_of_iterations = Some(num_of_iterations);
    self
  }

  pub fn max_iteration_till_reset(mut self, max_iteration_till_reset: usize) -> Self {
    self.max_iteration_till_reset = Some(max_iteration_till_reset);
    self
  }

  pub fn save_all_iterations(mut self, save_all_iterations: bool) -> Self {
    self.save_all_iterations = Some(save_all_iterations);
    self
  }

  pub fn one_frame_duration(mut self, one_frame_duration: f64) -> Self {
    self.one_frame_duration = Some(one_frame_duration);
    self
  }

  pub fn save_all_iterations_laamps(mut self, save_all_iterations_laamps: bool) -> Self {
    self.save_all_iterations_laamps = Some(save_all_iterations_laamps);
    self
  }

  pub fn one_frame_duration_laamps(mut self, one_frame_duration_laamps: f64) -> Self {
    self.one_frame_duration_laamps = Some(one_frame_duration_laamps);
    self
  }

  pub fn save_all_iterations_energy(mut self, save_all_iterations_energy: bool) -> Self {
    self.save_all_iterations_energy = Some(save_all_iterations_energy);
    self
  }

  pub fn one_frame_duration_energy(mut self, one_frame_duration_energy: f64) -> Self {
    self.one_frame_duration_energy = Some(one_frame_duration_energy);
    self
  }

  pub fn save_options(mut self, save_options: SaveOptions) -> Self {
    self.save_options = Some(save_options);
    self
  }

  pub fn integration_algorithm(mut self, integration_algorithm: IntegrationAlgorithm) -> Self {
    self.integration_algorithm = Some(integration_algorithm);
    self
  }

  pub fn world_type(mut self, world_type: WorldType) -> Self {
    self.world_type = Some(world_type);
    self
  }

  pub fn edge_condition(mut self, edge_condition: EdgeCondition) -> Self {
    self.edge_condition = Some(edge_condition);
    self
  }

  pub fn assert_all_set(self) -> Result<Self, String> {
    let mut missing_fields = Vec::new();

    if self.world_size.is_none() {
      missing_fields.push("world_size");
    }
    if self.time_step.is_none() {
      missing_fields.push("time_step");
    }
    if self.num_of_iterations.is_none() {
      missing_fields.push("num_of_iterations");
    }
    if self.max_iteration_till_reset.is_none() {
      missing_fields.push("max_iteration_till_reset");
    }
    if self.save_options.is_none() {
      missing_fields.push("save_options");
    }
    if self.integration_algorithm.is_none() {
      missing_fields.push("integration_algorithm");
    }
    if self.world_type.is_none() {
      missing_fields.push("world_type");
    }
    if self.edge_condition.is_none() {
      missing_fields.push("edge_condition");
    }

    if missing_fields.is_empty() {
      Ok(self)
    } else {
      Err(format!(
        "The following fields are not set: {}",
        missing_fields.join(", ")
      ))
    }
  }

  pub fn get_missing_fields(&self) -> Vec<&'static str> {
    let mut missing = Vec::new();

    if self.world_size.is_none() {
      missing.push("world_size");
    }
    if self.time_step.is_none() {
      missing.push("time_step");
    }
    if self.num_of_iterations.is_none() {
      missing.push("num_of_iterations");
    }
    if self.max_iteration_till_reset.is_none() {
      missing.push("max_iteration_till_reset");
    }
    if self.save_options.is_none() {
      missing.push("save_options");
    }
    if self.integration_algorithm.is_none() {
      missing.push("integration_algorithm");
    }
    if self.world_type.is_none() {
      missing.push("world_type");
    }
    if self.edge_condition.is_none() {
      missing.push("edge_condition");
    }

    missing
  }

  pub fn build(self) -> Result<SimulationConfig, String> {
    let world_size = self.world_size.ok_or("world_size must be specified")?;
    let time_step = self.time_step.unwrap_or(1e-17);
    let default_save_all_iterations = self.save_all_iterations.unwrap_or(false);
    let default_one_frame_duration = self.one_frame_duration.unwrap_or(1e-16);

    let save_options_provided = self.save_options.is_some();
    let mut save_options = self.save_options.unwrap_or_default();

    // Priority: explicit builder field > SaveOptions field (if SaveOptions was provided) > default
    let save_all_iterations_laamps = self
      .save_all_iterations_laamps
      .unwrap_or_else(|| {
        if save_options_provided { save_options.save_all_iterations_laamps } else { default_save_all_iterations }
      });
    let save_all_iterations_energy = self
      .save_all_iterations_energy
      .unwrap_or_else(|| {
        if save_options_provided { save_options.save_all_iterations_energy } else { default_save_all_iterations }
      });
    save_options.save_all_iterations_laamps = save_all_iterations_laamps;
    save_options.save_all_iterations_energy = save_all_iterations_energy;
    save_options.laamps_sampling = FrameSamplingConfig::from_duration(
      time_step,
      self
        .one_frame_duration_laamps
        .unwrap_or(default_one_frame_duration),
      save_all_iterations_laamps,
    );
    save_options.energy_sampling = FrameSamplingConfig::from_duration(
      time_step,
      self
        .one_frame_duration_energy
        .unwrap_or(default_one_frame_duration),
      save_all_iterations_energy,
    );

    let gravity_schedule = self
      .gravity_schedule
      .unwrap_or_else(|| vec![(0, 1.0)]);

    Ok(SimulationConfig::new(
      world_size,
      gravity_schedule,
      time_step,
      self.num_of_iterations.unwrap_or(1000),
      self.max_iteration_till_reset.unwrap_or(1000),
      save_options,
      self.integration_algorithm.unwrap_or(
        IntegrationAlgorithm::NoseHooverVerlet {
          q_effective_mass: 1.,
          desired_temperature: vec![TemperatureInfo::new(
            2000.,
            TimeIterationDistance::Iteration { value: 20000 },
          )],
        },
      ),
      self.world_type.unwrap_or(WorldType::BoxedWorld {
        task_manager_config: TaskManagerConfig {
          debug: false,
          task_worker_multiplier: 4.0,
          split: TaskSplitVariant::Floor,
        },
      }),
      self.edge_condition.unwrap_or(EdgeCondition::Periodic {
        trigger_small_subtask_size: 1,
        split: EdgeCondition::DEFAULT_SPLIT,
      }),
      false,
    ))
  }

  pub fn build_all(self) -> Result<ConfigAll, String> {
    let atoms = self.atoms.clone().ok_or("atoms must be specified")?;
    let simulation_config = self.build()?;
    let particle_config = ParticleConfig::new(atoms);

    Ok(ConfigAll {
      simulation_config,
      particle_config,
    })
  }
}
