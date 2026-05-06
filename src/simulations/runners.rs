use nalgebra::Vector3;
use crate::data::units::{TEMPERATURE_U, TIME_U};
use crate::sim_core::world::boundary_constraint::EdgeCondition;
use crate::sim_core::world::integration::{IntegrationAlgorithm, TemperatureInfo, TimeIterationDistance};
use crate::sim_core::world::WorldType;
use crate::utils::units::celcius_to_kelvin;
use crate::utils::logging;
use crate::utils::logging::get_save_path;
use crate::simulations::examples::{dense_particles, one_particle_edge, sphere_particles, symmetric_triangle_test, triangle, two_particles_edge};

const TIME_STEP: f64 = 1e-17 / TIME_U;
const TEMPERATURE_CELCIUS: f64 = 1600.0;

const TEMPERATURE_KELVIN: f64 = 1600.0;
const Q_EFFECTIVE_MASS: f64 = 1000.;

pub fn dense_runner() {
  let save_path = get_save_path("/media/washindeiru/EE366BA9366B718F/md/output/".to_string());
  logging::init_logging(save_path.clone());

  let desired_temperatures = vec![
    TemperatureInfo::new(1600., TimeIterationDistance::Iteration { value: 20000 }),
    TemperatureInfo::new(200., TimeIterationDistance::Iteration { value: 2000 }),
  ].into_iter().map(
    |i| TemperatureInfo{
      desired_temperature: i.desired_temperature / TEMPERATURE_U,
      ..i
    }).collect();

  let num_iterations = (8e4) as usize;
  let save = true;

  let integration_algorithm = IntegrationAlgorithm::NoseHooverVerlet {
    desired_temperature: desired_temperatures,
    q_effective_mass: Q_EFFECTIVE_MASS,
  };

  let world_type = WorldType::BoxedWorld { task_worker_multiplier: 4.0 };
  let edge_condition = EdgeCondition::Periodic;

  let world_size = Vector3::new(50., 50., 50.);
  let offset = Vector3::new(1.7, 1.7, 1.7);

  dense_particles(TIME_STEP, save, save_path, num_iterations, 7., world_size, offset, integration_algorithm, world_type, edge_condition);
}

pub fn one_particle_edge_runner() {
  let save_path = get_save_path("../output/".to_string());
  logging::init_logging(save_path.clone());

  let desired_temperatures = vec![
    TemperatureInfo::new(600., TimeIterationDistance::Iteration { value: 1000 }),
  ].into_iter().map(
    |i| TemperatureInfo{
      desired_temperature: i.desired_temperature / TEMPERATURE_U,
      ..i
    }).collect();

  let num_iterations = (4e4) as usize;
  let save = true;

  let integration_algorithm = IntegrationAlgorithm::NoseHooverVerlet {
    desired_temperature: desired_temperatures,
    q_effective_mass: Q_EFFECTIVE_MASS,
  };

  let world_type = WorldType::BoxedWorld { task_worker_multiplier: 2.0 };
  let edge_condition = EdgeCondition::Periodic;

  one_particle_edge(TIME_STEP, save, save_path, num_iterations,
                  integration_algorithm, world_type, edge_condition);
}

pub fn two_particles_edge_runner() {
  let save_path = get_save_path("../output/".to_string());
  logging::init_logging(save_path.clone());

  let desired_temperatures = vec![
    TemperatureInfo::new(600., TimeIterationDistance::Iteration { value: 1000 }),
  ].into_iter().map(
    |i| TemperatureInfo{
      desired_temperature: i.desired_temperature / TEMPERATURE_U,
      ..i
    }).collect();

  let num_iterations = (4e4) as usize;
  let save = true;

  let integration_algorithm = IntegrationAlgorithm::NoseHooverVerlet {
    desired_temperature: desired_temperatures,
    q_effective_mass: Q_EFFECTIVE_MASS,
  };

  let world_type = WorldType::BoxedWorld { task_worker_multiplier: 2.0 };
  let edge_condition = EdgeCondition::Periodic;

  two_particles_edge(TIME_STEP, save, save_path, num_iterations,
                    integration_algorithm, world_type, edge_condition);
}

pub fn sphere_runner() {
  let save_path = get_save_path("/media/washindeiru/EE366BA9366B718F/md/output/".to_string());
  logging::init_logging(save_path.clone());

  let temperature_kelvin_unitless = celcius_to_kelvin(TEMPERATURE_CELCIUS) / TEMPERATURE_U;

  let num_iterations = (2e4) as usize;
  let save = true;

  let integration_algorithm = IntegrationAlgorithm::NoseHooverVerlet {
    desired_temperature: vec![TemperatureInfo::new(temperature_kelvin_unitless, TimeIterationDistance::Iteration { value: num_iterations })],
    q_effective_mass: Q_EFFECTIVE_MASS,
  };
  let world_type = WorldType::BoxedWorld { task_worker_multiplier: 2.0 };

  sphere_particles(TIME_STEP, save, save_path, num_iterations, 30, integration_algorithm, world_type);
}

pub fn triangle_runner() {
  let save_path = get_save_path("/media/washindeiru/EE366BA9366B718F/md/output/".to_string());
  logging::init_logging(save_path.clone());

  let temperature_kelvin_unitless = celcius_to_kelvin(TEMPERATURE_CELCIUS) / TEMPERATURE_U;

  let num_iterations = (6e4) as usize;
  let save = true;

  let integration_algorithm = IntegrationAlgorithm::NoseHooverVerlet {
    desired_temperature: vec![TemperatureInfo::new(temperature_kelvin_unitless, TimeIterationDistance::Iteration { value: num_iterations })],
    q_effective_mass: Q_EFFECTIVE_MASS,
  };
  let world_type = WorldType::BoxedWorld { task_worker_multiplier: 2.0 };

  triangle(TIME_STEP, save, save_path, num_iterations, integration_algorithm, world_type);
}

pub fn symmetric_triangle_test_runner() {
  let save_path = get_save_path("/media/washindeiru/EE366BA9366B718F/md/output/".to_string());
  logging::init_logging(save_path.clone());

  let temperature_kelvin_unitless = celcius_to_kelvin(TEMPERATURE_CELCIUS) / TEMPERATURE_U;

  let num_iterations = (1e3) as usize;
  let save = true;

  let integration_algorithm = IntegrationAlgorithm::NoseHooverVerlet {
    desired_temperature: vec![TemperatureInfo::new(temperature_kelvin_unitless, TimeIterationDistance::Iteration { value: num_iterations })],
    q_effective_mass: Q_EFFECTIVE_MASS,
  };

  // let integration_algorithm = IntegrationAlgorithm::VelocityVerlet;

  let world_type = WorldType::BoxedWorld {task_worker_multiplier: 4.0};

  symmetric_triangle_test(TIME_STEP, save, save_path, num_iterations, integration_algorithm, world_type);
}

// fn two_particles_runner() {
//   logging::init_logging("two_particles_runner");
//   ting simulation...");
//
//   let num_iterations = 1e4 as usize;
//   let save = true;
//   let use_thermostat = false;
//   let verbose = true;
//
//   two_particles(save, num_iterations, use_thermostat, verbose);
// }
