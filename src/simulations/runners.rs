use nalgebra::Vector3;
use carbon_nanotube::data::units::{TEMPERATURE_U, TIME_U};
use carbon_nanotube::sim_core::world::integration::{IntegrationAlgorithm, TemperatureInfo, TimeIterationDistance};
use carbon_nanotube::sim_core::world::WorldType;
use carbon_nanotube::utils::units::celcius_to_kelvin;
use carbon_nanotube::utils::logging;
use carbon_nanotube::utils::logging::get_save_path;
use crate::simulations::examples::{dense_particles, sphere_particles, symmetric_triangle_test, triangle};

const TIME_STEP: f64 = 1e-18 / TIME_U;
const TEMPERATURE_CELCIUS: f64 = 1600.0;

const TEMPERATURE_KELVIN: f64 = 1600.0;
const Q_EFFECTIVE_MASS: f64 = 0.1;

pub fn dense_runner() {
  let save_path = get_save_path();
  logging::init_logging(save_path.clone());

  let desired_temperatures = vec![
    TemperatureInfo{desired_temperature: 200., distance: TimeIterationDistance::Iteration(30000)},
    TemperatureInfo{desired_temperature: 2400., distance: TimeIterationDistance::Iteration(30000)},
    ].into_iter().map(
      |i| TemperatureInfo{
        desired_temperature: i.desired_temperature / TEMPERATURE_U,
        distance: i.distance}).collect();

  // let desired_temperatures = vec![
  //   TemperatureInfo{desired_temperature: 1000., distance: TimeIterationDistance::Iteration(500)}
  // ].into_iter().map(
  //   |i| TemperatureInfo{
  //     desired_temperature: i.desired_temperature / TEMPERATURE_U,
  //     distance: i.distance}).collect();

  let num_iterations = (8e4) as usize;
  let save = true;

  let integration_algorithm = IntegrationAlgorithm::NoseHooverVerlet {
    desired_temperature: desired_temperatures,
    q_effective_mass: Q_EFFECTIVE_MASS,
  };

  let world_type = WorldType::BoxedWorld;

  let world_size = Vector3::new(30., 30., 30.);
  let offset = Vector3::new(3.4, 3.4, 3.4);

  dense_particles(TIME_STEP, save, save_path, num_iterations, 3.1, world_size, offset,
                  integration_algorithm, world_type);
}

pub fn sphere_runner() {
  let save_path = get_save_path();
  logging::init_logging(save_path.clone());

  let temperature_kelvin_unitless = celcius_to_kelvin(TEMPERATURE_CELCIUS) / TEMPERATURE_U;

  let num_iterations = (2e4) as usize;
  let save = true;

  let integration_algorithm = IntegrationAlgorithm::NoseHooverVerlet {
    desired_temperature: vec![TemperatureInfo {
      desired_temperature: temperature_kelvin_unitless,
      distance: TimeIterationDistance::Iteration(num_iterations),
    }],
    q_effective_mass: Q_EFFECTIVE_MASS,
  };
  let world_type = WorldType::BoxedWorld;

  sphere_particles(TIME_STEP, save, save_path, num_iterations, 30, integration_algorithm, world_type);
}

pub fn triangle_runner() {
  let save_path = get_save_path();
  logging::init_logging(save_path.clone());

  let temperature_kelvin_unitless = celcius_to_kelvin(TEMPERATURE_CELCIUS) / TEMPERATURE_U;

  let num_iterations = (6e4) as usize;
  let save = true;

  let integration_algorithm = IntegrationAlgorithm::NoseHooverVerlet {
    desired_temperature: vec![TemperatureInfo {
      desired_temperature: temperature_kelvin_unitless,
      distance: TimeIterationDistance::Iteration(num_iterations),
    }],
    q_effective_mass: Q_EFFECTIVE_MASS,
  };
  let world_type = WorldType::BoxedWorld;

  triangle(TIME_STEP, save, save_path, num_iterations, integration_algorithm, world_type);
}

pub fn symmetric_triangle_test_runner() {
  let save_path = get_save_path();
  logging::init_logging(save_path.clone());

  let temperature_kelvin_unitless = celcius_to_kelvin(TEMPERATURE_CELCIUS) / TEMPERATURE_U;

  let num_iterations = (1e3) as usize;
  let save = true;

  let integration_algorithm = IntegrationAlgorithm::NoseHooverVerlet {
    desired_temperature: vec![TemperatureInfo {
      desired_temperature: temperature_kelvin_unitless,
      distance: TimeIterationDistance::Iteration(num_iterations),
    }],
    q_effective_mass: Q_EFFECTIVE_MASS,
  };

  // let integration_algorithm = IntegrationAlgorithm::VelocityVerlet;

  let world_type = WorldType::BoxedWorld;

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