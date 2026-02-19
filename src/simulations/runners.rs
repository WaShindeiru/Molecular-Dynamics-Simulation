use log::info;
use nalgebra::Vector3;
use rand::distr::Uniform;
use carbon_nanotube::data::units::{TEMPERATURE_U, TIME_U};
use carbon_nanotube::sim_core::world::integration::{IntegrationAlgorithm, IntegrationAlgorithmParams};
use carbon_nanotube::sim_core::world::WorldType;
use carbon_nanotube::utils::units::celcius_to_kelvin;
use crate::simulations::examples::{dense_particles, sphere_particles, symmetric_triangle_test, triangle};

const TIME_STEP: f64 = 1e-18 / TIME_U;
const TEMPERATURE_CELCIUS: f64 = 1600.0;

const TEMPERATURE_KELVIN: f64 = 1600.0;
const Q_EFFECTIVE_MASS: f64 = 0.1;

pub fn dense_runner() {
  env_logger::init();
  info!("Starting simulation...");

  // let temperature_kelvin_unitless = celcius_to_kelvin(TEMPERATURE_CELCIUS) / TEMPERATURE_U;
  let temperature_kelvin_unitless = TEMPERATURE_KELVIN / TEMPERATURE_U;

  let num_iterations = (3.2e4) as usize;
  let save = true;

  let params = IntegrationAlgorithmParams::NoseHooverVerlet {
    desired_temperature: temperature_kelvin_unitless,
    q_effective_mass: Q_EFFECTIVE_MASS,
  };

  let integration_algorithm = IntegrationAlgorithm::NoseHooverVerlet;
  let world_type = WorldType::BoxedWorld;

  let world_size = Vector3::new(27., 27., 27.);
  let offset = Vector3::new(3., 3., 3.);

  dense_particles(TIME_STEP, save, num_iterations, 3.1, world_size, offset,
                  integration_algorithm, params, world_type);
}

pub fn sphere_runner() {
  env_logger::init();
  info!("Starting simulation...");

  let temperature_kelvin_unitless = celcius_to_kelvin(TEMPERATURE_CELCIUS) / TEMPERATURE_U;

  let num_iterations = (2e4) as usize;
  let save = true;

  let params = IntegrationAlgorithmParams::NoseHooverVerlet {
    desired_temperature: temperature_kelvin_unitless,
    q_effective_mass: Q_EFFECTIVE_MASS,
  };

  let integration_algorithm = IntegrationAlgorithm::NoseHooverVerlet;
  let world_type = WorldType::BoxedWorld;

  sphere_particles(TIME_STEP, save, num_iterations, 30, integration_algorithm, params, world_type);
}

pub fn triangle_runner() {
  env_logger::init();
  info!("Starting simulation...");

  let temperature_kelvin_unitless = celcius_to_kelvin(TEMPERATURE_CELCIUS) / TEMPERATURE_U;

  let num_iterations = (3e4) as usize;
  let save = true;

  let params = IntegrationAlgorithmParams::NoseHooverVerlet {
    desired_temperature: temperature_kelvin_unitless,
    q_effective_mass: Q_EFFECTIVE_MASS,
  };

  let integration_algorithm = IntegrationAlgorithm::NoseHooverVerlet;
  let world_type = WorldType::BoxedWorld;

  triangle(TIME_STEP, save, num_iterations, integration_algorithm, params, world_type);
}

pub fn symmetric_triangle_test_runner() {
  env_logger::init();
  info!("Starting simulation...");

  let temperature_kelvin_unitless = celcius_to_kelvin(TEMPERATURE_CELCIUS) / TEMPERATURE_U;

  let num_iterations = (1e3) as usize;
  let save = true;

  let params = IntegrationAlgorithmParams::NoseHooverVerlet {
    desired_temperature: temperature_kelvin_unitless,
    q_effective_mass: Q_EFFECTIVE_MASS,
  };

  let integration_algorithm = IntegrationAlgorithm::NoseHooverVerlet;

  // let params = IntegrationAlgorithmParams::VelocityVerlet;
  // let integration_algorithm = IntegrationAlgorithm::VelocityVerlet;

  let world_type = WorldType::BoxedWorld;

  symmetric_triangle_test(TIME_STEP, save, num_iterations, integration_algorithm, params, world_type);
}

// fn two_particles_runner() {
//   env_logger::init();
//   info!("Starting simulation...");
//
//   let num_iterations = 1e4 as usize;
//   let save = true;
//   let use_thermostat = false;
//   let verbose = true;
//
//   two_particles(save, num_iterations, use_thermostat, verbose);
// }