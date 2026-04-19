mod common;

use carbon_nanotube::sim_core::world::{WorldType, boundary_constraint::EdgeCondition};

use common::box_helpers::{
  test_reset_world_no_missing_iterations_runner, 
  test_reset_world_with_thermostat_runner, 
  test_save_files_completeness_runner, 
  test_single_reset_runner, 
  test_no_reset_runner
};

#[test]
fn test_reset_world_no_missing_iterations_boxed_simple() {
  test_reset_world_no_missing_iterations_runner(WorldType::BoxedWorld, EdgeCondition::Simple)
}

#[test]
fn test_reset_world_no_missing_iterations_boxed_periodic() {
  test_reset_world_no_missing_iterations_runner(WorldType::BoxedWorld, EdgeCondition::Periodic)
}

#[test]
fn test_reset_world_with_thermostat_boxed_simple() {
  test_reset_world_with_thermostat_runner(WorldType::BoxedWorld, EdgeCondition::Simple);
}

#[test]
fn test_reset_world_with_thermostat_boxed_periodic() {
  test_reset_world_with_thermostat_runner(WorldType::BoxedWorld, EdgeCondition::Periodic);
}

#[test]
fn test_save_files_completeness_boxed_simple() {
  test_save_files_completeness_runner(WorldType::BoxedWorld, EdgeCondition::Simple);
}

#[test]
fn test_save_files_completeness_boxed_periodic() {
  test_save_files_completeness_runner(WorldType::BoxedWorld, EdgeCondition::Periodic);
}

#[test]
fn test_single_reset_boxed_simple() {
  test_single_reset_runner(WorldType::BoxedWorld, EdgeCondition::Simple);
}

#[test]
fn test_single_reset_boxed_periodic() {
  test_single_reset_runner(WorldType::BoxedWorld, EdgeCondition::Periodic);
}

#[test]
fn test_no_reset_boxed_simple() {
  test_no_reset_runner(WorldType::BoxedWorld, EdgeCondition::Simple);
}

#[test]
fn test_no_reset_boxed_periodic() {
  test_no_reset_runner(WorldType::BoxedWorld, EdgeCondition::Periodic);
}

