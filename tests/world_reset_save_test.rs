mod common;

use carbon_nanotube::particle::SafeAtomFactory;
use carbon_nanotube::sim_core::world::boxed_world::box_task::task_manager::TaskManagerConfig;
use carbon_nanotube::sim_core::world::cell::TaskSplitVariant;
use carbon_nanotube::sim_core::world::{WorldType, boundary_constraint::EdgeCondition};

use common::box_helpers::{
  test_no_reset_runner, test_reset_world_no_missing_iterations_runner,
  test_reset_world_with_thermostat_runner, test_save_files_completeness_runner,
  test_single_reset_runner,
};

static FACTORY_TEST_LOCK: std::sync::Mutex<()> = std::sync::Mutex::new(());

#[test]
#[ignore]
fn test_reset_world_no_missing_iterations_boxed_simple() {
  test_reset_world_no_missing_iterations_runner(
    WorldType::BoxedWorld {
      task_manager_config: TaskManagerConfig {
        debug: false,
        task_worker_multiplier: 2.0,
        split: TaskSplitVariant::Floor,
      },
    },
    EdgeCondition::Simple {
      trigger_small_subtask_size: 1,
      split: EdgeCondition::DEFAULT_SPLIT,
    },
  )
}

#[test]
fn test_reset_world_no_missing_iterations_boxed_periodic() {
  let _lock = FACTORY_TEST_LOCK.lock().unwrap();
  SafeAtomFactory::reset_for_testing();
  test_reset_world_no_missing_iterations_runner(
    WorldType::BoxedWorld {
      task_manager_config: TaskManagerConfig {
        debug: false,
        task_worker_multiplier: 2.0,
        split: TaskSplitVariant::Floor,
      },
    },
    EdgeCondition::Periodic {
      trigger_small_subtask_size: 1,
      split: EdgeCondition::DEFAULT_SPLIT,
    },
  )
}

#[test]
#[ignore]
fn test_reset_world_with_thermostat_boxed_simple() {
  test_reset_world_with_thermostat_runner(
    WorldType::BoxedWorld {
      task_manager_config: TaskManagerConfig {
        debug: false,
        task_worker_multiplier: 2.0,
        split: TaskSplitVariant::Floor,
      },
    },
    EdgeCondition::Simple {
      trigger_small_subtask_size: 1,
      split: EdgeCondition::DEFAULT_SPLIT,
    },
  );
}

#[test]
fn test_reset_world_with_thermostat_boxed_periodic() {
  let _lock = FACTORY_TEST_LOCK.lock().unwrap();
  SafeAtomFactory::reset_for_testing();
  test_reset_world_with_thermostat_runner(
    WorldType::BoxedWorld {
      task_manager_config: TaskManagerConfig {
        debug: false,
        task_worker_multiplier: 2.0,
        split: TaskSplitVariant::Floor,
      },
    },
    EdgeCondition::Periodic {
      trigger_small_subtask_size: 1,
      split: EdgeCondition::DEFAULT_SPLIT,
    },
  );
}

#[test]
#[ignore]
fn test_save_files_completeness_boxed_simple() {
  test_save_files_completeness_runner(
    WorldType::BoxedWorld {
      task_manager_config: TaskManagerConfig {
        debug: false,
        task_worker_multiplier: 2.0,
        split: TaskSplitVariant::Floor,
      },
    },
    EdgeCondition::Simple {
      trigger_small_subtask_size: 1,
      split: EdgeCondition::DEFAULT_SPLIT,
    },
  );
}

#[test]
fn test_save_files_completeness_boxed_periodic() {
  let _lock = FACTORY_TEST_LOCK.lock().unwrap();
  SafeAtomFactory::reset_for_testing();
  test_save_files_completeness_runner(
    WorldType::BoxedWorld {
      task_manager_config: TaskManagerConfig {
        debug: false,
        task_worker_multiplier: 2.0,
        split: TaskSplitVariant::Floor,
      },
    },
    EdgeCondition::Periodic {
      trigger_small_subtask_size: 1,
      split: EdgeCondition::DEFAULT_SPLIT,
    },
  );
}

#[test]
#[ignore]
fn test_single_reset_boxed_simple() {
  test_single_reset_runner(
    WorldType::BoxedWorld {
      task_manager_config: TaskManagerConfig {
        debug: false,
        task_worker_multiplier: 2.0,
        split: TaskSplitVariant::Floor,
      },
    },
    EdgeCondition::Simple {
      trigger_small_subtask_size: 1,
      split: EdgeCondition::DEFAULT_SPLIT,
    },
  );
}

#[test]
fn test_single_reset_boxed_periodic() {
  let _lock = FACTORY_TEST_LOCK.lock().unwrap();
  SafeAtomFactory::reset_for_testing();
  test_single_reset_runner(
    WorldType::BoxedWorld {
      task_manager_config: TaskManagerConfig {
        debug: false,
        task_worker_multiplier: 2.0,
        split: TaskSplitVariant::Floor,
      },
    },
    EdgeCondition::Periodic {
      trigger_small_subtask_size: 1,
      split: EdgeCondition::DEFAULT_SPLIT,
    },
  );
}

#[test]
#[ignore]
fn test_no_reset_boxed_simple() {
  test_no_reset_runner(
    WorldType::BoxedWorld {
      task_manager_config: TaskManagerConfig {
        debug: false,
        task_worker_multiplier: 2.0,
        split: TaskSplitVariant::Floor,
      },
    },
    EdgeCondition::Simple {
      trigger_small_subtask_size: 1,
      split: EdgeCondition::DEFAULT_SPLIT,
    },
  );
}

#[test]
fn test_no_reset_boxed_periodic() {
  let _lock = FACTORY_TEST_LOCK.lock().unwrap();
  SafeAtomFactory::reset_for_testing();
  test_no_reset_runner(
    WorldType::BoxedWorld {
      task_manager_config: TaskManagerConfig {
        debug: false,
        task_worker_multiplier: 2.0,
        split: TaskSplitVariant::Floor,
      },
    },
    EdgeCondition::Periodic {
      trigger_small_subtask_size: 1,
      split: EdgeCondition::DEFAULT_SPLIT,
    },
  );
}
