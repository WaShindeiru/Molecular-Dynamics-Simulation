use crate::sim_core::world::WorldType;
use crate::sim_core::world::boundary_constraint::EdgeCondition;
use crate::sim_core::world::boxed_world::box_task::task_manager::TaskManagerConfig;

#[derive(Debug, Clone, Copy, serde::Serialize, serde::Deserialize)]
#[serde(tag = "type")]
pub enum EdgeConditionFile {
  Simple {
    split: bool,
    trigger_small_subtask_size: usize,
  },
  Periodic {
    split: bool,
    trigger_small_subtask_size: usize,
  },
  PeriodicAll,
}

impl EdgeConditionFile {
  pub fn from_runtime(value: EdgeCondition) -> Self {
    match value {
      EdgeCondition::Simple { trigger_small_subtask_size, split } => EdgeConditionFile::Simple {
        trigger_small_subtask_size,
        split,
      },
      EdgeCondition::Periodic { trigger_small_subtask_size, split } => EdgeConditionFile::Periodic {
        trigger_small_subtask_size,
        split,
      },
      EdgeCondition::PeriodicAll => EdgeConditionFile::PeriodicAll,
    }
  }

  pub fn to_runtime(&self) -> EdgeCondition {
    match self {
      EdgeConditionFile::Simple { trigger_small_subtask_size, split } => EdgeCondition::Simple {
        trigger_small_subtask_size: *trigger_small_subtask_size,
        split: *split,
      },
      EdgeConditionFile::Periodic { trigger_small_subtask_size, split } => EdgeCondition::Periodic {
        trigger_small_subtask_size: *trigger_small_subtask_size,
        split: *split,
      },
      EdgeConditionFile::PeriodicAll => EdgeCondition::PeriodicAll,
    }
  }
}

#[derive(Debug, Clone, Copy, serde::Serialize, serde::Deserialize)]
#[serde(tag = "type")]
pub enum WorldTypeFile {
  SimpleWorld,
  BoxedWorld {
    task_manager_config: TaskManagerConfigFile,
  },
  LinkedCellWorld {
    task_manager_config: TaskManagerConfigFile,
  },
}

#[derive(Debug, Clone, Copy, serde::Serialize, serde::Deserialize)]
pub struct TaskManagerConfigFile {
  pub debug: bool,
  pub task_worker_multiplier: f64,
}

impl TaskManagerConfigFile {
  pub fn from_runtime(value: TaskManagerConfig) -> Self {
    Self {
      debug: value.debug,
      task_worker_multiplier: value.task_worker_multiplier,
    }
  }

  pub fn to_runtime(&self) -> TaskManagerConfig {
    TaskManagerConfig {
      debug: self.debug,
      task_worker_multiplier: self.task_worker_multiplier,
    }
  }
}

impl WorldTypeFile {
  pub fn from_runtime(value: WorldType) -> Self {
    match value {
      WorldType::SimpleWorld => WorldTypeFile::SimpleWorld,
      WorldType::BoxedWorld { task_manager_config } => WorldTypeFile::BoxedWorld {
        task_manager_config: TaskManagerConfigFile::from_runtime(task_manager_config),
      },
      WorldType::LinkedCellWorld { task_manager_config } => WorldTypeFile::LinkedCellWorld {
        task_manager_config: TaskManagerConfigFile::from_runtime(task_manager_config),
      },
    }
  }

  pub fn to_runtime(&self) -> WorldType {
    match self {
      WorldTypeFile::SimpleWorld => WorldType::SimpleWorld,
      WorldTypeFile::BoxedWorld { task_manager_config } => WorldType::BoxedWorld {
        task_manager_config: task_manager_config.to_runtime(),
      },
      WorldTypeFile::LinkedCellWorld { task_manager_config } => WorldType::LinkedCellWorld {
        task_manager_config: task_manager_config.to_runtime(),
      },
    }
  }
}
