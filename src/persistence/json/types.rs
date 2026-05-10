use crate::sim_core::world::WorldType;
use crate::sim_core::world::boundary_constraint::EdgeCondition;

#[derive(Debug, Clone, Copy, serde::Serialize, serde::Deserialize)]
pub enum EdgeConditionFile {
  Simple,
  Periodic,
}

impl EdgeConditionFile {
  pub fn from_runtime(value: EdgeCondition) -> Self {
    match value {
      EdgeCondition::Simple => EdgeConditionFile::Simple,
      EdgeCondition::Periodic => EdgeConditionFile::Periodic,
    }
  }

  pub fn to_runtime(&self) -> EdgeCondition {
    match self {
      EdgeConditionFile::Simple => EdgeCondition::Simple,
      EdgeConditionFile::Periodic => EdgeCondition::Periodic,
    }
  }
}

#[derive(Debug, Clone, Copy, serde::Serialize, serde::Deserialize)]
#[serde(tag = "type")]
pub enum WorldTypeFile {
  SimpleWorld,
  BoxedWorld { task_worker_multiplier: f64 },
}

impl WorldTypeFile {
  pub fn from_runtime(value: WorldType) -> Self {
    match value {
      WorldType::SimpleWorld => WorldTypeFile::SimpleWorld,
      WorldType::BoxedWorld {
        task_worker_multiplier,
      } => WorldTypeFile::BoxedWorld {
        task_worker_multiplier,
      },
    }
  }

  pub fn to_runtime(&self) -> WorldType {
    match self {
      WorldTypeFile::SimpleWorld => WorldType::SimpleWorld,
      WorldTypeFile::BoxedWorld {
        task_worker_multiplier,
      } => WorldType::BoxedWorld {
        task_worker_multiplier: *task_worker_multiplier,
      },
    }
  }
}
