use crate::data::units::{TIME_U, ValueUnits};

#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(tag = "type")]
pub enum TimeIterationDistance {
  Time { value: f64 },
  Iteration { value: usize },
}

impl TimeIterationDistance {
  pub fn to_iteration(self, time_step: f64) -> usize {
    match self {
      TimeIterationDistance::Iteration { value } => value,
      TimeIterationDistance::Time { value } => {
        let ratio = value / time_step;
        if !ratio.is_finite() || ratio <= 0.0 {
          return 1;
        }
        let rounded = ratio.round();
        if !rounded.is_finite() || rounded <= 0.0 {
          return 1;
        }
        let capped = rounded.min(usize::MAX as f64);
        (capped as usize).max(1)
      }
    }
  }

  pub fn to_value_units(&self, source: ValueUnits, target: ValueUnits) -> Self {
    match self {
      TimeIterationDistance::Time { value } => TimeIterationDistance::Time {
        value: value * ValueUnits::scale_between(source, target, TIME_U),
      },
      TimeIterationDistance::Iteration { value } => {
        TimeIterationDistance::Iteration { value: *value }
      }
    }
  }
}

#[derive(Debug, PartialEq, Eq, Clone, Copy, serde::Serialize, serde::Deserialize)]
#[repr(u64)]
pub enum AtomType {
  C = 0,
  Fe = 1,
  #[allow(non_camel_case_types)]
  C_nanotube = 2,
}

#[derive(Debug, Eq, Hash, PartialEq, Clone, Copy)]
pub enum InteractionType {
  FeFe,
  CC,
  FeC,
}

pub fn get_interaction_type(atom1: &AtomType, atom2: &AtomType) -> InteractionType {
  match (atom1, atom2) {
    (AtomType::Fe, AtomType::Fe) => InteractionType::FeFe,
    (AtomType::C | AtomType::C_nanotube, AtomType::C | AtomType::C_nanotube) => InteractionType::CC,
    (AtomType::Fe, AtomType::C | AtomType::C_nanotube)
    | (AtomType::C | AtomType::C_nanotube, AtomType::Fe) => InteractionType::FeC,
  }
}

#[derive(Eq, Hash, PartialEq)]
pub enum Constant {
  D0,
  r0,
  Beta,
  S,
  Gamma,
  c,
  d,
  h,
  R,
  D,
  rf,
  bf,
}
