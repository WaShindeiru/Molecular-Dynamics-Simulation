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
