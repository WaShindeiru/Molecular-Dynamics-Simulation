use super::types::{InteractionType};

pub const ATOMIC_MASS_FE: f64 = 55.845;
pub const ATOMIC_MASS_C: f64 = 12.;

pub struct InteractionConstants {
  pub D0: f64,
  pub r0: f64,
  pub Beta: f64,
  pub S: f64,
  pub Gamma: f64,
  pub c: f64,
  pub d: f64,
  pub h: f64,
  pub R: f64,
  pub D: f64,
  pub rf: f64,
  pub bf: f64,
}

static FEFE: InteractionConstants = InteractionConstants {
  D0: 1.5,
  r0: 2.29,
  Beta: 1.4,
  S: 2.0693109,
  Gamma: 0.0115751,
  c: 1.2898716,
  d: 0.3413219,
  h: -0.26,
  R: 3.15,
  D: 0.2,
  rf: 0.95,
  bf: 2.9,
};

static CC: InteractionConstants = InteractionConstants {
  D0: 6.,
  r0: 1.39,
  Beta: 2.1,
  S: 1.22,
  Gamma: 2.0813e-4,
  c: 330.,
  d: 3.5,
  h: 1.,
  R: 1.85,
  D: 0.15,
  rf: 0.6,
  bf: 8.,
};

static FEC: InteractionConstants = InteractionConstants {
  D0: 4.82645134,
  r0: 1.47736510,
  Beta: 1.63208170,
  S: 1.43134755,
  Gamma: 0.00205862,
  c: 8.95583221,
  d: 0.72062047,
  h: 0.87099874,
  R: 2.5,
  D: 0.2,
  rf: 1.,
  bf: 10.,
};

pub fn get_constants(interaction: &InteractionType) -> &'static InteractionConstants {
  match interaction {
    InteractionType::FeFe => &FEFE,
    InteractionType::CC => &CC,
    InteractionType::FeC => &FEC,
  }
}

pub fn get_box_size(interaction: &InteractionType) -> f64 {
  match interaction {
    InteractionType::CC => 2.,
    InteractionType::FeFe => 3.35,
    InteractionType::FeC => 2.7,
  }
}
