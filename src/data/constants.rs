use std::collections::HashMap;
use std::sync::OnceLock;

use super::types::{Constant, InteractionType};

pub const ATOMIC_MASS_FE: f64 = 55.845;
pub const ATOMIC_MASS_C: f64 = 12.;

fn fefe_constants() -> &'static HashMap<Constant, f64> {
  static MAP_FEFE: OnceLock<HashMap<Constant, f64>> = OnceLock::new();
  MAP_FEFE.get_or_init(|| {
    let mut map : HashMap<Constant, f64> = HashMap::new();
    map.insert(Constant::D0, 1.5);
    map.insert(Constant::r0, 2.29);
    map.insert(Constant::Beta, 1.4);
    map.insert(Constant::S, 2.0693109);
    map.insert(Constant::Gamma, 0.0115751);
    map.insert(Constant::c, 1.2898716);
    map.insert(Constant::d, 0.3413219);
    map.insert(Constant::h, -0.26);
    map.insert(Constant::R, 3.15);
    map.insert(Constant::D, 0.2);
    map.insert(Constant::rf, 0.95);
    map.insert(Constant::bf, 2.9);

    map
  })
}

fn cc_constants() -> &'static HashMap<Constant, f64> {
  static MAP_CC: OnceLock<HashMap<Constant, f64>> = OnceLock::new();
  MAP_CC.get_or_init(|| {
    let mut map : HashMap<Constant, f64> = HashMap::new();
    map.insert(Constant::D0, 6.);
    map.insert(Constant::r0, 1.39);
    map.insert(Constant::Beta, 2.1);
    map.insert(Constant::S, 1.22);
    map.insert(Constant::Gamma, 2.0813e-4);
    map.insert(Constant::c, 330.);
    map.insert(Constant::d, 3.5);
    map.insert(Constant::h, 1.);
    map.insert(Constant::R, 1.85);
    map.insert(Constant::D, 0.15);
    map.insert(Constant::rf, 0.6);
    map.insert(Constant::bf, 8.);

    map
  })
}

fn fec_constants() -> &'static HashMap<Constant, f64> {
  static MAP_FEC: OnceLock<HashMap<Constant, f64>> = OnceLock::new();
  MAP_FEC.get_or_init(|| {
    let mut map : HashMap<Constant, f64> = HashMap::new();
    map.insert(Constant::D0, 4.82645134);
    map.insert(Constant::r0, 1.47736510);
    map.insert(Constant::Beta, 1.63208170);
    map.insert(Constant::S, 1.43134755);
    map.insert(Constant::Gamma, 0.00205862);
    map.insert(Constant::c, 8.95583221);
    map.insert(Constant::d, 0.72062047);
    map.insert(Constant::h, 0.87099874);
    map.insert(Constant::R, 2.5);
    map.insert(Constant::D, 0.2);
    map.insert(Constant::rf, 1.);
    map.insert(Constant::bf, 10.);

    map
  })
}

fn constants() -> &'static HashMap<InteractionType, &'static HashMap<Constant, f64>> {
  static MAP_: OnceLock<HashMap<InteractionType, &HashMap<Constant, f64>>> = OnceLock::new();
  MAP_.get_or_init(|| {
    let mut map: HashMap<InteractionType, &HashMap<Constant, f64>> = HashMap::new();
    map.insert(InteractionType::FeFe, fefe_constants());
    map.insert(InteractionType::CC, cc_constants());
    map.insert(InteractionType::FeC, fec_constants());

    map
  })
}

pub fn get_constant(interaction: &InteractionType, constant: Constant) -> f64 {
  let map = constants();
  let inner_map = map.get(interaction).expect("Wrong interaction type somehow!");
  *inner_map.get(&constant).expect("Wrong constant type somehow!")
}

pub fn get_box_size(interaction: &InteractionType) -> f64 {
  match interaction {
    InteractionType::CC => 2.,
    InteractionType::FeFe => 3.35,
    InteractionType::FeC => 2.7,
  }
}