use crate::simulations::runners::dense_runner;
use crate::simulations::various::{see_config_json, see_dense_generator_configuration, see_nanotube_generator_configuration};

mod cmd;
mod data;
mod particle;
mod persistence;
mod sim_core;
mod simulations;
mod utils;

fn main() {
  // dense_runner()

  if let Err(e) = cmd::run() {
    eprintln!("{e}");
    std::process::exit(1);
  }

  // see_config_json()

  // see_nanotube_generator_configuration()
}
