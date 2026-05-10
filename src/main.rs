use crate::simulations::various::see_dense_generator_configuration;
use crate::simulations::runners::dense_runner;

mod data;
mod particle;
mod persistence;
mod sim_core;
mod simulations;
mod utils;
mod cmd;

fn main() {
  // dense_runner()
  
  if let Err(e) = cmd::run() {
    eprintln!("{e}");
    std::process::exit(1);
  }

  // see_dense_generator_configuration()
}
