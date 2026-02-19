use crate::simulations::runners::{dense_runner};

mod data;
mod particle;
mod utils;
mod sim_core;
mod output;
mod simulations;

fn main() {
  dense_runner()
}