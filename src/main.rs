use crate::simulations::runners::{dense_runner, dense_runner_old, one_particle_edge_runner, triangle_runner, two_particles_edge_runner};
use crate::simulations::various::see_config_json;

mod data;
mod particle;
mod utils;
mod sim_core;
mod output;
mod simulations;

fn main() {
  dense_runner_old();
}