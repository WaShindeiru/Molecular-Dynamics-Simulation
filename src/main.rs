use log::info;
use nalgebra::Vector3;
use crate::data::types::AtomType;
use crate::data::units::TIME_U;
use crate::particle::{ParticleOperations, SafeAtomFactory};
use crate::sim_core::Engine;

mod data;
mod particle;
mod utils;
mod sim_core;
mod output;

const TIME_STEP: f64 = 1e-17 / TIME_U;

fn main() {
  env_logger::init();
  info!("Starting simulation...");

  let num_iterations = 10000;
  let save = true;
  many_particles(save, num_iterations);
}

fn single_atom(save: bool, num_iterations: usize) {
  let simulation_size = Vector3::new(50., 50., 50.);

  let custom_path: Vec<Vector3<f64>> = vec![Vector3::new(10., 12.5, 10.)];
  
  let atom_factory = SafeAtomFactory::new();
  let atom_0 = Box::new(atom_factory.get_atom(AtomType::Fe, Vector3::new(10., 10., 10.), Vector3::new(0., 0., 0.)));
  let static_atom = Box::new(atom_factory.get_atom_custom_path(AtomType::Fe, custom_path));
  let atoms: Vec<Box<dyn ParticleOperations>> = vec![atom_0, static_atom];

  let mut engine = Engine::new_from_atoms(atoms, simulation_size, TIME_STEP, num_iterations);

  engine.run(save);
}

fn symmetric_triangle_test(save: bool, num_iterations: usize) {
  let simulation_size = Vector3::new(50., 50., 50.);
  let atom_factory = SafeAtomFactory::new();

  let atom_0_path = vec![Vector3::new(10., 10., 10.)];
  let atom_1_path = vec![Vector3::new(10., 12., 10.)];

  let mut atom_2_path: Vec<Vector3<f64>> = Vec::new();
  let start_x = 11.73205;
  let end_x = 12.26795;
  for i in 0..100 {
    let t = i as f64 / 99.;
    let x = start_x + t * (end_x - start_x);
    let y = 11.;
    let z = 10.;
    atom_2_path.push(Vector3::new(x, y, z));
  }

  let atom_0 = Box::new(atom_factory.get_atom_custom_path(AtomType::Fe, atom_0_path));
  let atom_1 = Box::new(atom_factory.get_atom_custom_path(AtomType::Fe, atom_1_path));
  let atom_2 = Box::new(atom_factory.get_atom_custom_path(AtomType::Fe, atom_2_path));
  let atoms: Vec<Box<dyn ParticleOperations>> = vec![atom_0, atom_1, atom_2];

  let mut engine = Engine::new_from_atoms(atoms, simulation_size, TIME_STEP, num_iterations);

  engine.run(save);
}

fn triangle(save: bool, num_iterations: usize) {
  let simulation_size = Vector3::new(50., 50., 50.);

  let atom_factory = SafeAtomFactory::new();
  let atom_0 = Box::new(atom_factory.get_atom(AtomType::Fe, Vector3::new(10., 10., 10.), Vector3::new(0., 0., 0.)));
  let atom_1 = Box::new(atom_factory.get_atom(AtomType::Fe, Vector3::new(10., 12.8, 10.), Vector3::new(0., 0., 0.)));
  let atom_2 = Box::new(atom_factory.get_atom(AtomType::Fe, Vector3::new(12.4248, 11.4, 10.), Vector3::new(0., 0., 0.)));
  // let atom_3 = atom_factory.get_atom(AtomType::Fe, Vector3::new(12.8, 10., 12.8), Vector3::new(0., 0., 0.));
  let atoms: Vec<Box<dyn ParticleOperations>> = vec![atom_0, atom_1, atom_2];

  let mut engine = Engine::new_from_atoms(atoms, simulation_size, TIME_STEP, num_iterations);

  engine.run(save);
}

fn many_particles(save: bool, num_iterations: usize) {
  let simulation_size = Vector3::new(100., 100., 100.);

  let atom_factory = SafeAtomFactory::new();

  let number_of_atoms = 10;
  let lower_bound = 46.;
  let upper_bound = 52.;

  let mut atoms: Vec<Box<dyn ParticleOperations>> = Vec::new();

  for i in 0..number_of_atoms {
    let atom = atom_factory.get_atom_random(AtomType::Fe, lower_bound, upper_bound);
    atoms.push(Box::new(atom));
  }

  let mut engine = Engine::new_from_atoms(atoms, simulation_size, TIME_STEP, num_iterations);

  engine.run(save);
}

// use macroquad::prelude::*;
// 
// #[macroquad::main("BasicShapes")]
// async fn main() {
//     loop {
//         clear_background(RED);
// 
//         draw_line(40.0, 40.0, 100.0, 200.0, 15.0, BLUE);
//         draw_rectangle(screen_width() / 2.0 - 60.0, 100.0, 120.0, 60.0, GREEN);
//         draw_circle(screen_width() - 30.0, screen_height() - 30.0, 15.0, YELLOW);
//         draw_text("HELLO", 20.0, 20.0, 20.0, DARKGRAY);
// 
//         next_frame().await
//     }
// }

