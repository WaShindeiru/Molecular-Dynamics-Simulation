use nalgebra::Vector3;
use crate::data::types::AtomType;
use crate::data::units::TIME_U;
use crate::particle::{Atom, SafeAtomFactory};
use crate::sim_core::Engine;

mod data;
mod particle;
mod utils;
mod sim_core;
mod output;

const TIME_STEP: f64 = 10e-16 / TIME_U;

fn main() {
  test1();
  // println!("{}", (2.0_f64.sqrt() / 2.0));
  // println!("{}", (PI/4.0).sin());

  // let a = Vector3::new(-1, 0, 2);
  // let b = Vector3::new(1, 3, 4);
  //
  // println!("{}", a.dot(&b));
}

fn test1() {
  let simulation_size = Vector3::new(50., 50., 50.);

  let atom_factory = SafeAtomFactory::new();
  let atom_0 = atom_factory.get_atom(AtomType::Fe, Vector3::new(10., 10., 10.), Vector3::new(0., 0., 0.));
  let atom_1 = atom_factory.get_atom(AtomType::Fe, Vector3::new(10., 12.8, 10.), Vector3::new(0., 0., 0.));
  let atom_2 = atom_factory.get_atom(AtomType::Fe, Vector3::new(12.4248, 11.4, 10.), Vector3::new(0., 0., 0.));
  // let atom_3 = atom_factory.get_atom(AtomType::Fe, Vector3::new(12.8, 10., 12.8), Vector3::new(0., 0., 0.));
  let atoms: Vec<Atom> = vec![atom_0, atom_1, atom_2];

  let num_of_iterations = 10000;

  let mut engine = Engine::new_from_atoms(atoms, simulation_size, TIME_STEP, num_of_iterations);

  engine.run();
}

// fn test2() {
//   let simulation_size = Vector3::new(100., 100., 100.);
//
//   let atom_factory = SafeAtomFactory::new();
//
//   let number_of_atoms = 10;
//   let lower_bound = 46.;
//   let upper_bound = 52.;
//
//   let mut atoms: Vec<Atom> = Vec::new();
//
//   for i in 0..number_of_atoms {
//     let atom = atom_factory.get_atom_random(AtomType::Fe, lower_bound, upper_bound);
//     atoms.push(atom);
//   }
//
//   let num_of_iterations = 2000;
//
//   let mut engine = Engine::new_from_atoms(atoms, simulation_size, TIME_STEP, num_of_iterations);
//
//   engine.run();
// }

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

