use log::info;
use nalgebra::Vector3;
use crate::data::types::AtomType;
use crate::data::units::{TEMPERATURE_U, TIME_U};
use crate::particle::{Particle, SafeAtomFactory};
use crate::sim_core::Engine;
use crate::sim_core::world::integration::IntegrationAlgorithm;
use crate::utils::units::celcius_to_kelvin;

mod data;
mod particle;
mod utils;
mod sim_core;
mod output;

const TIME_STEP: f64 = 1e-18 / TIME_U;
const TEMPERATURE_CELCIUS: f64 = 1000.0;
const Q_EFFECTIVE_MASS: f64 = 1000.0;


fn main() {
  sphere_runner();
}

// fn two_particles_runner() {
//   env_logger::init();
//   info!("Starting simulation...");
//
//   let num_iterations = 1e4 as usize;
//   let save = true;
//   let use_thermostat = false;
//   let verbose = true;
//
//   two_particles(save, num_iterations, use_thermostat, verbose);
// }

fn sphere_runner() {
  env_logger::init();
  info!("Starting simulation...");

  let temperature_kelvin_unitless = celcius_to_kelvin(TEMPERATURE_CELCIUS) / TEMPERATURE_U;

  let num_iterations = (2e4 + 4.0) as usize;
  let save = true;
  sphere_particles(save, num_iterations, 5, temperature_kelvin_unitless,
                              Q_EFFECTIVE_MASS);
}

// fn two_particles(save: bool, num_iterations: usize, use_thermostat: bool, verbose: bool) {
//   let simulation_size = Vector3::new(50., 50., 50.);
//   let atom_factory = SafeAtomFactory::new();
//
//   let custom_path: Vec<Vector3<f64>> = vec![Vector3::new(10., 10., 10.)];
//
//   let mut custom_path_2: Vec<Vector3<f64>> = Vec::new();
//
//   let step_size = 1e-2;
//   let num_of_steps = 1e3 as usize;
//
//   for i in 0..num_of_steps {
//     let position = 10. + step_size * i as f64;
//     let position_vector = Vector3::new(position, 10., 10.);
//     custom_path_2.push(position_vector);
//   }
//
//   let atom_0 = Particle::CustomPathAtom(atom_factory.get_atom_custom_path(AtomType::Fe, custom_path_2));
//   let static_atom = Particle::CustomPathAtom(atom_factory.get_atom_custom_path(AtomType::Fe, custom_path));
//
//   let atoms: Vec<Particle> = vec![atom_0, static_atom];
//
//   let mut engine = Engine::new_from_atoms(
//     atoms, simulation_size, TIME_STEP,
//     num_iterations);
//
//   engine.run(save, use_thermostat, TEMPERATURE_CELCIUS, Q_EFFECTIVE_MASS, verbose);
// }

fn sphere_particles(save: bool, num_iterations: usize, num_particles: usize, temperature_kelvin: f64,
                    q_effective_mass: f64) {
  let simulation_size = Vector3::new(16., 16., 16.);

  let atom_factory = SafeAtomFactory::new();

  let center = Vector3::new(8., 8., 8.); // Center of the simulation space
  let sphere_radius = 3.4; // Radius of the sphere (leaving some margin from boundaries)

  let mut atoms: Vec<Particle> = Vec::new();

  for _ in 0..num_particles {
    // Generate random point in unit sphere using rejection sampling
    let mut position;
    loop {
      // Generate random coordinates in [-1, 1] range
      let x = rand::random::<f64>() * 2.0 - 1.0;
      let y = rand::random::<f64>() * 2.0 - 1.0;
      let z = rand::random::<f64>() * 2.0 - 1.0;

      // Check if point is within unit sphere
      if x * x + y * y + z * z <= 1.0 {
        // Scale to sphere radius and translate to center
        position = Vector3::new(
          center.x + x * sphere_radius,
          center.y + y * sphere_radius,
          center.z + z * sphere_radius,
        );
        break;
      }
    }

    let atom_type = match rand::random::<f64>() {
      v if v < 0.5 => AtomType::Fe,
      _ => AtomType::C,
    };

    let atom = atom_factory.get_atom(atom_type, position, Vector3::new(0., 0., 0.));
    atoms.push(Particle::Atom(atom));
  }

  let max_iteration_till_reset = 1e4 as usize;
  let save_laamps = true;
  let save_verbose = true;
  let save_all_iterations = false;
  let one_frame_duration = 1e-16 / TIME_U;

  let mut engine = Engine::new_from_atoms(
    atoms, simulation_size, TIME_STEP,
    num_iterations,
    max_iteration_till_reset,
    save,
    save_laamps,
    save_verbose,
    save_all_iterations,
    one_frame_duration,
  );

  let params = IntegrationAlgorithm::NoseHooverVerlet {
    desired_temperature: temperature_kelvin,
    q_effective_mass: q_effective_mass,
  };

  engine.run(&params, TIME_STEP);
}

// fn small_box_one_particle(save: bool, num_iterations: usize) {
//   let simulation_size = Vector3::new(8., 8., 8.);
//
//   let atom_factory = SafeAtomFactory::new();
//
//   let atom_0 = Particle::Atom(atom_factory.get_atom(AtomType::Fe, Vector3::new(2., 2., 2.), Vector3::new(2., 0., 0.)));
//   let atoms: Vec<Particle> = vec![atom_0];
//
//   let mut engine = Engine::new_from_atoms(atoms, simulation_size, TIME_STEP, num_iterations);
//
//   engine.run(save);
// }
//
// fn small_box(save: bool, num_iterations: usize, num_particles: usize) {
//   let simulation_size = Vector3::new(8., 8., 8.);
//
//   let atom_factory = SafeAtomFactory::new();
//
//   let mut atoms: Vec<Particle> = Vec::new();
//
//   for _ in 0..num_particles {
//     let x = rand::random::<f64>() * 6.0 + 1.0;
//     let y = rand::random::<f64>() * 6.0 + 1.0;
//     let z = rand::random::<f64>() * 6.0 + 1.0;
//
//     let atom_type = match rand::random::<f64>() {
//       v if v < 0.5 => AtomType::Fe,
//       _ => AtomType::C,
//     };
//
//     let atom_0 = Particle::Atom(atom_factory.get_atom(atom_type, Vector3::new(x, y, z), Vector3::new(0., 0., 0.)));
//
//
//   }
//
//   let atom_0 = Particle::Atom(atom_factory.get_atom(AtomType::Fe, Vector3::new(10., 10., 10.), Vector3::new(0., 0., 0.)));
//   let atom_1 = Particle::Atom(atom_factory.get_atom(AtomType::C, Vector3::new(12., 10., 10.), Vector3::new(0., 0., 0.)));
//   let atoms: Vec<Particle> = vec![atom_0, atom_1];
//
//   let mut engine = Engine::new_from_atoms(atoms, simulation_size, TIME_STEP, num_iterations);
//
//   engine.run(save);
// }
//
// fn sphere_particles(save: bool, num_iterations: usize, num_particles: usize) {
//   let simulation_size = Vector3::new(16., 16., 16.);
//
//   let atom_factory = SafeAtomFactory::new();
//
//   let center = Vector3::new(8., 8., 8.); // Center of the simulation space
//   let sphere_radius = 3.4; // Radius of the sphere (leaving some margin from boundaries)
//
//   let mut atoms: Vec<Particle> = Vec::new();
//
//   for _ in 0..num_particles {
//     // Generate random point in unit sphere using rejection sampling
//     let mut position;
//     loop {
//       // Generate random coordinates in [-1, 1] range
//       let x = rand::random::<f64>() * 2.0 - 1.0;
//       let y = rand::random::<f64>() * 2.0 - 1.0;
//       let z = rand::random::<f64>() * 2.0 - 1.0;
//
//       // Check if point is within unit sphere
//       if x * x + y * y + z * z <= 1.0 {
//         // Scale to sphere radius and translate to center
//         position = Vector3::new(
//           center.x + x * sphere_radius,
//           center.y + y * sphere_radius,
//           center.z + z * sphere_radius,
//         );
//         break;
//       }
//     }
//
//     let atom_type = match rand::random::<f64>() {
//       v if v < 0.5 => AtomType::Fe,
//       _ => AtomType::C,
//     };
//
//     let atom = atom_factory.get_atom(atom_type, position, Vector3::new(0., 0., 0.));
//     atoms.push(Particle::Atom(atom));
//   }
//
//   let mut engine = Engine::new_from_atoms(atoms, simulation_size, TIME_STEP, num_iterations);
//
//   engine.run(save);
// }
//
//
// fn single_atom(save: bool, num_iterations: usize) {
//   let simulation_size = Vector3::new(50., 50., 50.);
//
//   let custom_path: Vec<Vector3<f64>> = vec![Vector3::new(10., 12.5, 10.)];
//
//   let atom_factory = SafeAtomFactory::new();
//   let atom_0 = Particle::Atom(atom_factory.get_atom(AtomType::Fe, Vector3::new(10., 10., 10.), Vector3::new(0., 0., 0.)));
//   let static_atom = Particle::CustomPathAtom(atom_factory.get_atom_custom_path(AtomType::Fe, custom_path));
//   let atoms: Vec<Particle> = vec![atom_0, static_atom];
//
//   let mut engine = Engine::new_from_atoms(atoms, simulation_size, TIME_STEP, num_iterations);
//
//   engine.run(save);
// }
//
// fn symmetric_triangle_test(save: bool, num_iterations: usize) {
//   let simulation_size = Vector3::new(50., 50., 50.);
//   let atom_factory = SafeAtomFactory::new();
//
//   let atom_0_path = vec![Vector3::new(10., 10., 10.)];
//   let atom_1_path = vec![Vector3::new(10., 12., 10.)];
//
//   let mut atom_2_path: Vec<Vector3<f64>> = Vec::new();
//   let start_x = 11.73205;
//   let end_x = 12.26795;
//   for i in 0..100 {
//     let t = i as f64 / 99.;
//     let x = start_x + t * (end_x - start_x);
//     let y = 11.;
//     let z = 10.;
//     atom_2_path.push(Vector3::new(x, y, z));
//   }
//
//   let atom_0 = Particle::CustomPathAtom(atom_factory.get_atom_custom_path(AtomType::Fe, atom_0_path));
//   let atom_1 = Particle::CustomPathAtom(atom_factory.get_atom_custom_path(AtomType::Fe, atom_1_path));
//   let atom_2 = Particle::CustomPathAtom(atom_factory.get_atom_custom_path(AtomType::Fe, atom_2_path));
//   let atoms: Vec<Particle> = vec![atom_0, atom_1, atom_2];
//
//   let mut engine = Engine::new_from_atoms(atoms, simulation_size, TIME_STEP, num_iterations);
//
//   engine.run(save);
// }
//
// fn triangle(save: bool, num_iterations: usize) {
//   let simulation_size = Vector3::new(50., 50., 50.);
//
//   let atom_factory = SafeAtomFactory::new();
//   let atom_0 = Particle::Atom(atom_factory.get_atom(AtomType::Fe, Vector3::new(10., 10., 10.), Vector3::new(0., 0., 0.)));
//   let atom_1 = Particle::Atom(atom_factory.get_atom(AtomType::Fe, Vector3::new(10., 12.8, 10.), Vector3::new(0., 0., 0.)));
//   let atom_2 = Particle::Atom(atom_factory.get_atom(AtomType::Fe, Vector3::new(12.4248, 11.4, 10.), Vector3::new(0., 0., 0.)));
//   // let atom_3 = atom_factory.get_atom(AtomType::Fe, Vector3::new(12.8, 10., 12.8), Vector3::new(0., 0., 0.));
//   let atoms: Vec<Particle> = vec![atom_0, atom_1, atom_2];
//
//   let mut engine = Engine::new_from_atoms(atoms, simulation_size, TIME_STEP, num_iterations);
//
//   engine.run(save);
// }
//
// fn many_particles(save: bool, num_iterations: usize, num_particles: usize) {
//   let simulation_size = Vector3::new(100., 100., 100.);
//
//   let atom_factory = SafeAtomFactory::new();
//
//   let number_of_atoms = num_particles;
//   let lower_bound = 46.;
//   let upper_bound = 52.;
//
//   let mut atoms: Vec<Particle> = Vec::new();
//
//   for i in 0..number_of_atoms {
//     let atom = atom_factory.get_atom_random(AtomType::Fe, lower_bound, upper_bound);
//     atoms.push(Particle::Atom(atom));
//   }
//
//   let mut engine = Engine::new_from_atoms(atoms, simulation_size, TIME_STEP, num_iterations);
//
//   engine.run(save);
// }
