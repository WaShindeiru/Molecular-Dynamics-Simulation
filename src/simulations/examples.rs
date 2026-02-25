use log::info;
use nalgebra::Vector3;
use rand::distr::Uniform;
use carbon_nanotube::data::types::AtomType;
use carbon_nanotube::data::units::TIME_U;
use carbon_nanotube::particle::{Particle, SafeAtomFactory};
use carbon_nanotube::sim_core::Engine;
use carbon_nanotube::sim_core::world::integration::{IntegrationAlgorithm, IntegrationAlgorithmParams};
use carbon_nanotube::sim_core::world::saver::SaveOptions;
use carbon_nanotube::sim_core::world::WorldType;

use rand::prelude::*;
use rand_distr::Normal;
use carbon_nanotube::data::types::AtomType::{Fe, C};

pub fn symmetric_triangle_test(time_step: f64, save: bool, num_iterations: usize, integration_algorithm: IntegrationAlgorithm,
                               integration_algorithm_params: IntegrationAlgorithmParams, world_type: WorldType) {
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

  let atom_0 = atom_factory.get_atom_custom_path(AtomType::Fe, atom_0_path);
  let atom_1 = atom_factory.get_atom_custom_path(AtomType::Fe, atom_1_path);
  let atom_2 = atom_factory.get_atom_custom_path(AtomType::Fe, atom_2_path);
  let atoms: Vec<Particle> = vec![atom_0, atom_1, atom_2];

  let max_iteration_till_reset = 1e4 as usize;
  let save_laamps = true;
  let save_verbose = true;
  let save_all_iterations = true;
  let one_frame_duration = 1e-16 / TIME_U;

  let save_options = SaveOptions {
    save,
    save_laamps,
    save_verbose,
    save_path: "".to_string(),
  };

  let mut engine = Engine::new_from_atoms(
    atoms, simulation_size, time_step,
    num_iterations,
    max_iteration_till_reset,
    save_all_iterations,
    one_frame_duration,
    save_options,
    integration_algorithm,
    world_type
  );

  engine.run(&integration_algorithm_params, time_step);
}

pub fn triangle(time_step: f64, save: bool, num_iterations: usize, integration_algorithm: IntegrationAlgorithm,
            integration_algorithm_params: IntegrationAlgorithmParams, world_type: WorldType) {
  let simulation_size = Vector3::new(50., 50., 50.);

  let atom_factory = SafeAtomFactory::new();
  let atom_0 = atom_factory.get_atom(AtomType::Fe, Vector3::new(10., 10., 10.), Vector3::new(0., 0., 0.));
  let atom_1 = atom_factory.get_atom(AtomType::Fe, Vector3::new(10., 12.8, 10.), Vector3::new(0., 0., 0.));
  let atom_2 = atom_factory.get_atom(AtomType::Fe, Vector3::new(12.4248, 11.4, 10.), Vector3::new(0., 0., 0.));
  // let atom_3 = atom_factory.get_atom(AtomType::Fe, Vector3::new(12.8, 10., 12.8), Vector3::new(0., 0., 0.));
  let atoms: Vec<Particle> = vec![atom_0, atom_1, atom_2];

  let max_iteration_till_reset = 1e4 as usize;
  let save_laamps = true;
  let save_verbose = true;
  let save_all_iterations = true;
  let one_frame_duration = 1e-16 / TIME_U;

  let save_options = SaveOptions {
    save,
    save_laamps,
    save_verbose,
    save_path: "".to_string(),
  };

  let mut engine = Engine::new_from_atoms(
    atoms, simulation_size, time_step,
    num_iterations,
    max_iteration_till_reset,
    save_all_iterations,
    one_frame_duration,
    save_options,
    integration_algorithm,
    world_type
  );

  engine.run(&integration_algorithm_params, time_step);
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

pub fn sphere_particles(time_step: f64, save: bool, num_iterations: usize, num_particles: usize, integration_algorithm: IntegrationAlgorithm,
                    integration_algorithm_params: IntegrationAlgorithmParams, world_type: WorldType) {
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
    atoms.push(atom);
  }

  let max_iteration_till_reset = 1e4 as usize;
  let save_laamps = true;
  let save_verbose = true;
  let save_all_iterations = false;
  let one_frame_duration = 1e-16 / TIME_U;

  let save_options = SaveOptions {
    save,
    save_laamps,
    save_verbose,
    save_path: "".to_string(),
  };

  let mut engine = Engine::new_from_atoms(
    atoms, simulation_size, time_step,
    num_iterations,
    max_iteration_till_reset,
    save_all_iterations,
    one_frame_duration,
    save_options,
    integration_algorithm,
    world_type
  );

  engine.run(&integration_algorithm_params, time_step);
}

pub fn dense_particles(time_step: f64, save: bool, num_iterations: usize, particle_distance: f64, world_size: Vector3<f64>,
                       offset: Vector3<f64>, integration_algorithm: IntegrationAlgorithm,
                       integration_algorithm_params: IntegrationAlgorithmParams, world_type: WorldType) {
  let atom_factory = SafeAtomFactory::new();

  let mut atoms: Vec<Particle> = Vec::new();
  let mut count = 0;
  let mut fe_count = 0;
  let mut c_count = 0;

  let normal = Normal::new(0.0, 1.5).unwrap();
  let mut rng = rand::rng();

  let range = Uniform::new(0., 1.).unwrap();
  let mut rng_2 = rand::rng();

  let velocity_normal = Normal::new(0.0, 1e-2).unwrap();
  let mut rng_3 = rand::rng();

  let mut z = offset.z;
  while z + offset.z <= world_size.z {
    let mut y = offset.y;

    while y + offset.y <= world_size.y {
      let mut x = offset.x;

      while x + offset.x <= world_size.x {
        let x_ = loop {
          let x_candidate = x + normal.sample(&mut rng);
          if x_candidate >= 0.0 && x_candidate <= world_size.x {
            break x_candidate;
          }
        };

        let y_ = loop {
          let y_candidate = y + normal.sample(&mut rng);
          if y_candidate >= 0.0 && y_candidate <= world_size.y {
            break y_candidate;
          }
        };

        let z_ = loop {
          let z_candidate = z + normal.sample(&mut rng);
          if z_candidate >= 0.0 && z_candidate <= world_size.z {
            break z_candidate;
          }
        };

        let position = Vector3::new(x_, y_, z_);

        let atom_type = match range.sample(&mut rng_2) {
          x if x < 0.75 => {
            fe_count += 1;
            Fe
          },
          x if x >= 0.75 => {
            c_count += 1;
            C
          },
          _ => panic!("what's that?"),
        };

        let v_x = velocity_normal.sample(&mut rng_3);
        let v_y = velocity_normal.sample(&mut rng_3);
        let v_z = velocity_normal.sample(&mut rng_3);
        let particle = atom_factory.get_atom(atom_type, position, Vector3::new(v_x, v_y, v_z));
        atoms.push(particle);
        count += 1;

        x = x + particle_distance;
      }

      y = y + particle_distance;
    }

    z = z + particle_distance;
  }

  let max_iteration_till_reset = 2e3 as usize;
  let save_laamps = true;
  let save_verbose = false;
  let save_all_iterations = false;
  let one_frame_duration = 1e-16 / TIME_U;

  let save_options = SaveOptions {
    save,
    save_laamps,
    save_verbose,
    save_path: "".to_string(),
  };

  let mut engine = Engine::new_from_atoms(
    atoms, world_size, time_step,
    num_iterations,
    max_iteration_till_reset,
    save_all_iterations,
    one_frame_duration,
    save_options,
    integration_algorithm,
    world_type
  );

  info!("Starting simulation with {count} particles.");
  let fe_frac = fe_count as f64 / count as f64;
  let c_frac = c_count as f64 / count as f64;
  info!("Fe particles, count: {fe_count}, frac: {fe_frac}");
  info!("C particles, count: {c_count}, frac: {c_frac}");

  engine.run(&integration_algorithm_params, time_step);
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
//
//
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
