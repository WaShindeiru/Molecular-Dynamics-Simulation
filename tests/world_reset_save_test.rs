use nalgebra::Vector3;
use std::collections::HashSet;
use std::fs;
use std::path::Path;
use carbon_nanotube::output::WorldDTO;
use carbon_nanotube::data::types::AtomType;
use carbon_nanotube::particle::{Particle, SafeAtomFactory};
use carbon_nanotube::sim_core::Engine;
use carbon_nanotube::sim_core::world::integration::{IntegrationAlgorithm, IntegrationAlgorithmParams};
use carbon_nanotube::sim_core::world::saver::SaveOptions;
use carbon_nanotube::sim_core::world::WorldType;

const TIME_STEP: f64 = 1e-3;

fn validate_dto_type(world_type: &WorldType, world_dto: &WorldDTO) -> bool {
  match (world_type, world_dto) {
    (WorldType::SimpleWorld, WorldDTO::SimpleWorldDTO(_)) => true,
    (WorldType::BoxedWorld, WorldDTO::BoxedWorldDTO(_)) => true,
    _ => false,
  }
}

/// Test 1: Check that to_transfer_struct() returns all iterations without missing or repeating any
fn test_reset_world_no_missing_iterations_runner(world_type: WorldType) {
  // Setup: Create a small simulation
  let simulation_size = Vector3::new(10., 10., 10.);
  let atom_factory = SafeAtomFactory::new();

  let mut atoms: Vec<Particle> = Vec::new();
  for i in 0..5 {
    let position = Vector3::new(2.0 + i as f64 * 0.5, 5.0, 5.0);
    let velocity = Vector3::new(0.1, 0.0, 0.0);
    let atom = atom_factory.get_atom(AtomType::C, position, velocity);
    atoms.push(atom);
  }

  // Run simulation without resets to test to_transfer_struct
  let num_iterations = 45;
  let max_iteration_till_reset = 100; // Large enough to avoid resets
  let one_frame_duration = 1e-2;

  let save_options = SaveOptions {
    save: false,
    save_path: String::new(),
    save_laamps: false,
    save_verbose: false,
  };

  let integration_algorithm = IntegrationAlgorithm::NoseHooverVerlet;

  let mut engine = Engine::new_from_atoms(
    atoms,
    simulation_size,
    TIME_STEP,
    num_iterations,
    max_iteration_till_reset,
    true, // save all iterations
    one_frame_duration,
    save_options,
    integration_algorithm.clone(),
    world_type.clone(),
  );

  // TODO: change it back to verlet
  let params = IntegrationAlgorithmParams::NoseHooverVerlet {
    desired_temperature: 100.,
    q_effective_mass: 100.,
  };
  engine.run(&params, TIME_STEP);

  // Get the transfer struct
  let engine_dto = engine.to_transfer_struct();
  let world_dto = &engine_dto.world;

  validate_dto_type(&world_type, &world_dto);

  let atoms = match world_dto {
    WorldDTO::SimpleWorldDTO(simple_dto) => &simple_dto.atoms,
    WorldDTO::BoxedWorldDTO(boxed_dto) => &boxed_dto.box_container.atoms,
  };

  // Collect all iterations from all atoms
  let mut all_iterations: Vec<usize> = Vec::new();

  for iteration_atoms in atoms {
    if !iteration_atoms.is_empty() {
      // All atoms in the same iteration should have the same iteration number
      let iteration_number = iteration_atoms[0].iteration;

      // Verify all atoms in this timestep have the same iteration
      for atom in iteration_atoms {
        assert_eq!(
          atom.iteration, iteration_number,
          "Atoms in the same timestep have different iteration numbers"
        );
      }

      all_iterations.push(iteration_number);
    }
  }

  // Sort iterations
  all_iterations.sort();

  println!("Total iterations collected: {}", all_iterations.len());
  println!("Expected iterations: {}", num_iterations);
  println!("Iterations: {:?}", all_iterations);

  // Check no iterations are missing
  for i in 0..all_iterations.len() - 1 {
    let current = all_iterations[i];
    let next = all_iterations[i + 1];

    assert_eq!(
      next, current + 1,
      "Missing iteration detected: expected {} after {}, but got {}",
      current + 1, current, next
    );
  }

  // Check no iterations are repeated
  let unique_iterations: HashSet<usize> = all_iterations.iter().copied().collect();
  assert_eq!(
    unique_iterations.len(), all_iterations.len(),
    "Duplicate iterations detected: {} unique vs {} total",
    unique_iterations.len(), all_iterations.len()
  );

  // Check we have the right number of iterations (initial state + num_iterations steps)
  let expected_total = num_iterations + 1;
  assert_eq!(
    all_iterations.len(), expected_total,
    "Expected {} iterations (1 initial + {} steps), but got {}",
    expected_total, num_iterations, all_iterations.len()
  );

  println!("✓ Test passed: No missing or repeated iterations");
}

#[test]
fn test_reset_world_no_missing_iterations_simple_world() {
  test_reset_world_no_missing_iterations_runner(WorldType::SimpleWorld)
}

#[test]
fn test_reset_world_no_missing_iterations_boxed_world() {
  test_reset_world_no_missing_iterations_runner(WorldType::BoxedWorld)
}

/// Test 2: Test with Nose-Hoover thermostat and resets, verifying saved files
fn test_reset_world_with_thermostat_runner(world_type: WorldType) {
  let world_type_str = match world_type {
    WorldType::SimpleWorld => "simple",
    WorldType::BoxedWorld => "boxed",
  };
  let test_dir = format!("../test/test_reset_world_with_thermostat_{}", world_type_str);

  // Clean up test directory if it exists
  if Path::new(&test_dir).exists() {
    fs::remove_dir_all(&test_dir).expect("Failed to remove test directory");
  }
  fs::create_dir_all(&test_dir).expect("Failed to create test directory");

  let simulation_size = Vector3::new(10., 10., 10.);
  let atom_factory = SafeAtomFactory::new();

  let mut atoms: Vec<Particle> = Vec::new();
  for i in 0..3 {
    let position = Vector3::new(2.0 + i as f64 * 0.5, 5.0, 5.0);
    let velocity = Vector3::new(0.1, 0.0, 0.0);
    let atom = atom_factory.get_atom(AtomType::Fe, position, velocity);
    atoms.push(atom);
  }

  let num_iterations = 100;
  let max_iteration_till_reset = 8;
  let one_frame_duration = 1e-2;

  let save_options = SaveOptions {
    save: true,
    save_path: test_dir.clone(),
    save_laamps: false,
    save_verbose: true,
  };

  let integration_algorithm = IntegrationAlgorithm::NoseHooverVerlet;

  let mut engine = Engine::new_from_atoms_with_path(
    atoms,
    simulation_size,
    TIME_STEP,
    num_iterations,
    max_iteration_till_reset,
    false, // save all iterations, no I don't think so
    one_frame_duration,
    save_options,
    integration_algorithm.clone(),
    world_type,
  );

  let params = IntegrationAlgorithmParams::NoseHooverVerlet {
    desired_temperature: 300.0,
    q_effective_mass: 1000.0,
  };

  engine.run(&params, TIME_STEP);

  // Verify saved files contain all iterations despite resets
  let energy_iterations = parse_csv_iterations(&format!("{}/energy.csv", test_dir), 0);
  verify_iterations_complete(&energy_iterations, num_iterations, "energy.csv");

  println!("✓ Test with thermostat passed: Saved files contain all iterations despite resets");
}

#[test]
fn test_reset_world_with_thermostat_simple_world() {
  test_reset_world_with_thermostat_runner(WorldType::SimpleWorld);
}

#[test]
fn test_reset_world_with_thermostat_boxed_world() {
  test_reset_world_with_thermostat_runner(WorldType::BoxedWorld);
}

/// Test 3: Save to files and verify all CSV files contain correct iterations
fn test_save_files_completeness_runner(world_type: WorldType) {
  let world_type_str = match world_type {
    WorldType::SimpleWorld => "simple",
    WorldType::BoxedWorld => "boxed",
  };
  let test_dir = format!("../test/test_save_files_completeness_{}", world_type_str);

  // Clean up test directory if it exists
  if Path::new(&test_dir).exists() {
    fs::remove_dir_all(&test_dir).expect("Failed to remove test directory");
  }
  fs::create_dir_all(&test_dir).expect("Failed to create test directory");

  let simulation_size = Vector3::new(10., 10., 10.);
  let atom_factory = SafeAtomFactory::new();

  let atoms_count = 4;
  let mut atoms: Vec<Particle> = Vec::new();
  for i in 0..4 {
    let position = Vector3::new(3.0 + i as f64 * 0.3, 5.0, 5.0);
    let velocity = Vector3::new(0.05, 0.05, 0.0);
    let atom = atom_factory.get_atom(AtomType::C, position, velocity);
    atoms.push(atom);
  }

  let num_iterations = 50;
  let max_iteration_till_reset = 12; // Will have multiple save batches
  let one_frame_duration = 1e-2;

  let save_options = SaveOptions {
    save: true,
    save_path: test_dir.clone(),
    save_laamps: true,
    save_verbose: true,
  };

  // TODO: change back to verlet
  let integration_algorithm = IntegrationAlgorithm::NoseHooverVerlet;

  // Create engine with custom save path
  let mut engine = Engine::new_from_atoms_with_path(
    atoms,
    simulation_size,
    TIME_STEP,
    num_iterations,
    max_iteration_till_reset,
    true, // save all iterations
    one_frame_duration,
    save_options,
    integration_algorithm.clone(),
    world_type,
  );

  let params = IntegrationAlgorithmParams::NoseHooverVerlet {
    desired_temperature: 100.,
    q_effective_mass: 100.,
  };
  engine.run(&params, TIME_STEP);

  // Verify files exist
  assert!(Path::new(&format!("{}/energy.csv", test_dir)).exists(), "energy.csv not found");
  assert!(Path::new(&format!("{}/forces.csv", test_dir)).exists(), "forces.csv not found");
  assert!(Path::new(&format!("{}/positions.csv", test_dir)).exists(), "positions.csv not found");
  assert!(Path::new(&format!("{}/potential_energies.csv", test_dir)).exists(), "potential_energies.csv not found");

  // Check LAMMPS dump files exist
  let num_expected_dumps = num_iterations; // Since save_all_iterations is true
  for i in 0..num_expected_dumps {
    let dump_file = format!("{}/output_{}.dump", test_dir, i);
    assert!(Path::new(&dump_file).exists(), "LAMMPS dump file {} not found", dump_file);
  }

  println!("✓ All expected files exist");

  let num_atoms = atoms_count;

  // Parse energy.csv and check iterations (doesn't have per-atom data)
  let energy_iterations = parse_csv_iterations(&format!("{}/energy.csv", test_dir), 0);
  verify_iterations_complete(&energy_iterations, num_iterations, "energy.csv");

  // Use new validation function for files with per-atom data
  // These files have format: iteration,atom_id,...
  verify_iterations_complete_with_atoms(
    &format!("{}/forces.csv", test_dir),
    num_iterations,
    num_atoms,
    "forces.csv"
  );

  verify_iterations_complete_with_atoms(
    &format!("{}/positions.csv", test_dir),
    num_iterations,
    num_atoms,
    "positions.csv"
  );

  verify_iterations_complete_with_atoms(
    &format!("{}/potential_energies.csv", test_dir),
    num_iterations,
    num_atoms,
    "potential_energies.csv"
  );

  println!("✓ All CSV files contain complete iterations without gaps or duplicates");
  println!("✓ All per-atom CSV files have complete atom data for each iteration");
}

#[test]
fn test_save_files_completeness_simple_world() {
  test_save_files_completeness_runner(WorldType::SimpleWorld);
}

#[test]
fn test_save_files_completeness_boxed_world() {
  test_save_files_completeness_runner(WorldType::BoxedWorld);
}

/// Helper function to parse iterations from CSV file
fn parse_csv_iterations(file_path: &str, iteration_column: usize) -> Vec<usize> {
  let content = fs::read_to_string(file_path)
    .expect(&format!("Failed to read {}", file_path));

  let mut iterations = Vec::new();

  for line in content.lines() {
    let parts: Vec<&str> = line.split(',').collect();
    if parts.len() > iteration_column {
      if let Ok(iter) = parts[iteration_column].trim().parse::<usize>() {
        iterations.push(iter);
      }
    }
  }

  iterations
}

/// Helper function to verify iterations are complete (no gaps, no duplicates)
fn verify_iterations_complete(iterations: &Vec<usize>, expected_count: usize, file_name: &str) {
  let mut sorted_iterations = iterations.clone();
  sorted_iterations.sort();

  println!("Checking {}: found {} unique iterations out of {} total entries",
           file_name, sorted_iterations.len(), iterations.len());

  // Check for gaps
  if !sorted_iterations.is_empty() {
    for i in 0..sorted_iterations.len() - 1 {
      if sorted_iterations[i + 1] != sorted_iterations[i] + 1 {
        panic!(
          "{}: Gap detected between iterations {} and {}",
          file_name, sorted_iterations[i], sorted_iterations[i + 1]
        );
      }
    }
  }

  // Check the range covers expected iterations (including initial state at iteration 0)
  let unique_iterations: HashSet<usize> = sorted_iterations.iter().copied().collect();
  let expected_total = expected_count + 1; // +1 for initial state
  assert_eq!(
    unique_iterations.len(), expected_total,
    "{}: Expected {} unique iterations (1 initial + {} steps), found {}",
    file_name, expected_total, expected_count, unique_iterations.len()
  );
}

/// Helper function to verify all atoms are present for each iteration
/// Validates that for each iteration, all atom IDs (0 to num_atoms-1) are present
fn verify_iterations_complete_with_atoms(file_path: &str, expected_iterations: usize, num_atoms: usize, file_name: &str) {
  let content = fs::read_to_string(file_path)
    .expect(&format!("Failed to read {}", file_path));

  // Parse CSV: each line should have format: iteration,atom_id,...
  // Store data as: iteration -> set of atom_ids
  let mut iteration_atoms: std::collections::HashMap<usize, HashSet<usize>> = std::collections::HashMap::new();

  for line in content.lines() {
    let parts: Vec<&str> = line.split(',').collect();
    if parts.len() >= 2 {
      if let (Ok(iteration), Ok(atom_id)) = (parts[0].trim().parse::<usize>(), parts[1].trim().parse::<usize>()) {
        iteration_atoms.entry(iteration).or_insert_with(HashSet::new).insert(atom_id);
      }
    }
  }

  // Get sorted iterations
  let mut iterations: Vec<usize> = iteration_atoms.keys().copied().collect();
  iterations.sort();

  println!("Checking {}: found {} iterations with {} atoms expected per iteration",
           file_name, iterations.len(), num_atoms);

  // Check we have the right number of iterations (including initial state)
  let expected_total = expected_iterations + 1; // +1 for initial state
  assert_eq!(
    iterations.len(), expected_total,
    "{}: Expected {} iterations (1 initial + {} steps), found {}",
    file_name, expected_total, expected_iterations, iterations.len()
  );

  // Check for gaps in iterations
  if !iterations.is_empty() {
    for i in 0..iterations.len() - 1 {
      if iterations[i + 1] != iterations[i] + 1 {
        panic!(
          "{}: Gap detected between iterations {} and {}",
          file_name, iterations[i], iterations[i + 1]
        );
      }
    }
  }

  // Check that each iteration has all atoms (IDs 0 to num_atoms-1)
  for iteration in &iterations {
    let atom_ids = iteration_atoms.get(iteration).unwrap();

    // Check we have the right number of atoms
    assert_eq!(
      atom_ids.len(), num_atoms,
      "{}: Iteration {} has {} atoms, expected {}",
      file_name, iteration, atom_ids.len(), num_atoms
    );

    // Check all atom IDs from 0 to num_atoms-1 are present
    for expected_id in 0..num_atoms {
      assert!(
        atom_ids.contains(&expected_id),
        "{}: Iteration {} is missing atom with ID {}",
        file_name, iteration, expected_id
      );
    }

    // Check no invalid atom IDs (>= num_atoms)
    for &atom_id in atom_ids.iter() {
      assert!(
        atom_id < num_atoms,
        "{}: Iteration {} has invalid atom ID {} (should be < {})",
        file_name, iteration, atom_id, num_atoms
      );
    }
  }

  println!("✓ {}: All {} iterations have complete atom data (all atom IDs 0-{} present)",
           file_name, iterations.len(), num_atoms - 1);
}

/// Test 4: Test edge case with single reset - verifies to_transfer_struct returns only last batch
fn test_single_reset_runner(world_type: WorldType) {
  let simulation_size = Vector3::new(10., 10., 10.);
  let atom_factory = SafeAtomFactory::new();

  let mut atoms: Vec<Particle> = Vec::new();
  let atom = atom_factory.get_atom(
    AtomType::C,
    Vector3::new(5.0, 5.0, 5.0),
    Vector3::new(0.1, 0.0, 0.0)
  );
  atoms.push(atom);

  let num_iterations = 15;
  let max_iteration_till_reset = 10;
  let one_frame_duration = 1e-2;

  let save_options = SaveOptions {
    save: false,
    save_path: String::new(),
    save_laamps: false,
    save_verbose: false,
  };

  // TODO: go back to verlet
  let integration_algorithm = IntegrationAlgorithm::NoseHooverVerlet;

  let mut engine = Engine::new_from_atoms(
    atoms,
    simulation_size,
    TIME_STEP,
    num_iterations,
    max_iteration_till_reset,
    true, // save all iterations
    one_frame_duration,
    save_options,
    integration_algorithm.clone(),
    world_type.clone(),
  );

  let params = IntegrationAlgorithmParams::NoseHooverVerlet {
    desired_temperature: 100.,
    q_effective_mass: 100.,
  };
  engine.run(&params, TIME_STEP);

  let engine_dto = engine.to_transfer_struct();
  let world_dto = &engine_dto.world;

  assert!(validate_dto_type(&world_type, &world_dto), "WorldDTO type doesn't match WorldType");

  let atoms = match world_dto {
    WorldDTO::SimpleWorldDTO(simple_dto) => &simple_dto.atoms,
    WorldDTO::BoxedWorldDTO(boxed_dto) => &boxed_dto.box_container.atoms,
  };

  let mut all_iterations: Vec<usize> = Vec::new();

  for iteration_atoms in atoms {
    if !iteration_atoms.is_empty() {
      all_iterations.push(iteration_atoms[0].iteration);
    }
  }

  all_iterations.sort();

  // After reset at iteration 10, we should have iterations 10-15 (6 iterations including the duplicate at 10)
  // The first batch (0-9) was cleared during reset, only last atom from it (iteration 10) is kept
  let expected_iterations_in_memory = num_iterations - max_iteration_till_reset + 1; // 15 - 10 + 1 = 6
  assert_eq!(all_iterations.len(), expected_iterations_in_memory,
             "Wrong number of iterations in memory after reset");

  // Verify no gaps in the last batch
  for i in 0..all_iterations.len() - 1 {
    assert_eq!(all_iterations[i + 1], all_iterations[i] + 1, "Gap in iterations");
  }

  // Verify no duplicates
  let unique: HashSet<usize> = all_iterations.iter().copied().collect();
  assert_eq!(unique.len(), all_iterations.len(), "Duplicate iterations");

  // Verify the iterations start from the reset point
  assert_eq!(all_iterations[0], max_iteration_till_reset,
             "First iteration in memory should be at reset point");

  println!("✓ Single reset test passed: to_transfer_struct correctly returns only last batch");
}

#[test]
fn test_single_reset_simple_world() {
  test_single_reset_runner(WorldType::SimpleWorld);
}

#[test]
fn test_single_reset_boxed_world() {
  test_single_reset_runner(WorldType::BoxedWorld);
}

/// Test 5: Test no reset scenario
fn test_no_reset_runner(world_type: WorldType) {
  let simulation_size = Vector3::new(10., 10., 10.);
  let atom_factory = SafeAtomFactory::new();

  let mut atoms: Vec<Particle> = Vec::new();
  let atom = atom_factory.get_atom(
    AtomType::Fe,
    Vector3::new(5.0, 5.0, 5.0),
    Vector3::new(0.0, 0.1, 0.0)
  );
  atoms.push(atom);

  let num_iterations = 8;
  let max_iteration_till_reset = 20; // Won't trigger reset
  let one_frame_duration = 1e-2;

  let save_options = SaveOptions {
    save: false,
    save_path: String::new(),
    save_laamps: false,
    save_verbose: false,
  };

  let integration_algorithm = IntegrationAlgorithm::NoseHooverVerlet;

  let mut engine = Engine::new_from_atoms(
    atoms,
    simulation_size,
    TIME_STEP,
    num_iterations,
    max_iteration_till_reset,
    true, // save all iterations
    one_frame_duration,
    save_options,
    integration_algorithm.clone(),
    world_type.clone(),
  );

  let params = IntegrationAlgorithmParams::NoseHooverVerlet {
    desired_temperature: 100.,
    q_effective_mass: 100.,
  };
  engine.run(&params, TIME_STEP);

  let engine_dto = engine.to_transfer_struct();
  let world_dto = &engine_dto.world;

  assert!(validate_dto_type(&world_type, &world_dto), "WorldDTO type doesn't match WorldType");

  let atoms = match world_dto {
    WorldDTO::SimpleWorldDTO(simple_dto) => &simple_dto.atoms,
    WorldDTO::BoxedWorldDTO(boxed_dto) => &boxed_dto.box_container.atoms,
  };

  let mut all_iterations: Vec<usize> = Vec::new();

  for iteration_atoms in atoms {
    if !iteration_atoms.is_empty() {
      all_iterations.push(iteration_atoms[0].iteration);
    }
  }

  all_iterations.sort();

  let expected_total = num_iterations + 1;
  assert_eq!(all_iterations.len(), expected_total, "Wrong number of iterations");

  for i in 0..all_iterations.len() - 1 {
    assert_eq!(all_iterations[i + 1], all_iterations[i] + 1, "Gap in iterations");
  }

  let unique: HashSet<usize> = all_iterations.iter().copied().collect();
  assert_eq!(unique.len(), all_iterations.len(), "Duplicate iterations");

  println!("✓ No reset test passed");
}

#[test]
fn test_no_reset_simple_world() {
  test_no_reset_runner(WorldType::SimpleWorld);
}

#[test]
fn test_no_reset_boxed_world() {
  test_no_reset_runner(WorldType::BoxedWorld);
}

