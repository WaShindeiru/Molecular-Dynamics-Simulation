use carbon_nanotube::persistence::json::particle_config::read_particle_config_from_json_file;
use carbon_nanotube::persistence::json::simulation_config::read_simulation_config_from_json_file;
use carbon_nanotube::sim_core::world::boxed_world::BoxedWorld;
use carbon_nanotube::sim_core::world::thermostat::IntegrationAlgorithm;

/// Loads `tests/fixtures/thermostat/v1` (a real `BoxedWorld` + Nose-Hoover thermostat config
/// with `acceptance_distance: 0` and `achieved_distance: 3`) and asserts that the desired
/// temperature switches every 5th iteration.
///
/// `acceptance_distance: 0` means the simulated temperature is never compared against the
/// desired one - the stage moves to `TemperatureAchieved` immediately. From there,
/// `achieved_distance: 3` (`consecutive > 3`) takes 4 more iterations to fire. The first stage
/// starts already in `StabilizeTemperature`, so its switch lands on iteration 5 (1 to become
/// achieved + 4 to count down). Every subsequent stage spends 1 iteration re-entering
/// `StabilizeTemperature` before immediately becoming achieved again, so it also takes 5
/// iterations - giving a steady 5-iteration cadence between all switches.
#[test]
fn desired_temperature_changes_every_fifth_iteration() {
  let manifest_dir = env!("CARGO_MANIFEST_DIR");
  let parameters_path = format!("{}/tests/fixtures/thermostat/v1/parameters.json", manifest_dir);
  let particles_path = format!(
    "{}/tests/fixtures/thermostat/v1/particles_initial.json",
    manifest_dir
  );

  let mut simulation_config = read_simulation_config_from_json_file(&parameters_path)
    .expect("failed to parse thermostat parameters.json fixture");
  let particle_config = read_particle_config_from_json_file(&particles_path)
    .expect("failed to parse thermostat particles_initial.json fixture");

  // The fixture's save_path points at a machine-specific location; disable saving so the
  // test does not depend on (or write to) that path.
  simulation_config.save_options.save = false;

  let num_temperature_steps = match &simulation_config.integration_algorithm {
    IntegrationAlgorithm::NoseHooverVerlet {
      desired_temperature,
      ..
    } => desired_temperature.len(),
    other => panic!("expected NoseHooverVerlet integration algorithm, got {:?}", other),
  };

  let num_of_iterations = simulation_config.num_of_iterations;
  let mut world = BoxedWorld::with_config(simulation_config.clone(), particle_config);

  for iteration in 1..=num_of_iterations {
    world
      .update(
        &simulation_config.integration_algorithm,
        simulation_config.time_step,
        iteration,
      )
      .expect("world update failed");
  }

  let history = world
    .integration_algorithm_state()
    .temperature_history()
    .expect("expected NoseHooverVerlet temperature history");

  let switch_iterations: Vec<usize> = history
    .iter()
    .filter_map(|entry| entry.temperature_switched.map(|s| s.iteration))
    .collect();

  println!(
    "Desired temperature switched {} time(s) over {} iterations: {:?}",
    switch_iterations.len(),
    num_of_iterations,
    switch_iterations
  );

  // The last temperature entry never switches away, so at most `num_temperature_steps - 1`
  // switches can happen, each 5 iterations apart.
  let max_switches = num_temperature_steps - 1;
  let expected_switch_iterations: Vec<usize> = (1..=max_switches)
    .map(|n| n * 5)
    .take_while(|&iteration| iteration <= num_of_iterations)
    .collect();

  assert_eq!(
    switch_iterations, expected_switch_iterations,
    "expected desired temperature to switch every 5th iteration"
  );
}
