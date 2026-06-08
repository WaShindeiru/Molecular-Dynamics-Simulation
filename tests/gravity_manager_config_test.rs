use carbon_nanotube::persistence::json::simulation_config::read_simulation_config_from_json_file;

#[test]
fn gravity_schedule_parsed_from_json_config() {
    let manifest_dir = env!("CARGO_MANIFEST_DIR");
    let path = format!("{}/tests/fixtures/gravity_test_config.json", manifest_dir);

    let config = read_simulation_config_from_json_file(&path)
        .expect("failed to parse gravity_test_config.json");

    let mut gravity_manager = config.gravity_manager.lock().unwrap();

    // Iteration 0: initial gravity entry (type: Iteration, value: 0) → 0.001
    assert_eq!(
        gravity_manager.get_gravity(0),
        0.001,
        "expected gravity 0.001 at iteration 0"
    );

    // Time-based entry: 34e-13 s / 1.75e-18 s = 1942857.14... → rounds to 1942857
    assert_eq!(
        gravity_manager.get_gravity(1942857),
        0.0,
        "expected gravity 0.0 at iteration 1942857"
    );
}
