use carbon_nanotube::data::ParticleConfig;
use carbon_nanotube::data::config::builder::SimulationConfigBuilder;
use carbon_nanotube::data::types::AtomType;
use carbon_nanotube::persistence::json::particle_config::read_particle_config_from_json_file;
use carbon_nanotube::sim_core::world::WorldType;
use carbon_nanotube::sim_core::world::boundary_constraint::EdgeCondition;
use carbon_nanotube::sim_core::world::boxed_world::box_task::task_manager::TaskManagerConfig;
use carbon_nanotube::sim_core::world::linked_cell_world::{LinkedCellContainer, LinkedCellWorld};
use carbon_nanotube::sim_core::world::thermostat::{
    IntegrationAlgorithm, TemperatureInfo, TimeIterationDistance,
};
use nalgebra::Vector3;

const NUM_ITERATIONS: usize = 5;
const TIME_STEP: f64 = 1e-3;
const WORLD_SIZE: Vector3<f64> = Vector3::new(30.0, 30.0, 30.0);

fn load_world() -> (LinkedCellWorld, IntegrationAlgorithm, Vec<(usize, AtomType, f64)>) {
    let fixture_path =
        concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/particles_initial.json");
    let particle_config =
        read_particle_config_from_json_file(fixture_path).expect("failed to read fixture");

    let initial_data: Vec<(usize, AtomType, f64)> = particle_config
        .atoms
        .iter()
        .map(|p| (p.get_id(), p.get_type(), p.get_mass()))
        .collect();

    let integration_algorithm = IntegrationAlgorithm::NoseHooverVerlet {
        desired_temperature: vec![TemperatureInfo::new(
            300.0,
            TimeIterationDistance::Iteration {
                value: NUM_ITERATIONS,
            },
        )],
        q_effective_mass: 100.0,
    };

    let config = SimulationConfigBuilder::new()
        .atoms(particle_config.atoms.clone())
        .world_size(WORLD_SIZE)
        .time_step(TIME_STEP)
        .num_of_iterations(NUM_ITERATIONS)
        .max_iteration_till_reset(100)
        .integration_algorithm(integration_algorithm.clone())
        .world_type(WorldType::LinkedCellWorld {
            task_manager_config: TaskManagerConfig {
                debug: false,
                task_worker_multiplier: 2.0,
            },
        })
        .edge_condition(EdgeCondition::PeriodicAll)
        .build()
        .unwrap();

    let world = LinkedCellWorld::with_config(config, ParticleConfig::new(particle_config.atoms));

    (world, integration_algorithm, initial_data)
}

fn check_vec_index_matches_id(container: &LinkedCellContainer, label: &str) {
    for (i, slot) in container.particles().iter().enumerate() {
        if let Some(particle) = slot {
            assert_eq!(
                particle.get_id(),
                i,
                "{}: particles[{}] holds particle with id {}",
                label,
                i,
                particle.get_id()
            );
        }
    }
}

fn check_type_and_mass(
    container: &LinkedCellContainer,
    initial_data: &[(usize, AtomType, f64)],
    label: &str,
) {
    for &(id, initial_type, initial_mass) in initial_data {
        let particle = container.particles()[id].as_ref().unwrap_or_else(|| {
            panic!("{}: particle with id {} is missing from particles[{}]", label, id, id)
        });
        assert_eq!(
            particle.get_type(),
            initial_type,
            "{}: particle id {} type changed",
            label,
            id
        );
        assert_eq!(
            particle.get_mass(),
            initial_mass,
            "{}: particle id {} mass changed",
            label,
            id
        );
    }
}

#[test]
fn test_linked_cell_world_particle_invariants_at_checkpoints() {
    let (mut world, algorithm, initial_data) = load_world();

    for iteration in 0..NUM_ITERATIONS {
        let before_label = format!("iteration {} / before", iteration + 1);
        let half_vel_label = format!("iteration {} / after half-velocity step", iteration + 1);
        let end_label = format!("iteration {} / end", iteration + 1);

        let initial_data_ref = &initial_data;

        world.update_verlet_nose_hoover_instrumented(
            iteration + 1,
            |container| {
                check_vec_index_matches_id(container, &before_label);
                check_type_and_mass(container, initial_data_ref, &before_label);
            },
            |container| {
                check_vec_index_matches_id(container, &half_vel_label);
                check_type_and_mass(container, initial_data_ref, &half_vel_label);
            },
            |container| {
                check_vec_index_matches_id(container, &end_label);
                check_type_and_mass(container, initial_data_ref, &end_label);
            },
        );
    }
}
