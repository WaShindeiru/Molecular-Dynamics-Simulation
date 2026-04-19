pub mod builder;

use nalgebra::Vector3;
use crate::particle::Particle;
use crate::sim_core::world::integration::IntegrationAlgorithm;
use crate::sim_core::world::saver::SaveOptions;
use crate::sim_core::world::WorldType;
use crate::sim_core::world::boundary_constraint::EdgeCondition;
use crate::data::types::AtomType;

#[derive(Clone)]
pub struct SimulationConfig {
    pub atoms: Option<Vec<Particle>>,
    pub num_of_atoms: usize,
    pub num_of_carbon_atoms: usize,
    pub num_of_iron_atoms: usize,
    pub world_size: Vector3<f64>,
    pub potential_gravity_max: f64,
    pub time_step: f64,
    pub num_of_iterations: usize,
    pub max_iteration_till_reset: usize,
    pub save_all_iterations: bool,
    pub one_frame_duration: f64,
    pub save_options: SaveOptions,
    pub integration_algorithm: IntegrationAlgorithm,
    pub world_type: WorldType,
    pub edge_condition: EdgeCondition,
}

impl SimulationConfig {
    pub fn new(
        atoms: Vec<Particle>,
        world_size: Vector3<f64>,
        potential_gravity_max: f64,
        time_step: f64,
        num_of_iterations: usize,
        max_iteration_till_reset: usize,
        save_all_iterations: bool,
        one_frame_duration: f64,
        save_options: SaveOptions,
        integration_algorithm: IntegrationAlgorithm,
        world_type: WorldType,
        edge_condition: EdgeCondition,
    ) -> Self {

        let num_of_atoms = atoms.len();
        let mut num_of_carbon_atoms = 0;
        let mut num_of_iron_atoms = 0;

        for particle in &atoms {
            match particle.get_type() {
                AtomType::C => num_of_carbon_atoms += 1,
                AtomType::Fe => num_of_iron_atoms += 1,
            }
        }

        SimulationConfig {
            atoms: Some(atoms),
            num_of_atoms,
            num_of_carbon_atoms,
            num_of_iron_atoms,
            world_size,
            potential_gravity_max,
            time_step,
            num_of_iterations,
            max_iteration_till_reset,
            save_all_iterations,
            one_frame_duration,
            save_options,
            integration_algorithm,
            world_type,
            edge_condition,
        }
    }

    /// Create a config without atoms (for internal structures like BoxContainer)
    pub fn new_without_atoms(
        world_size: Vector3<f64>,
        num_of_atoms: usize,
        num_of_carbon_atoms: usize,
        num_of_iron_atoms: usize,
        potential_gravity_max: f64,
        time_step: f64,
        num_of_iterations: usize,
        max_iteration_till_reset: usize,
        save_all_iterations: bool,
        one_frame_duration: f64,
        save_options: SaveOptions,
        integration_algorithm: IntegrationAlgorithm,
        world_type: WorldType,
        edge_condition: EdgeCondition,
    ) -> Self {
        SimulationConfig {
            atoms: None,
            num_of_atoms,
            num_of_carbon_atoms,
            num_of_iron_atoms,
            world_size,
            potential_gravity_max,
            time_step,
            num_of_iterations,
            max_iteration_till_reset,
            save_all_iterations,
            one_frame_duration,
            save_options,
            integration_algorithm,
            world_type,
            edge_condition,
        }
    }

    pub fn count_particles_by_type(&self) -> (usize, usize) {
        (self.num_of_carbon_atoms, self.num_of_iron_atoms)
    }
}

