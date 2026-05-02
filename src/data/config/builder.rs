use nalgebra::Vector3;
use crate::data::SimulationConfig;
use crate::data::types::AtomType;
use crate::particle::Particle;
use crate::sim_core::world::boundary_constraint::EdgeCondition;
use crate::sim_core::world::integration::{IntegrationAlgorithm, TemperatureInfo, TimeIterationDistance};
use crate::sim_core::world::saver::SaveOptions;
use crate::sim_core::world::WorldType;

pub struct SimulationConfigBuilder {
	atoms: Option<Vec<Particle>>,
	world_size: Option<Vector3<f64>>,
	potential_gravity_max: Option<f64>,
	time_step: Option<f64>,
	num_of_iterations: Option<usize>,
	max_iteration_till_reset: Option<usize>,
	save_all_iterations: Option<bool>,
	one_frame_duration: Option<f64>,
	save_options: Option<SaveOptions>,
	integration_algorithm: Option<IntegrationAlgorithm>,
	world_type: Option<WorldType>,
	edge_condition: Option<EdgeCondition>,
}

impl SimulationConfigBuilder {
	pub fn new() -> Self {
		SimulationConfigBuilder {
			atoms: None,
			world_size: None,
			potential_gravity_max: None,
			time_step: None,
			num_of_iterations: None,
			max_iteration_till_reset: None,
			save_all_iterations: None,
			one_frame_duration: None,
			save_options: None,
			integration_algorithm: None,
			world_type: None,
			edge_condition: None,
		}
	}

	pub fn atoms(mut self, atoms: Vec<Particle>) -> Self {
		self.atoms = Some(atoms);
		self
	}

	pub fn world_size(mut self, world_size: Vector3<f64>) -> Self {
		self.world_size = Some(world_size);
		self
	}

	pub fn potential_gravity_max(mut self, potential_gravity_max: f64) -> Self {
		self.potential_gravity_max = Some(potential_gravity_max);
		self
	}

	pub fn time_step(mut self, time_step: f64) -> Self {
		self.time_step = Some(time_step);
		self
	}

	pub fn num_of_iterations(mut self, num_of_iterations: usize) -> Self {
		self.num_of_iterations = Some(num_of_iterations);
		self
	}

	pub fn max_iteration_till_reset(mut self, max_iteration_till_reset: usize) -> Self {
		self.max_iteration_till_reset = Some(max_iteration_till_reset);
		self
	}

	pub fn save_all_iterations(mut self, save_all_iterations: bool) -> Self {
		self.save_all_iterations = Some(save_all_iterations);
		self
	}

	pub fn one_frame_duration(mut self, one_frame_duration: f64) -> Self {
		self.one_frame_duration = Some(one_frame_duration);
		self
	}

	pub fn save_options(mut self, save_options: SaveOptions) -> Self {
		self.save_options = Some(save_options);
		self
	}

	pub fn integration_algorithm(mut self, integration_algorithm: IntegrationAlgorithm) -> Self {
		self.integration_algorithm = Some(integration_algorithm);
		self
	}

	pub fn world_type(mut self, world_type: WorldType) -> Self {
		self.world_type = Some(world_type);
		self
	}

	pub fn edge_condition(mut self, edge_condition: EdgeCondition) -> Self {
		self.edge_condition = Some(edge_condition);
		self
	}

	pub fn assert_all_set(self) -> Result<Self, String> {
		let mut missing_fields = Vec::new();

		if self.atoms.is_none() {
			missing_fields.push("atoms");
		}
		if self.world_size.is_none() {
			missing_fields.push("world_size");
		}
		if self.potential_gravity_max.is_none() {
			missing_fields.push("potential_gravity_max");
		}
		if self.time_step.is_none() {
			missing_fields.push("time_step");
		}
		if self.num_of_iterations.is_none() {
			missing_fields.push("num_of_iterations");
		}
		if self.max_iteration_till_reset.is_none() {
			missing_fields.push("max_iteration_till_reset");
		}
		if self.save_all_iterations.is_none() {
			missing_fields.push("save_all_iterations");
		}
		if self.one_frame_duration.is_none() {
			missing_fields.push("one_frame_duration");
		}
		if self.save_options.is_none() {
			missing_fields.push("save_options");
		}
		if self.integration_algorithm.is_none() {
			missing_fields.push("integration_algorithm");
		}
		if self.world_type.is_none() {
			missing_fields.push("world_type");
		}
		if self.edge_condition.is_none() {
			missing_fields.push("edge_condition");
		}

		if missing_fields.is_empty() {
			Ok(self)
		} else {
			Err(format!("The following fields are not set: {}", missing_fields.join(", ")))
		}
	}

	pub fn get_missing_fields(&self) -> Vec<&'static str> {
		let mut missing = Vec::new();

		if self.atoms.is_none() {
			missing.push("atoms");
		}
		if self.world_size.is_none() {
			missing.push("world_size");
		}
		if self.potential_gravity_max.is_none() {
			missing.push("potential_gravity_max");
		}
		if self.time_step.is_none() {
			missing.push("time_step");
		}
		if self.num_of_iterations.is_none() {
			missing.push("num_of_iterations");
		}
		if self.max_iteration_till_reset.is_none() {
			missing.push("max_iteration_till_reset");
		}
		if self.save_all_iterations.is_none() {
			missing.push("save_all_iterations");
		}
		if self.one_frame_duration.is_none() {
			missing.push("one_frame_duration");
		}
		if self.save_options.is_none() {
			missing.push("save_options");
		}
		if self.integration_algorithm.is_none() {
			missing.push("integration_algorithm");
		}
		if self.world_type.is_none() {
			missing.push("world_type");
		}
		if self.edge_condition.is_none() {
			missing.push("edge_condition");
		}

		missing
	}

	pub fn build(self) -> Result<SimulationConfig, String> {
		// Atoms are required - world_size also required to make sense
		let atoms = self.atoms.ok_or("atoms must be specified")?;
		let world_size = self.world_size.ok_or("world_size must be specified")?;

		let num_of_atoms = atoms.len();
		let mut num_of_carbon_atoms = 0;
		let mut num_of_iron_atoms = 0;

		for particle in &atoms {
			match particle.get_type() {
				AtomType::C => num_of_carbon_atoms += 1,
				AtomType::Fe => num_of_iron_atoms += 1,
			}
		}

		Ok(SimulationConfig {
			atoms: Some(atoms),
			num_of_atoms,
			num_of_carbon_atoms,
			num_of_iron_atoms,
			world_size,
			potential_gravity_max: self.potential_gravity_max.unwrap_or(1.0),
			time_step: self.time_step.unwrap_or(1e-17),
			num_of_iterations: self.num_of_iterations.unwrap_or(1000),
			max_iteration_till_reset: self.max_iteration_till_reset.unwrap_or(1000),
			save_all_iterations: self.save_all_iterations.unwrap_or(false),
			one_frame_duration: self.one_frame_duration.unwrap_or(1e-16),
			save_options: self.save_options.unwrap_or(SaveOptions {
				save: false,
				save_laamps: false,
				save_verbose: false,
				save_path: String::from("output"),
			}),
			integration_algorithm: self.integration_algorithm.unwrap_or(IntegrationAlgorithm::NoseHooverVerlet {
				q_effective_mass: 1.,
				desired_temperature: vec![TemperatureInfo{desired_temperature: 2000., distance: TimeIterationDistance::Iteration { value: 20000 }}]
			}),
			world_type: self.world_type.unwrap_or(WorldType::BoxedWorld { task_worker_multiplier: 4.0 }),
			edge_condition: self.edge_condition.unwrap_or(EdgeCondition::Periodic),
		})
	}
}