use nalgebra::Vector3;

use crate::data::config::builder::SimulationConfigBuilder;
use crate::data::types::AtomType;
use crate::data::units::{TEMPERATURE_U, ValueUnits};
use crate::persistence::json::{GeneratorConfigFile, SimulationConfigFile};

use crate::particle::SafeAtomFactory;

use crate::sim_core::world::boundary_constraint::EdgeCondition;
use crate::sim_core::world::integration::{TemperatureInfo, TimeIterationDistance};
use crate::simulations::generators::generator_config::dense::DenseGeneratorConfig;
use crate::simulations::generators::generator_config::GeneratorConfig;

pub fn see_config_json() {
  let simulation_size = Vector3::new(10., 10., 10.);
  let atom_factory = SafeAtomFactory::new(1., 10.);
  let atoms = vec![atom_factory.get_atom(
    AtomType::Fe,
    Vector3::new(5., 5., 5.),
    Vector3::new(0.01, 0.01, 0.01),
  )];

  let desired_temperatures = vec![
    TemperatureInfo::new(2000., TimeIterationDistance::Iteration { value: 20000 }),
    TemperatureInfo::new(200., TimeIterationDistance::Iteration { value: 2000 }),
  ]
  .into_iter()
  .map(|i| TemperatureInfo {
    desired_temperature: i.desired_temperature / TEMPERATURE_U,
    ..i
  })
  .collect();

  let config = SimulationConfigBuilder::new()
    .atoms(atoms)
    .world_size(simulation_size)
    .edge_condition(EdgeCondition::Periodic)
    .integration_algorithm(
      crate::sim_core::world::integration::IntegrationAlgorithm::NoseHooverVerlet {
        desired_temperature: (desired_temperatures),
        q_effective_mass: (1.0),
      },
    )
    .build()
    .unwrap();

  let result = SimulationConfigFile::from_runtime(&config, ValueUnits::Unitless)
    .to_json_string()
    .unwrap();

  println!("{}", result);
}

pub fn see_dense_generator_configuration() {
  let particle_distance = 7.;
  let offset = Vector3::new(1.7, 1.7, 1.7);
  let config_local = GeneratorConfig::Dense(DenseGeneratorConfig::new(particle_distance, offset));
  let config_file = GeneratorConfigFile::new(config_local, ValueUnits::Unitless).to_value_units(ValueUnits::Si);
  config_file.to_json_file("/media/washindeiru/7E442D59442D1585/md/temp/dense.json").expect("should work!");
}
