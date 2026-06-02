use nalgebra::Vector3;

use crate::data::config::builder::SimulationConfigBuilder;
use crate::data::types::AtomType;
use crate::data::units::{TEMPERATURE_U, ValueUnits};
use crate::persistence::json::{GeneratorConfigFile, SimulationConfigFile};

use crate::particle::SafeAtomFactory;

use crate::sim_core::world::boundary_constraint::EdgeCondition;
use crate::sim_core::world::thermostat::{TemperatureInfo, TimeIterationDistance};
use crate::persistence::json::particle_config::Vector3Record;
use crate::simulations::generators::core::generator_config::GeneratorConfig;
use crate::simulations::generators::core::generator_config::dense::DenseGeneratorConfig;
use crate::simulations::generators::core::generator_config::nanotube::{NanotubeGeneratorConfig, NanotubeGeneratorParticleFile};

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
    .edge_condition(EdgeCondition::Periodic {
      trigger_small_subtask_size: 1,
      split: EdgeCondition::DEFAULT_SPLIT,
    })
    .integration_algorithm(
      crate::sim_core::world::thermostat::IntegrationAlgorithm::NoseHooverVerlet {
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
  let config_file =
    GeneratorConfigFile::new(config_local, ValueUnits::Unitless).to_value_units(ValueUnits::Si);
  config_file
    .to_json_file("/media/washindeiru/7E442D59442D1585/md/temp/dense.json")
    .expect("should work!");
}

pub fn see_nanotube_generator_configuration() {
  let particles = vec![
    NanotubeGeneratorParticleFile {
      position: Vector3Record::from(Vector3::new(3., 4., 5.)),
      particle_type: AtomType::C_nanotube,
    },
  ];

  let offset = Vector3Record::from(Vector3::new(1., 1., 1.));
  let vel_mean = Vector3Record::from(Vector3::new(1., 1., 1.));
  let vel_std_dev = Vector3Record::from(Vector3::new(1., 1., 1.));

  let config_local = GeneratorConfig::Nanotube(NanotubeGeneratorConfig::new(
    particles, offset, vel_mean, vel_std_dev
  ));
  let config_file = GeneratorConfigFile::new(config_local, ValueUnits::Unitless).to_value_units(ValueUnits::Unitless);
  let result = config_file.to_json_string().unwrap();
  println!("{result}");
}
