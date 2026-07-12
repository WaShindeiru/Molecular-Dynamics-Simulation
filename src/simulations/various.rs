use nalgebra::Vector3;
use crate::sim_core::world::boxed_world::box_task::task_manager::TaskManagerConfig;
use crate::sim_core::world::WorldType;
use crate::data::ValueUnits::{Si, Unitless};
use crate::data::config::builder::SimulationConfigBuilder;
use crate::data::types::AtomType;
use crate::data::units::{TEMPERATURE_U, ValueUnits};
use crate::persistence::json::integration_algorithm_file::IntegrationAlgorithmFile;
use crate::persistence::json::simple_temperature_info_generator_file::SimpleTemperatureInfoGeneratorFile;
use crate::persistence::json::temperature_info_file::TemperatureInfoFile;
use crate::persistence::json::temperature_info_source_file::TemperatureInfoSourceFile;
use crate::persistence::json::{GeneratorConfigFile, SimulationConfigFile};

use crate::particle::SafeAtomFactory;

use crate::sim_core::world::boundary_constraint::EdgeCondition;
use crate::sim_core::world::thermostat::{TemperatureInfo, TimeIterationDistance};
use crate::persistence::json::particle_config::Vector3Record;
use crate::sim_core::world::cell::TaskSplitVariant::FloorBox;
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
    .world_type(
      WorldType::LinkedCellWorld {
        task_manager_config: TaskManagerConfig{
          debug: false,
          split: FloorBox {
            x: 2,
            y: 2,
          },
          task_worker_multiplier: 3.0
        }
      }
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

pub fn see_temperature_info_generator_config() {
  let to_unitless = |info: TemperatureInfo| TemperatureInfo {
    desired_temperature: info.desired_temperature / TEMPERATURE_U,
    ..info
  };

  let high_temperature = to_unitless(TemperatureInfo::new(
    2000.,
    TimeIterationDistance::Iteration { value: 20_000 },
  ));
  let low_temperature = to_unitless(TemperatureInfo::new(
    200.,
    TimeIterationDistance::Iteration { value: 2_000 },
  ));

  let desired_temperature = vec![
    TemperatureInfoSourceFile::Direct(TemperatureInfoFile::from_runtime(&high_temperature)),
    TemperatureInfoSourceFile::Simple(SimpleTemperatureInfoGeneratorFile {
      start_temperature: 2000. / TEMPERATURE_U,
      end_temperature: 500. / TEMPERATURE_U,
      acceptance_distance: TimeIterationDistance::Iteration { value: 0 },
      achieved_distance: TimeIterationDistance::Iteration { value: 5_000 },
      temperature_step: Some(-500.0 / TEMPERATURE_U),
      threshold: Some(60.0 / TEMPERATURE_U),
      save_step: Some(3),
    }),
    TemperatureInfoSourceFile::Direct(TemperatureInfoFile::from_runtime(&low_temperature)),
  ];

  let algorithm = IntegrationAlgorithmFile::NoseHooverVerlet {
    desired_temperature,
    q_effective_mass: 1.0,
  };

  let json = serde_json::to_string_pretty(&algorithm.to_value_units(Unitless, Si)).unwrap();
  println!("{json}");

  let expanded = algorithm.to_runtime().unwrap();
  if let crate::sim_core::world::thermostat::IntegrationAlgorithm::NoseHooverVerlet {
    desired_temperature,
    ..
  } = expanded
  {
    println!("\nExpanded temperature schedule:");
    for (index, info) in desired_temperature.iter().enumerate() {
      println!(
        "  [{index}] {:.1} K (save={})",
        info.desired_temperature * TEMPERATURE_U,
        info.save,
      );
    }
  }
}