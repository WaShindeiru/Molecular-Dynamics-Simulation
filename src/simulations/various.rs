use crate::{data::{SimulationConfig, config::builder::SimulationConfigBuilder, types::AtomType, units::TEMPERATURE_U}, particle::SafeAtomFactory, sim_core::world::integration::{TemperatureInfo, TimeIterationDistance}};
use nalgebra::Vector3;
use crate::sim_core::world::boundary_constraint::EdgeCondition;

pub fn see_config_json() {
  let simulation_size = Vector3::new(10., 10., 10.);
  let atom_factory = SafeAtomFactory::new(1., 10.);
  let atoms = vec![atom_factory.get_atom(AtomType::Fe, Vector3::new(5., 5., 5.), Vector3::new(0.01, 0.01, 0.01))];

  let desired_temperatures = vec![
    TemperatureInfo::new(2000., TimeIterationDistance::Iteration { value: 20000 }),
    TemperatureInfo::new(200., TimeIterationDistance::Iteration { value: 2000 }),
  ].into_iter().map(
    |i| TemperatureInfo{
      desired_temperature: i.desired_temperature / TEMPERATURE_U,
      ..i
    }).collect();

  let config = SimulationConfigBuilder::new()
    .atoms(atoms)
    .world_size(simulation_size)
    .edge_condition(EdgeCondition::Periodic)
    .integration_algorithm(crate::sim_core::world::integration::IntegrationAlgorithm::NoseHooverVerlet { desired_temperature: (desired_temperatures), q_effective_mass: (1.0) })
    .build().unwrap();
  
  let result = config.to_json_string().unwrap();

  println!("{}", result);
}
