use std::path::Path;
use std::{fs, io};

use log::info;
use nalgebra::Vector3;

use crate::data::types::AtomType;
use crate::data::ValueUnits;
use crate::particle::particle::ParticleKind;
use crate::persistence::dto::world::boxed::BoxedWorldDTO;
use crate::persistence::json::particle_config::particle_type_file::ParticleTypeFile;
use crate::persistence::json::particle_config::{ParticleConfigFile, ParticleInitialState};
use crate::sim_core::world::saver::{PartialWorldSaver, PeriodicSave, csv_writer_with_header};
use crate::sim_core::world::thermostat::IntegrationAlgorithm;

impl PartialWorldSaver {
  pub fn persist_boxed_world(&mut self, world: &BoxedWorldDTO) -> io::Result<()> {
    info!("Saving batch: {} ...", world.number_of_resets);
    let save_dir = Path::new(&self.save_options.save_path);
    fs::create_dir_all(save_dir)?;

    self.append_energies_boxed_world(world)?;

    if self.save_options.save_laamps {
      self.save_laamps_boxed_world(world)?;
    }

    let (forces, potential_energies) = self.compute_force_and_potential_boxed_world(world);

    if self.save_options.save_verbose {
      self.append_forces_boxed_world(world, &forces)?;
      self.append_potential_energies_boxed_world(world, &potential_energies)?;
      self.append_positions_in_one_file_boxed_world(world)?;
    }

    self.save_periodic_particles_boxed_world(world)?;

    info!("Batch {} saved", world.number_of_resets);

    Ok(())
  }

  fn save_periodic_particles_boxed_world(&mut self, world: &BoxedWorldDTO) -> io::Result<()> {
    let iteration_distance = match self.save_options.periodic_save {
      PeriodicSave::Disabled => return Ok(()),
      PeriodicSave::Enabled { iteration_distance } => iteration_distance,
    };

    let dir = Path::new(&self.save_options.save_path)
      .join("particles")
      .join("periodic");
    fs::create_dir_all(&dir)?;

    for (i, box_container) in world.history.box_container.iter().enumerate() {
      let should_save = self.periodic_save_iteration_count % iteration_distance == 0;
      self.periodic_save_iteration_count += 1;

      if !should_save {
        continue;
      }

      let global_iteration = world.number_of_resets * world.max_iteration_till_reset + i;

      let mut c_count = 0usize;
      let mut fe_count = 0usize;
      let mut particles = Vec::with_capacity(box_container.atoms.len());

      for atom in box_container.atoms.iter() {
        match atom.atom_type {
          AtomType::C | AtomType::C_nanotube | AtomType::C_nanotube_static => c_count += 1,
          AtomType::Fe => fe_count += 1,
        }
        particles.push(ParticleInitialState::new(
          atom.id,
          atom.atom_type,
          ParticleTypeFile::from(atom.kind),
          atom.position,
          atom.velocity,
          atom.velocity_manager_id,
          atom.control_velocity_manager_id,
        ));
      }

      let config_file = ParticleConfigFile {
        value_units: ValueUnits::Unitless,
        particles,
        num_of_atoms: c_count + fe_count,
        num_of_carbon_atoms: c_count,
        num_of_iron_atoms: fe_count,
        velocity_managers: self.velocity_managers_file.clone(),
        control_velocity_managers: self.control_velocity_managers_file.clone(),
      }
      .to_value_units(ValueUnits::Si);

      let json = serde_json::to_string_pretty(&config_file)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

      fs::write(dir.join(format!("{global_iteration}.json")), json)?;
    }

    Ok(())
  }

  pub fn append_energies_boxed_world(&mut self, world: &BoxedWorldDTO) -> io::Result<()> {
    let box_containers = &world.history.box_container;
    let num_of_iterations = world.num_of_world_iterations;

    let mut kinetic_energy_atom: Vec<f64> = Vec::with_capacity(num_of_iterations);
    let mut kinetic_energy_other: Vec<f64> = Vec::with_capacity(num_of_iterations);
    let mut potential_energy: Vec<f64> = Vec::with_capacity(num_of_iterations);
    let mut potential_gravity_energy: Vec<f64> = Vec::with_capacity(num_of_iterations);
    let mut total_energy: Vec<f64> = Vec::with_capacity(num_of_iterations);
    let mut thermostat_work: Vec<f64> = Vec::with_capacity(num_of_iterations);
    let mut p_control_energy: Vec<f64> = Vec::with_capacity(num_of_iterations);

    for (i, box_container_dto) in box_containers.iter().enumerate() {
      let mut kinetic_energy_atom_i = 0.;
      let mut kinetic_energy_other_i = 0.;
      let mut thermostat_work_i = 0.;
      let mut potential_energy_i = 0.;
      let mut potential_gravity_energy_i = 0.;
      let mut p_control_energy_i = 0.;

      for atom_dto in box_container_dto.atoms.iter() {
        if atom_dto.kind == ParticleKind::Atom {
          kinetic_energy_atom_i += atom_dto.kinetic_energy;
        } else {
          kinetic_energy_other_i += atom_dto.kinetic_energy;
        }
        thermostat_work_i += atom_dto.thermostat_work;
        potential_energy_i += atom_dto.potential_energy;
        potential_gravity_energy_i += atom_dto.potential_gravity_energy;
        p_control_energy_i += atom_dto.p_control_energy;
      }

      self.thermostat_work_total += thermostat_work_i;
      self.p_control_energy_total += p_control_energy_i;

      kinetic_energy_atom.push(kinetic_energy_atom_i);
      kinetic_energy_other.push(kinetic_energy_other_i);
      potential_energy.push(potential_energy_i);
      let potential_energy_i = *potential_energy.get(i).unwrap();
      potential_gravity_energy.push(potential_gravity_energy_i);
      total_energy.push(
        kinetic_energy_atom_i + kinetic_energy_other_i + potential_energy_i + potential_gravity_energy_i,
      );
      thermostat_work.push(self.thermostat_work_total);
      p_control_energy.push(self.p_control_energy_total);
    }

    let save_dir = Path::new(&self.save_options.save_path);
    let energy_header = match world.integration_algorithm {
      IntegrationAlgorithm::NoseHooverVerlet { .. } => &[
        "iteration",
        "kinetic_energy_atom",
        "kinetic_energy_other",
        "potential_energy",
        "potential_gravity_energy",
        "total_energy",
        "p_control_energy_total",
        "thermostat_work_total",
        "thermostat_epsilon",
      ][..],
      _ => &[
        "iteration",
        "kinetic_energy_atom",
        "kinetic_energy_other",
        "potential_energy",
        "potential_gravity_energy",
        "total_energy",
        "p_control_energy_total",
      ][..],
    };
    let mut wtr = csv_writer_with_header(&save_dir.join("energy.csv"), energy_header)?;

    assert!(
      kinetic_energy_atom.len() == potential_energy.len()
        && kinetic_energy_atom.len() == total_energy.len()
    );

    for i in 0..kinetic_energy_atom.len() {
      let should_save = self.energy_frame_iteration_count_current_iteration
        % world.energy_frame_iteration_count
        == 0;
      self.energy_frame_iteration_count_current_iteration += 1;

      if !should_save {
        continue;
      }

      let iteration = world.number_of_resets * world.max_iteration_till_reset + i;

      match world.integration_algorithm {
        IntegrationAlgorithm::NoseHooverVerlet { .. } => {
          wtr.write_record(&[
            format!("{}", iteration),
            format!("{}", kinetic_energy_atom.get(i).unwrap()),
            format!("{}", kinetic_energy_other.get(i).unwrap()),
            format!("{}", potential_energy.get(i).unwrap()),
            format!("{}", potential_gravity_energy.get(i).unwrap()),
            format!("{}", total_energy.get(i).unwrap()),
            format!("{}", p_control_energy.get(i).unwrap()),
            format!("{}", thermostat_work.get(i).unwrap()),
            format!("{}", world.history.thermostat_epsilon.get(i).unwrap()),
          ])?;
        }
        _ => {
          wtr.write_record(&[
            format!("{}", iteration),
            format!("{}", kinetic_energy_atom.get(i).unwrap()),
            format!("{}", kinetic_energy_other.get(i).unwrap()),
            format!("{}", potential_energy.get(i).unwrap()),
            format!("{}", potential_gravity_energy.get(i).unwrap()),
            format!("{}", total_energy.get(i).unwrap()),
            format!("{}", p_control_energy.get(i).unwrap()),
          ])?;
        }
      }
    }

    wtr.flush()?;

    Ok(())
  }

  fn save_laamps_boxed_world(&mut self, world: &BoxedWorldDTO) -> io::Result<()> {
    let box_containers = &world.history.box_container;
    let mut files_saved = 0;
    let save_dir = Path::new(&self.save_options.save_path);
    let laamps_dir = save_dir.join("laamps");
    fs::create_dir_all(&laamps_dir)?;

    for i in 0..box_containers.len() {
      if self.laamps_frame_iteration_count_current_iteration % world.laamps_frame_iteration_count
        == 0
      {
        let iteration_number = i + world.number_of_resets * world.max_iteration_till_reset;

        let mut result_string = String::new();
        result_string.push_str(&"ITEM: TIMESTEP\n".to_string());
        result_string.push_str(&format!("{}\n", iteration_number));

        result_string.push_str(&"ITEM: NUMBER OF ATOMS\n".to_string());
        let atom_container = &box_containers.get(i).unwrap().atoms;
        let num_of_atoms = atom_container.len();
        result_string.push_str(&format!("{}\n", num_of_atoms));

        result_string.push_str(&"ITEM: BOX BOUNDS pp pp pp\n".to_string());
        result_string.push_str(&format!("0.0 {}\n", world.size.x));
        result_string.push_str(&format!("0.0 {}\n", world.size.y));
        result_string.push_str(&format!("0.0 {}\n", world.size.z));

        result_string.push_str(&"ITEM: ATOMS id type x y z\n".to_string());

        for atom_dto in atom_container.iter() {
          result_string.push_str(&format!(
            "{} {} {} {} {}\n",
            atom_dto.id,
            atom_dto.atom_type as u64,
            atom_dto.position.x,
            atom_dto.position.y,
            atom_dto.position.z
          ));
        }

        fs::write(
          laamps_dir.join(format!("output_{}.dump", iteration_number)),
          result_string,
        )?;
        files_saved += 1;
      }

      self.laamps_frame_iteration_count_current_iteration += 1
    }

    info!("Saved {} laamps output files.", { files_saved });
    Ok(())
  }

  fn compute_force_and_potential_boxed_world(
    &self,
    world: &BoxedWorldDTO,
  ) -> (Vec<Vec<Vector3<f64>>>, Vec<Vec<f64>>) {
    let box_containers = &world.history.box_container;
    let mut forces: Vec<Vec<Vector3<f64>>> = Vec::new();
    let mut potential_energies: Vec<Vec<f64>> = Vec::new();

    for box_container_dto in box_containers.iter() {
      let num_atoms = box_container_dto.atoms.len();
      let mut current_forces: Vec<Vector3<f64>> = vec![Vector3::new(0., 0., 0.); num_atoms];
      let mut current_potential_energies: Vec<f64> = vec![0.; num_atoms];

      // Sort by atom ID to get a consistent local ordering that is independent
      // of the global ID values (which may not start at 0 across test runs).
      let mut sorted_atoms: Vec<_> = box_container_dto.atoms.iter().collect();
      sorted_atoms.sort_by_key(|a| a.id);

      for (local_idx, atom_dto) in sorted_atoms.iter().enumerate() {
        let force_i = current_forces.get_mut(local_idx).unwrap();
        force_i.x = atom_dto.force.x;
        force_i.y = atom_dto.force.y;
        force_i.z = atom_dto.force.z;

        *current_potential_energies.get_mut(local_idx).unwrap() = atom_dto.potential_energy;
      }

      forces.push(current_forces);
      potential_energies.push(current_potential_energies);
    }

    (forces, potential_energies)
  }

  // TODO: this method is exactly the same for simple world, refactor this please
  fn append_forces_boxed_world(
    &self,
    world: &BoxedWorldDTO,
    forces: &Vec<Vec<Vector3<f64>>>,
  ) -> io::Result<()> {
    let save_dir = Path::new(&self.save_options.save_path);
    let mut wtr = csv_writer_with_header(
      &save_dir.join("forces.csv"),
      &[
        "iteration",
        "particle_index",
        "force_x",
        "force_y",
        "force_z",
      ],
    )?;

    for i in 0..forces.len() {
      let iteration = world.number_of_resets * world.max_iteration_till_reset + i;

      let force_map = forces.get(i).unwrap();

      for (id, force_vector) in force_map.iter().enumerate() {
        wtr.write_record(&[
          format!("{}", iteration),
          format!("{}", id),
          format!("{}", force_vector.x),
          format!("{}", force_vector.y),
          format!("{}", force_vector.z),
        ])?;
      }
    }

    wtr.flush()?;

    Ok(())
  }

  // TODO: this method is the same as well
  fn append_potential_energies_boxed_world(
    &self,
    world: &BoxedWorldDTO,
    potential_energies: &Vec<Vec<f64>>,
  ) -> io::Result<()> {
    let save_dir = Path::new(&self.save_options.save_path);
    let mut wtr = csv_writer_with_header(
      &save_dir.join("potential_energies.csv"),
      &["iteration", "particle_index", "potential_energy"],
    )?;

    for i in 0..potential_energies.len() {
      let iteration = world.number_of_resets * world.max_iteration_till_reset + i;

      let potential_energy_container = potential_energies.get(i).unwrap();

      for (id, potential_energy) in potential_energy_container.iter().enumerate() {
        wtr.write_record(&[
          format!("{}", iteration),
          format!("{}", id),
          format!("{}", potential_energy),
        ])?;
      }
    }

    wtr.flush()?;

    Ok(())
  }

  fn append_positions_in_one_file_boxed_world(&self, world: &BoxedWorldDTO) -> io::Result<()> {
    let box_containers = &world.history.box_container;
    let save_dir = Path::new(&self.save_options.save_path);
    let mut wtr = csv_writer_with_header(
      &save_dir.join("positions.csv"),
      &[
        "iteration",
        "particle_index",
        "position_x",
        "position_y",
        "position_z",
      ],
    )?;

    for i in 0..box_containers.len() {
      let iteration = world.number_of_resets * world.max_iteration_till_reset + i;

      let atom_container = &box_containers.get(i).unwrap().atoms;

      for (id, atom) in atom_container.iter().enumerate() {
        wtr.write_record(&[
          format!("{}", iteration),
          format!("{}", id),
          format!("{}", atom.position.x),
          format!("{}", atom.position.y),
          format!("{}", atom.position.z),
        ])?;
      }
    }

    wtr.flush()?;

    Ok(())
  }
}
