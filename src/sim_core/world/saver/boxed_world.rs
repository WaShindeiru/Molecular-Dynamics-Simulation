use std::{fs, io};
use std::fs::OpenOptions;
use csv::Writer;
use log::info;
use nalgebra::Vector3;
use crate::sim_core::world::integration::IntegrationAlgorithm;
use crate::output::BoxedWorldDTO;
use crate::sim_core::world::saver::PartialWorldSaver;

impl PartialWorldSaver {
  pub fn persist_boxed_world(&mut self, world: &BoxedWorldDTO) -> io::Result<()> {
    info!("Saving batch: {} ...", world.number_of_resets);
    fs::create_dir_all(&format!("./{}", self.save_options.save_path))?;

    self.append_energies_boxed_world(world)?;

    if self.save_options.save_laamps {
      self.save_laamps_boxed_world(world)?;
    }

    let (forces, potential_energies) =
      self.compute_force_and_potential_boxed_world(world);

    if self.save_options.save_verbose {
      self.append_forces_boxed_world(world, &forces)?;
      self.append_potential_energies_boxed_world(world, &potential_energies)?;
      self.append_positions_in_one_file_boxed_world(world)?;
    }

    info!("Batch {} saved", world.number_of_resets);

    Ok(())
  }

  pub fn append_energies_boxed_world(&mut self, world: &BoxedWorldDTO) -> io::Result<()> {
    let atoms = &world.box_container.atoms;
    let num_of_iterations = world.num_of_world_iterations;

    let mut kinetic_energy: Vec<f64> = Vec::with_capacity(num_of_iterations);
    let mut potential_energy: Vec<f64> = Vec::with_capacity(num_of_iterations);
    let mut total_energy: Vec<f64> = Vec::with_capacity(num_of_iterations);
    let mut thermostat_work: Vec<f64> = Vec::with_capacity(num_of_iterations);

    for (i, atom_container) in atoms.iter().enumerate() {
      let mut kinetic_energy_i = 0.;
      let mut thermostat_work_i = 0.;
      let mut potential_energy_i = 0.;

      for atom_dto in atom_container.iter() {
        kinetic_energy_i += atom_dto.kinetic_energy;
        thermostat_work_i += atom_dto.thermostat_work;
        potential_energy_i += atom_dto.potential_energy;
      }

      self.thermostat_work_total += thermostat_work_i;

      kinetic_energy.push(kinetic_energy_i);
      potential_energy.push(potential_energy_i);
      let potential_energy_i = *potential_energy.get(i).unwrap();
      total_energy.push(kinetic_energy_i + potential_energy_i);
      thermostat_work.push(self.thermostat_work_total);
    }

    let file = OpenOptions::new()
      .create(true)
      .append(true)
      .open(&format!("./{}/energy.csv", self.save_options.save_path))?;

    let mut wtr = Writer::from_writer(file);

    assert!(kinetic_energy.len() == potential_energy.len() && kinetic_energy.len() == total_energy.len());

    for i in 0..kinetic_energy.len() {
      let iteration = world.number_of_resets * world.max_iteration_till_reset + i;

      match world.integration_algorithm {
        IntegrationAlgorithm::NoseHooverVerlet => {
          wtr.write_record(&[
            format!("{}", iteration),
            format!("{}", kinetic_energy.get(i).unwrap()),
            format!("{}", potential_energy.get(i).unwrap()),
            format!("{}", total_energy.get(i).unwrap()),
            format!("{}", thermostat_work.get(i).unwrap()),
            format!("{}", world.box_container.thermostat_epsilon.get(i).unwrap()),
          ])?;
        }
        _ => {
          wtr.write_record(&[
            format!("{}", iteration),
            format!("{}", kinetic_energy.get(i).unwrap()),
            format!("{}", potential_energy.get(i).unwrap()),
            format!("{}", total_energy.get(i).unwrap()),
          ])?;
        }
      }
    }

    wtr.flush()?;

    Ok(())
  }

  fn save_laamps_boxed_world(&mut self, world: &BoxedWorldDTO) -> io::Result<()> {
    let atoms = &world.box_container.atoms;
    let mut files_saved = 0;

    for i in 0..atoms.len() {
      if self.frame_iteration_count_current_iteration % world.frame_iteration_count == 0 {
        let mut result_string = String::new();
        result_string.push_str(&"ITEM: TIMESTEP\n".to_string());
        result_string.push_str(&format!("{}\n", i));

        result_string.push_str(&"ITEM: NUMBER OF ATOMS\n".to_string());
        let atom_container = atoms.get(i).unwrap();
        let num_of_atoms = atom_container.len();
        result_string.push_str(&format!("{}\n", num_of_atoms));

        result_string.push_str(&"ITEM: BOX BOUNDS pp pp pp\n".to_string());
        result_string.push_str(&format!("0.0 {}\n", world.size.x));
        result_string.push_str(&format!("0.0 {}\n", world.size.y));
        result_string.push_str(&format!("0.0 {}\n", world.size.z));

        result_string.push_str(&"ITEM: ATOMS id type x y z\n".to_string());

        for atom_dto in atom_container.iter() {
          result_string.push_str(&format!("{} {} {} {} {}\n", atom_dto.id, atom_dto.atom_type, atom_dto.x, atom_dto.y, atom_dto.z));
        }

        let file_number = i + world.number_of_resets * world.max_iteration_till_reset;
        fs::write(&format!("./{}/output_{}.dump", self.save_options.save_path, file_number), result_string)?;
        files_saved += 1;
      }

      self.frame_iteration_count_current_iteration += 1
    }

    info!("Saved {} laamps output files.", {files_saved});
    Ok(())
  }

  fn compute_force_and_potential_boxed_world(&self, world: &BoxedWorldDTO) -> (Vec<Vec<Vector3<f64>>>, Vec<Vec<f64>>) {
    let atoms = &world.box_container.atoms;
    let mut forces: Vec<Vec<Vector3<f64>>> = Vec::new();
    let mut potential_energies: Vec<Vec<f64>> = Vec::new();

    for (i, atom_container) in atoms.iter().enumerate() {
      let mut current_forces: Vec<Vector3<f64>> = vec![Vector3::new(0., 0., 0.); world.num_of_atoms];
      let mut current_potential_energies: Vec<f64> = vec![0.; world.num_of_atoms];

      for atom_dto in atom_container.iter() {
        let force_i = current_forces.get_mut(atom_dto.id as usize).unwrap();
        force_i.x = atom_dto.force_x;
        force_i.y = atom_dto.force_y;
        force_i.z = atom_dto.force_z;

        *current_potential_energies.get_mut(atom_dto.id as usize).unwrap() = atom_dto.potential_energy;
      }

      forces.push(current_forces);
      potential_energies.push(current_potential_energies);
    }

    (forces, potential_energies)
  }

  // TODO: this method is exactly the same for simple world, refactor this please
  fn append_forces_boxed_world(&self, world: &BoxedWorldDTO, forces: &Vec<Vec<Vector3<f64>>>) -> io::Result<()> {
    let file = OpenOptions::new()
      .create(true)
      .append(true)
      .open(&format!("./{}/forces.csv", self.save_options.save_path))?;

    let mut wtr = Writer::from_writer(file);

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
  fn append_potential_energies_boxed_world(&self, world: &BoxedWorldDTO, potential_energies: &Vec<Vec<f64>>) -> io::Result<()> {
    let file = OpenOptions::new()
      .create(true)
      .append(true)
      .open(&format!("./{}/potential_energies.csv", self.save_options.save_path))?;

    let mut wtr = Writer::from_writer(file);

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
    let atoms = &world.box_container.atoms;
    let file = OpenOptions::new()
      .create(true)
      .append(true)
      .open(&format!("./{}/positions.csv", self.save_options.save_path))?;

    let mut wtr = Writer::from_writer(file);

    for i in 0..atoms.len() {
      let iteration = world.number_of_resets * world.max_iteration_till_reset + i;

      let atom_container = atoms.get(i).unwrap();

      for (id, atom) in atom_container.iter().enumerate() {
        wtr.write_record(&[
          format!("{}", iteration),
          format!("{}", id),
          format!("{}", atom.x),
          format!("{}", atom.y),
          format!("{}", atom.z),
        ])?;
      }
    }

    wtr.flush()?;

    Ok(())
  }
}