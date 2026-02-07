use std::fs::OpenOptions;
use std::{fs, io};
use csv::Writer;
use log::info;
use nalgebra::Vector3;
use crate::output::{AtomDTO, WorldDTO};
use super::World;

impl World {

  pub fn to_transfer_struct(&self) -> WorldDTO {
    let mut all_atoms_dto: Vec<Vec<AtomDTO>> = Vec::with_capacity(self.atoms.len());
    let mut potential_energies: Vec<f64> = Vec::with_capacity(self.atoms.len());

    let lower_index: usize;

    if self.number_of_resets > 0 {
      lower_index = 1;
    } else {
      lower_index = 0;
    }

    for iteration in lower_index..self.atoms.len() {
      let atom_container = self.atoms.get(iteration).unwrap();
      let mut atoms_dto: Vec<AtomDTO> = Vec::with_capacity(atom_container.len());

      for atom in atom_container.get_atoms().iter() {
        atoms_dto.push(atom.to_transfer_struct());
      }

      all_atoms_dto.push(atoms_dto);
      potential_energies.push(atom_container.get_potential_energy());
    }

    WorldDTO {
      num_of_atoms: self.atom_count,
      atoms: all_atoms_dto,
      potential_energy: potential_energies,
      thermostat_epsilon: self.thermostat_epsilon.clone(),
      box_x: self.size.x,
      box_y: self.size.y,
      box_z: self.size.z,
    }
  }

  fn save_energies(&self, world: &WorldDTO, atoms: &Vec<Vec<AtomDTO>>, use_thermostat: bool) -> io::Result<(Vec<Vec<Vector3<f64>>>, Vec<Vec<f64>>)> {
    let num_of_iterations = self.atoms.len() - 1;

    let mut kinetic_energy: Vec<f64> = Vec::with_capacity(num_of_iterations);
    // let mut potential_energy: Vec<f64> = Vec::with_capacity(engine_dto.num_of_iterations);
    let potential_energy = &world.potential_energy;
    let mut total_energy: Vec<f64> = Vec::with_capacity(num_of_iterations);
    let mut thermostat_work: Vec<f64> = Vec::with_capacity(num_of_iterations);

    let mut forces: Vec<Vec<Vector3<f64>>> = Vec::new();
    let mut potential_energies: Vec<Vec<f64>> = Vec::new();

    let mut thermostat_work_sum = 0.;

    for (i, atom_container) in atoms.iter().enumerate() {
      let mut current_forces: Vec<Vector3<f64>> = vec![Vector3::new(0., 0., 0.); world.num_of_atoms];
      let mut current_potential_energies: Vec<f64> = vec![0.; world.num_of_atoms];
      let mut kinetic_energy_i = 0.;
      let mut thermostat_work_i = 0.;
      // let mut potential_energy_i = 0.;

      for atom_dto in atom_container.iter() {
        kinetic_energy_i += atom_dto.kinetic_energy;
        thermostat_work_i += atom_dto.thermostat_work;
        // potential_energy_i += atom_dto.potential_energy;

        let force_i = current_forces.get_mut(atom_dto.id as usize).unwrap();
        force_i.x = atom_dto.force_x;
        force_i.y = atom_dto.force_y;
        force_i.z = atom_dto.force_z;

        *current_potential_energies.get_mut(atom_dto.id as usize).unwrap() = atom_dto.potential_energy;
      }

      thermostat_work_sum += thermostat_work_i;

      kinetic_energy.push(kinetic_energy_i);
      // potential_energy.push(potential_energy_i);
      let potential_energy_i = *potential_energy.get(i).unwrap();
      total_energy.push(kinetic_energy_i + potential_energy_i);
      thermostat_work.push(thermostat_work_sum);
      forces.push(current_forces);
      potential_energies.push(current_potential_energies);
    }

    let mut wtr = Writer::from_path(&format!("./{}/energy.csv", self.save_path))?;

    assert!(kinetic_energy.len() == potential_energy.len() && kinetic_energy.len() == total_energy.len());

    for i in 0..kinetic_energy.len() {
      if use_thermostat {
        wtr.write_record(&[
          format!("{}", i+1),
          format!("{}", kinetic_energy.get(i).unwrap()),
          format!("{}", potential_energy.get(i).unwrap()),
          format!("{}", total_energy.get(i).unwrap()),
          format!("{}", thermostat_work.get(i).unwrap()),
          format!("{}", world.thermostat_epsilon.get(i).unwrap()),
        ])?;
      } else {
        wtr.write_record(&[
          format!("{}", i+1),
          format!("{}", kinetic_energy.get(i).unwrap()),
          format!("{}", potential_energy.get(i).unwrap()),
          format!("{}", total_energy.get(i).unwrap()),
        ])?;
      }

    }

    wtr.flush()?;

    Ok((forces, potential_energies))
  }

  fn save_laamps(&mut self, atoms: &Vec<Vec<AtomDTO>>) -> io::Result<()> {
    let mut files_saved = 0;

    for i in 0..atoms.len() {
      if self.frame_iteration_count_current_iteration % self.frame_iteration_count == 0 {
        let mut result_string = String::new();
        result_string.push_str(&"ITEM: TIMESTEP\n".to_string());
        result_string.push_str(&format!("{}\n", i+1));

        result_string.push_str(&"ITEM: NUMBER OF ATOMS\n".to_string());
        let atom_container = atoms.get(i).unwrap();
        let num_of_atoms = atom_container.len();
        result_string.push_str(&format!("{}\n", num_of_atoms));

        result_string.push_str(&"ITEM: BOX BOUNDS pp pp pp\n".to_string());
        result_string.push_str(&format!("0.0 {}\n", self.size.x));
        result_string.push_str(&format!("0.0 {}\n", self.size.y));
        result_string.push_str(&format!("0.0 {}\n", self.size.z));

        result_string.push_str(&"ITEM: ATOMS id type x y z\n".to_string());

        for atom_dto in atom_container.iter() {
          result_string.push_str(&format!("{} {} {} {} {}\n", atom_dto.id, atom_dto.atom_type, atom_dto.x, atom_dto.y, atom_dto.z));
        }

        let file_number = i + self.number_of_resets * self.max_iteration_till_reset;
        fs::write(&format!("./{}/output_{}.dump", self.save_path, file_number), result_string)?;
        files_saved += 1;
      }
    }

    info!("Saved {} laamps output files.", {files_saved});
    Ok(())
  }

  fn append_forces(&self, forces: &Vec<Vec<Vector3<f64>>>) -> io::Result<()> {
    let file = OpenOptions::new()
      .create(true)
      .append(true)
      .open(&format!("./{}/forces.csv", self.save_path))?;

    let mut wtr = Writer::from_writer(file);

    for i in 0..forces.len() {
      let force_map = forces.get(i).unwrap();

      for (id, force_vector) in force_map.iter().enumerate() {
        wtr.write_record(&[
          format!("{}", i+1),
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

  fn append_potential_energies(&self, potential_energies: &Vec<Vec<f64>>) -> io::Result<()> {
    let file = OpenOptions::new()
      .create(true)
      .append(true)
      .open(&format!("./{}/potential_energies.csv", self.save_path))?;

    let mut wtr = Writer::from_writer(file);

    for i in 0..potential_energies.len() {
      let potential_energy_container = potential_energies.get(i).unwrap();

      for (id, potential_energy) in potential_energy_container.iter().enumerate() {
        wtr.write_record(&[
          format!("{}", i+1),
          format!("{}", id),
          format!("{}", potential_energy),
        ])?;
      }
    }

    wtr.flush()?;

    Ok(())
  }

  fn append_positions_in_one_file(&self, atoms: &Vec<Vec<AtomDTO>>) -> io::Result<()> {
    let file = OpenOptions::new()
      .create(true)
      .append(true)
      .open(&format!("./{}/positions.csv", self.save_path))?;

    let mut wtr = Writer::from_writer(file);

    for i in 0..atoms.len() {
      let atom_container = atoms.get(i).unwrap();

      for (id, atom) in atom_container.iter().enumerate() {
        wtr.write_record(&[
          format!("{}", i+1),
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

  pub fn save(&mut self, use_thermostat: bool) -> io::Result<()> {
    info!("Saving batch: {} ...", self.number_of_resets);
    fs::create_dir_all(&format!("./{}", self.save_path))?;

    let world = self.to_transfer_struct();
    let atoms = &world.atoms;

    if self.save_laamps {
      self.save_laamps(&atoms)?;
    }

    let (forces, potential_energies) = self.save_energies(&world, &atoms, use_thermostat)?;

    if self.save_verbose {
      self.append_forces(&forces)?;
      self.append_potential_energies(&potential_energies)?;
      self.append_positions_in_one_file(&atoms)?;
    }

    info!("Batch {} saved", self.number_of_resets);
    Ok(())
  }
}