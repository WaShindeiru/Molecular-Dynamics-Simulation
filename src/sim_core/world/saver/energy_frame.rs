use std::io;

use csv::Writer;

/// How to turn a frame window (`frame_iteration_count` steps) into one `energy.csv` row.
///
/// LAMMPS dumps always stay snapshots; this only applies to energy/temperature series.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum FrameReduce {
  /// Write the first iteration of each window (historical behavior).
  #[default]
  Snapshot,
  /// Write the arithmetic mean of every iteration in the window.
  Mean,
}

impl FrameReduce {
  pub fn is_snapshot(&self) -> bool {
    matches!(self, FrameReduce::Snapshot)
  }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(super) struct EnergyCsvLayout {
  pub nose_hoover: bool,
  pub nanotube: bool,
}

#[derive(Clone, Copy, Debug, Default)]
pub(super) struct EnergyRow {
  pub kinetic_energy_atom: f64,
  pub kinetic_energy_other: f64,
  pub potential_energy: f64,
  pub potential_gravity_energy: f64,
  pub total_energy: f64,
  pub p_control_energy_total: f64,
  pub thermostat_work_total: f64,
  pub thermostat_epsilon: f64,
  pub temperature: f64,
  pub nanotube_thermostat_epsilon: f64,
  pub nanotube_temperature: f64,
}

impl EnergyRow {
  fn add_assign(&mut self, other: &EnergyRow) {
    self.kinetic_energy_atom += other.kinetic_energy_atom;
    self.kinetic_energy_other += other.kinetic_energy_other;
    self.potential_energy += other.potential_energy;
    self.potential_gravity_energy += other.potential_gravity_energy;
    self.total_energy += other.total_energy;
    self.p_control_energy_total += other.p_control_energy_total;
    self.thermostat_work_total += other.thermostat_work_total;
    self.thermostat_epsilon += other.thermostat_epsilon;
    self.temperature += other.temperature;
    self.nanotube_thermostat_epsilon += other.nanotube_thermostat_epsilon;
    self.nanotube_temperature += other.nanotube_temperature;
  }

  fn mean(self, count: usize) -> EnergyRow {
    let n = count as f64;
    EnergyRow {
      kinetic_energy_atom: self.kinetic_energy_atom / n,
      kinetic_energy_other: self.kinetic_energy_other / n,
      potential_energy: self.potential_energy / n,
      potential_gravity_energy: self.potential_gravity_energy / n,
      total_energy: self.total_energy / n,
      p_control_energy_total: self.p_control_energy_total / n,
      thermostat_work_total: self.thermostat_work_total / n,
      thermostat_epsilon: self.thermostat_epsilon / n,
      temperature: self.temperature / n,
      nanotube_thermostat_epsilon: self.nanotube_thermostat_epsilon / n,
      nanotube_temperature: self.nanotube_temperature / n,
    }
  }
}

#[derive(Default)]
pub(super) struct EnergyFrameAccumulator {
  count: usize,
  start_iteration: usize,
  sums: EnergyRow,
  layout: Option<EnergyCsvLayout>,
}

impl EnergyFrameAccumulator {
  pub(super) fn set_layout(&mut self, layout: EnergyCsvLayout) {
    debug_assert!(
      self.layout.is_none() || self.layout == Some(layout),
      "energy.csv layout must stay the same for the whole run"
    );
    self.layout = Some(layout);
  }

  pub(super) fn layout(&self) -> Option<EnergyCsvLayout> {
    self.layout
  }

  /// Accumulate `row` and, when `window` samples are in, return `(window_start_iteration, mean)`.
  pub(super) fn push(
    &mut self,
    iteration: usize,
    row: EnergyRow,
    window: usize,
  ) -> Option<(usize, EnergyRow)> {
    if self.count == 0 {
      self.start_iteration = iteration;
    }
    self.sums.add_assign(&row);
    self.count += 1;
    if self.count >= window.max(1) {
      self.take_mean()
    } else {
      None
    }
  }

  pub(super) fn take_partial(&mut self) -> Option<(usize, EnergyRow)> {
    if self.count == 0 {
      None
    } else {
      self.take_mean()
    }
  }

  fn take_mean(&mut self) -> Option<(usize, EnergyRow)> {
    if self.count == 0 {
      return None;
    }
    let mean = self.sums.mean(self.count);
    let start = self.start_iteration;
    self.count = 0;
    self.sums = EnergyRow::default();
    Some((start, mean))
  }
}

pub(super) fn energy_csv_header(layout: EnergyCsvLayout) -> Vec<&'static str> {
  let mut header: Vec<&str> = if layout.nose_hoover {
    vec![
      "iteration",
      "kinetic_energy_atom",
      "kinetic_energy_other",
      "potential_energy",
      "potential_gravity_energy",
      "total_energy",
      "p_control_energy_total",
      "thermostat_work_total",
      "thermostat_epsilon",
      "temperature",
    ]
  } else {
    vec![
      "iteration",
      "kinetic_energy_atom",
      "kinetic_energy_other",
      "potential_energy",
      "potential_gravity_energy",
      "total_energy",
      "p_control_energy_total",
      "temperature",
    ]
  };
  if layout.nanotube {
    header.push("nanotube_thermostat_epsilon");
    header.push("nanotube_temperature");
  }
  header
}

pub(super) fn write_energy_record<W: io::Write>(
  wtr: &mut Writer<W>,
  iteration: usize,
  row: &EnergyRow,
  layout: EnergyCsvLayout,
) -> csv::Result<()> {
  match (layout.nose_hoover, layout.nanotube) {
    (true, true) => wtr.write_record(&[
      format!("{}", iteration),
      format!("{}", row.kinetic_energy_atom),
      format!("{}", row.kinetic_energy_other),
      format!("{}", row.potential_energy),
      format!("{}", row.potential_gravity_energy),
      format!("{}", row.total_energy),
      format!("{}", row.p_control_energy_total),
      format!("{}", row.thermostat_work_total),
      format!("{}", row.thermostat_epsilon),
      format!("{}", row.temperature),
      format!("{}", row.nanotube_thermostat_epsilon),
      format!("{}", row.nanotube_temperature),
    ]),
    (true, false) => wtr.write_record(&[
      format!("{}", iteration),
      format!("{}", row.kinetic_energy_atom),
      format!("{}", row.kinetic_energy_other),
      format!("{}", row.potential_energy),
      format!("{}", row.potential_gravity_energy),
      format!("{}", row.total_energy),
      format!("{}", row.p_control_energy_total),
      format!("{}", row.thermostat_work_total),
      format!("{}", row.thermostat_epsilon),
      format!("{}", row.temperature),
    ]),
    (false, true) => wtr.write_record(&[
      format!("{}", iteration),
      format!("{}", row.kinetic_energy_atom),
      format!("{}", row.kinetic_energy_other),
      format!("{}", row.potential_energy),
      format!("{}", row.potential_gravity_energy),
      format!("{}", row.total_energy),
      format!("{}", row.p_control_energy_total),
      format!("{}", row.temperature),
      format!("{}", row.nanotube_thermostat_epsilon),
      format!("{}", row.nanotube_temperature),
    ]),
    (false, false) => wtr.write_record(&[
      format!("{}", iteration),
      format!("{}", row.kinetic_energy_atom),
      format!("{}", row.kinetic_energy_other),
      format!("{}", row.potential_energy),
      format!("{}", row.potential_gravity_energy),
      format!("{}", row.total_energy),
      format!("{}", row.p_control_energy_total),
      format!("{}", row.temperature),
    ]),
  }
}

#[cfg(test)]
mod tests {
  use super::*;

  fn row(temperature: f64) -> EnergyRow {
    EnergyRow {
      temperature,
      kinetic_energy_atom: temperature,
      ..EnergyRow::default()
    }
  }

  #[test]
  fn mean_emits_after_full_window_with_start_iteration() {
    let mut acc = EnergyFrameAccumulator::default();
    assert!(acc.push(10, row(1.0), 3).is_none());
    assert!(acc.push(11, row(2.0), 3).is_none());
    let (iteration, mean) = acc.push(12, row(3.0), 3).unwrap();
    assert_eq!(iteration, 10);
    assert!((mean.temperature - 2.0).abs() < 1e-12);
    assert!(acc.take_partial().is_none());
  }

  #[test]
  fn partial_window_flushes_remaining_mean() {
    let mut acc = EnergyFrameAccumulator::default();
    assert!(acc.push(0, row(4.0), 4).is_none());
    assert!(acc.push(1, row(6.0), 4).is_none());
    let (iteration, mean) = acc.take_partial().unwrap();
    assert_eq!(iteration, 0);
    assert!((mean.temperature - 5.0).abs() < 1e-12);
  }
}
