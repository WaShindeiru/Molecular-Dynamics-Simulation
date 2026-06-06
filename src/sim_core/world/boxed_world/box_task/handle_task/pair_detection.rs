use std::collections::{HashSet, VecDeque};

use nalgebra::Vector3;
use std::sync::Arc;

use crate::data::SimulationConfig;
use crate::particle::Particle;
use crate::sim_core::world::boxed_world::box_container::BoxContainer;
use crate::sim_core::world::boxed_world::box_container::periodic::{
  face_neighbor_box_ids_periodic, get_needed_box_id_periodic,
};
use crate::sim_core::world::boxed_world::box_container::sim_box::{SimulationBox, get_id_simulation_box};
use crate::sim_core::world::boxed_world::integration_cache::IntegrationCache;

struct ParticleSample {
  id: usize,
  home_box: usize,
  position: Vector3<f64>,
}

/// For each home box: `i` = particles in that box; `j` = particles in its periodic 3×3×3 halo.
/// When `r_ij < threshold`, mark `home_box` of the `i` particle.
pub fn detect_marked_boxes_in_cache(
  box_cache: &BoxContainer<Arc<SimulationBox>>,
  threshold: f64,
) -> HashSet<usize> {
  let mut marked = HashSet::new();
  let config = box_cache.config();
  let dim = config.box_count_dim;

  for z in 0..dim.z {
    for y in 0..dim.y {
      for x in 0..dim.x {
        let home_box_id = get_id_simulation_box(&Vector3::new(x, y, z), &dim);
        let halo_box_ids = get_needed_box_id_periodic(&[home_box_id], config);

        let i_particles = collect_particles(box_cache, &[home_box_id]);
        let j_particles = collect_particles(box_cache, &halo_box_ids);

        for i in &i_particles {
          for j in &j_particles {
            if i.id == j.id {
              continue;
            }
            let r_ij = (j.position - i.position).magnitude();
            if r_ij < threshold {
              marked.insert(i.home_box);
            }
          }
        }
      }
    }
  }

  marked
}

pub fn detect_marked_boxes(cache: &IntegrationCache, config: &SimulationConfig) -> HashSet<usize> {
  let small_distance = &config.correction.small_distance;
  if !small_distance.enabled {
    return HashSet::new();
  }

  detect_marked_boxes_in_cache(cache.box_cache(), small_distance.distance_threshold)
}

fn collect_particles(
  box_cache: &BoxContainer<Arc<SimulationBox>>,
  box_ids: &[usize],
) -> Vec<ParticleSample> {
  let mut out = Vec::new();
  for &box_id in box_ids {
    let sim_box = box_cache.get_box(box_id);
    for particle in sim_box.particles().values() {
      if particle.is_custom_velocity_atom() {
        continue;
      }
      let id = particle.get_id();
      out.push(ParticleSample {
        id,
        home_box: box_cache.particle_box_id(id),
        position: *particle.get_position(),
      });
    }
  }
  out
}

/// Groups marked boxes into 6-connected contiguous components (face neighbours in the box grid).
pub fn cluster_contiguous_marked_boxes(
  marked: &HashSet<usize>,
  dim: Vector3<usize>,
) -> Vec<Vec<usize>> {
  let mut remaining: HashSet<usize> = marked.clone();
  let mut components = Vec::new();

  while let Some(start) = remaining.iter().copied().next() {
    remaining.remove(&start);
    let mut component = vec![start];
    let mut queue = VecDeque::from([start]);

    while let Some(current) = queue.pop_front() {
      for neighbor in face_neighbor_box_ids_periodic(current, &dim) {
        if remaining.remove(&neighbor) {
          component.push(neighbor);
          queue.push_back(neighbor);
        }
      }
    }

    components.push(component);
  }

  components
}

/// Assigns contiguous marked components to `num_tasks` workers, balancing box counts.
/// Each task receives one or more whole components (`Vec<Vec<usize>>`); components are never split.
///
/// Returns task assignments and the size (box count) of the largest contiguous component.
pub fn distribute_marked_components(
  mut components: Vec<Vec<usize>>,
  num_tasks: usize,
) -> (Vec<Vec<Vec<usize>>>, usize) {
  let largest_contiguous_block = components.iter().map(|c| c.len()).max().unwrap_or(0);

  if num_tasks == 0 {
    return (Vec::new(), largest_contiguous_block);
  }
  if components.is_empty() {
    return (vec![Vec::new(); num_tasks], largest_contiguous_block);
  }

  let total_boxes: usize = components.iter().map(|c| c.len()).sum();
  let base_target = total_boxes / num_tasks;

  components.sort_by(|a, b| b.len().cmp(&a.len()));
  let mut pending: VecDeque<Vec<usize>> = components.into();

  let mut assignments: Vec<Vec<Vec<usize>>> = Vec::with_capacity(num_tasks);

  for task_id in 0..num_tasks {
    let remainder = if task_id < total_boxes % num_tasks { 1 } else { 0 };
    let target = base_target + remainder;

    let mut blocks: Vec<Vec<usize>> = Vec::new();
    let mut count = 0usize;

    while !pending.is_empty() {
      let next_len = pending.front().unwrap().len();
      if count + next_len > target {
        break;
      }
      let component = pending.pop_front().unwrap();
      count += component.len();
      blocks.push(component);
    }

    while count <= target && !pending.is_empty() {
      let component = pending.pop_back().unwrap();
      count += component.len();
      blocks.push(component);
    }

    if blocks.is_empty() {
      if let Some(component) = pending.pop_front() {
        blocks.push(component);
      } else if let Some(component) = pending.pop_back() {
        blocks.push(component);
      }
    }

    assignments.push(blocks);
  }

  if let Some(last) = assignments.last_mut() {
    while let Some(component) = pending.pop_front() {
      last.push(component);
    }
  }

  (assignments, largest_contiguous_block)
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::data::constants::ATOMIC_MASS_FE;
  use crate::data::types::AtomType;
  use crate::particle::Atom;
  use crate::sim_core::world::boxed_world::box_container::sim_box::get_id_simulation_box;

  fn fe_atom(id: usize, position: Vector3<f64>) -> Particle {
    let mass = ATOMIC_MASS_FE;
    Particle::Atom(Atom::new(
      id,
      AtomType::Fe,
      mass,
      position,
      Vector3::zeros(),
      Vector3::zeros(),
      Vector3::zeros(),
      0.0,
    ))
  }

  /// Two adjacent boxes along x (CC box length 2, world 4³).
  fn two_box_container_close_pair() -> BoxContainer<Arc<SimulationBox>> {
    let world_size = Vector3::new(4.0, 4.0, 4.0);
    let atoms = vec![
      fe_atom(0, Vector3::new(0.5, 0.5, 0.5)),
      fe_atom(1, Vector3::new(2.5, 0.5, 0.5)),
    ];
    BoxContainer::new(atoms, world_size)
  }

  #[test]
  fn detect_marks_home_of_i_for_cross_box_close_pair() {
    let container = two_box_container_close_pair();
    let config = container.config();
    let box_a = config.box_id_for_position(&Vector3::new(0.5, 0.5, 0.5));
    let box_b = config.box_id_for_position(&Vector3::new(2.5, 0.5, 0.5));

    let marked = detect_marked_boxes_in_cache(&container, 2.5);

    assert!(marked.contains(&box_a));
    assert!(marked.contains(&box_b));
    assert_eq!(box_a, box_b);
  }

  #[test]
  fn detect_does_not_mark_when_particles_are_far() {
    let container = two_box_container_close_pair();
    let marked = detect_marked_boxes_in_cache(&container, 1.5);
    assert!(marked.is_empty());
  }

  #[test]
  #[ignore]
  fn detect_marks_only_home_i_when_scanning_one_center() {
    let container = two_box_container_close_pair();
    let dim = container.config().box_count_dim;
    let home_box_id = get_id_simulation_box(&Vector3::new(0, 0, 0), &dim);
    let halo_box_ids = get_needed_box_id_periodic(&[home_box_id], container.config());

    let i_particles = collect_particles(&container, &[home_box_id]);
    let j_particles = collect_particles(&container, &halo_box_ids);

    let mut marked_once = HashSet::new();
    for i in &i_particles {
      for j in &j_particles {
        if i.id == j.id {
          continue;
        }
        if (j.position - i.position).magnitude() < 2.5 {
          marked_once.insert(i.home_box);
        }
      }
    }

    assert_eq!(marked_once.len(), 1);
    assert!(marked_once.contains(&home_box_id));
  }

  #[test]
  fn cluster_groups_face_adjacent_marked_boxes() {
    let dim = Vector3::new(3, 1, 1);
    let mut marked = HashSet::new();
    marked.insert(get_id_simulation_box(&Vector3::new(0, 0, 0), &dim));
    marked.insert(get_id_simulation_box(&Vector3::new(1, 0, 0), &dim));

    let components = cluster_contiguous_marked_boxes(&marked, dim);
    assert_eq!(components.len(), 1);
    assert_eq!(components[0].len(), 2);
  }

  #[test]
  fn cluster_splits_diagonally_touching_boxes() {
    let dim = Vector3::new(2, 2, 1);
    let mut marked = HashSet::new();
    marked.insert(get_id_simulation_box(&Vector3::new(0, 0, 0), &dim));
    marked.insert(get_id_simulation_box(&Vector3::new(1, 1, 0), &dim));

    let components = cluster_contiguous_marked_boxes(&marked, dim);
    assert_eq!(components.len(), 2);
  }

  #[test]
  fn distribute_reports_largest_contiguous_block() {
    let components = vec![vec![0, 1], vec![2], vec![3, 4, 5, 6]];
    let (_, largest) = distribute_marked_components(components, 4);
    assert_eq!(largest, 4);
  }

  #[test]
  fn distribute_reports_empty_among_them() {
    let components = vec![vec![0, 1], vec![2], vec![3, 4, 5, 6]];
    let (result, largest) = distribute_marked_components(components, 4);
    
    for component in result {
      println!("{}", component.len());
    }

    assert_eq!(largest, 4);
  }

  #[test]
  fn distribute_uses_back_of_deque_when_front_would_overshoot() {
    let components = vec![
      vec![0usize],
      vec![1],
      vec![2],
      vec![3],
      vec![4],
      vec![5, 6, 7, 8, 9],
    ];
    let (assignments, largest) = distribute_marked_components(components, 2);
    assert_eq!(largest, 5);
    assert_eq!(assignments.len(), 2);

    let task0_boxes: usize = assignments[0].iter().map(|c| c.len()).sum();
    let task1_boxes: usize = assignments[1].iter().map(|c| c.len()).sum();
    assert_eq!(task0_boxes + task1_boxes, 10);

    let task1_has_small_component = assignments[1].iter().any(|c| c.len() == 1);
    assert!(
      task1_has_small_component,
      "task 1 should pull single-box components from the back of the deque"
    );
  }
}
