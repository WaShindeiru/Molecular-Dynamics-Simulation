use approx::assert_abs_diff_eq;
use carbon_nanotube::sim_core::world::boxed_world::box_container::BoxContainer;
use carbon_nanotube::particle::{Particle, Atom, SafeAtomFactory};
use carbon_nanotube::data::types::{AtomType, InteractionType};
use nalgebra::Vector3;
use carbon_nanotube::sim_core::world::boxed_world::box_container::sim_box::get_id_simulation_box;

#[test]
fn test_box_container_simple_partition() {
    // Create a 10x10x10 world with CC interaction type
    // CC box size is 2.0, so we should get 5x5x5 = 125 boxes
    let size = Vector3::new(10.0, 10.0, 10.0);
    let box_type = InteractionType::CC;

    // Create empty atom list
    let atoms: Vec<Particle> = vec![];

    let container = BoxContainer::new(atoms, size, box_type, 100);

    // Verify container size
    assert_eq!(container.container_size(), &size);
    assert_eq!(container.box_type(), &InteractionType::CC);

    // With CC interaction (box_size = 2.0), we should have:
    // x: floor(10.0 / 2.0) = 5 boxes
    // y: floor(10.0 / 2.0) = 5 boxes
    // z: floor(10.0 / 2.0) = 5 boxes
    // Total: 5 * 5 * 5 = 125 boxes

    println!("Test completed: simple partition with CC interaction type");
}

#[test]
fn test_box_container_non_uniform_partition() {
    // Create a 10x6x8 world with FeFe interaction type
    // FeFe box size is 3.35, so we should get different dimensions
  let size = Vector3::new(100.0, 60.0, 80.0);
  let box_type = InteractionType::FeFe;

  let atom_factory = SafeAtomFactory::new();
  let atom0 = atom_factory.get_atom(AtomType::Fe, Vector3::new(1.0, 1.0, 1.0), Vector3::new(0., 0., 0.));
  let atom1 = atom_factory.get_atom(AtomType::Fe, Vector3::new(5.0, 3.0, 4.0), Vector3::new(0., 0., 0.));
  let atom2 = atom_factory.get_atom(AtomType::Fe, Vector3::new(9.0, 5.5, 7.5), Vector3::new(0., 0., 0.));
  let atom3 = atom_factory.get_atom(AtomType::Fe, Vector3::new(35.0, 21.37, 67.67), Vector3::new(0., 0., 0.));
  let atom4 = atom_factory.get_atom(AtomType::Fe, Vector3::new(11., 59.99, 0.), Vector3::new(0., 0., 0.));
  let atom5 = atom_factory.get_atom(AtomType::Fe, Vector3::new(100., 60., 80.), Vector3::new(0., 0., 0.));


  // Create a few test atoms
  let atoms: Vec<Particle> = vec![atom0.clone(), atom1.clone(), atom2.clone(), atom3.clone(), atom4.clone(), atom5.clone()];

  let container = BoxContainer::new(atoms, size, box_type, 100);

  // Verify container size
  assert_eq!(container.container_size(), &size);
  assert_eq!(container.box_type(), &InteractionType::FeFe);
  
  let expected_box_count_dim = Vector3::new(29, 17, 23);
  assert_eq!(container.box_count_dim(), &expected_box_count_dim);
  let expected_box_count = 29 + 17 + 23;
  assert_eq!(container.box_count(), expected_box_count);

  let expected_box_length = Vector3::new(3.44827586, 3.52941176, 3.47826087);
  assert_abs_diff_eq!(container.box_length(), &expected_box_length, epsilon =1e-6);

  // atom0
  let expected_box_coordinates: Vector3<usize> = Vector3::new(0, 0, 0);
  let expected_box_id = get_id_simulation_box(&expected_box_coordinates, &expected_box_count_dim);
  let actual_box_id = container.box_of_atom_given_index(atom0.get_id() as usize, 0).id();
  assert_eq!(actual_box_id, expected_box_id);

  // atom1
  let expected_box_coordinates: Vector3<usize> = Vector3::new(1, 0, 1);
  let expected_box_id = get_id_simulation_box(&expected_box_coordinates, &expected_box_count_dim);
  let actual_box_id = container.box_of_atom_given_index(atom1.get_id() as usize, 0).id();
  assert_eq!(actual_box_id, expected_box_id);

  // atom2
  let expected_box_coordinates: Vector3<usize> = Vector3::new(2, 1, 2);
  let expected_box_id = get_id_simulation_box(&expected_box_coordinates, &expected_box_count_dim);
  let actual_box_id = container.box_of_atom_given_index(atom2.get_id() as usize, 0).id();
  assert_eq!(actual_box_id, expected_box_id);

  // atom3
  let expected_box_coordinates: Vector3<usize> = Vector3::new(10,  6, 19);
  let expected_box_id = get_id_simulation_box(&expected_box_coordinates, &expected_box_count_dim);
  let actual_box_id = container.box_of_atom_given_index(atom3.get_id() as usize, 0).id();
  assert_eq!(actual_box_id, expected_box_id);

  // atom4
  let expected_box_coordinates: Vector3<usize> = Vector3::new(3, 16,  0);
  let expected_box_id = get_id_simulation_box(&expected_box_coordinates, &expected_box_count_dim);
  let actual_box_id = container.box_of_atom_given_index(atom4.get_id() as usize, 0).id();
  assert_eq!(actual_box_id, expected_box_id);

  // atom5
  let expected_box_coordinates: Vector3<usize> = Vector3::new(28, 16, 22);
  let expected_box_id = get_id_simulation_box(&expected_box_coordinates, &expected_box_count_dim);
  let actual_box_id = container.box_of_atom_given_index(atom5.get_id() as usize, 0).id();
  assert_eq!(actual_box_id, expected_box_id);
  
  println!("Test completed: non-uniform partition with FeFe interaction type");
}

