use rand::distr::{Distribution, Uniform};
use rand_distr::Normal;

use nalgebra::Vector3;

use crate::data::ParticleConfig;
use crate::data::types::AtomType::{C, Fe};
use crate::particle::{Particle, SafeAtomFactory};
use crate::simulations::generators::generate::Generator;

pub struct DenseGenerator {
  potential_gravity_max: f64,
  world_size: Vector3<f64>,
  particle_distance: f64,
  offset: Vector3<f64>,
}

impl DenseGenerator {
  pub fn new(
    potential_gravity_max: f64,
    world_size: Vector3<f64>,
    particle_distance: f64,
    offset: Vector3<f64>,
  ) -> Self {
    Self {
      potential_gravity_max,
      world_size,
      particle_distance,
      offset,
    }
  }
}

impl Generator for DenseGenerator {
  fn generate(&self) -> ParticleConfig {
    let atom_factory = SafeAtomFactory::new(self.potential_gravity_max, self.world_size.z);

    let mut atoms: Vec<Particle> = Vec::new();
    let mut count = 0;
    let mut fe_count = 0;
    let mut c_count = 0;

    let normal = Normal::new(0.0, 1.5).unwrap();
    let mut rng = rand::rng();

    let range = Uniform::new(0., 1.).unwrap();
    let mut rng_2 = rand::rng();

    let velocity_normal = Normal::new(0.0, 1e-2).unwrap();
    let mut rng_3 = rand::rng();

    let mut z = self.offset.z;
    while z <= self.world_size.z - self.offset.z {
      let mut y = self.offset.y;

      while y <= self.world_size.y - self.offset.y {
        let mut x = self.offset.x;

        while x <= self.world_size.x - self.offset.x {
          let x_ = loop {
            let x_candidate = x + normal.sample(&mut rng);
            if x_candidate >= 0.0 && x_candidate <= self.world_size.x {
              break x_candidate;
            }
          };

          let y_ = loop {
            let y_candidate = y + normal.sample(&mut rng);
            if y_candidate >= 0.0 && y_candidate <= self.world_size.y {
              break y_candidate;
            }
          };

          let z_ = loop {
            let z_candidate = z + normal.sample(&mut rng);
            if z_candidate >= 0.0 && z_candidate <= self.world_size.z {
              break z_candidate;
            }
          };

          let position = Vector3::new(x_, y_, z_);

          let atom_type = match range.sample(&mut rng_2) {
            x if x < 0.75 => {
              fe_count += 1;
              Fe
            }
            x if x >= 0.75 => {
              c_count += 1;
              C
            }
            _ => panic!("what's that?"),
          };

          let v_x = velocity_normal.sample(&mut rng_3);
          let v_y = velocity_normal.sample(&mut rng_3);
          let v_z = velocity_normal.sample(&mut rng_3);
          let particle = atom_factory.get_atom(atom_type, position, Vector3::new(v_x, v_y, v_z));
          atoms.push(particle);
          count += 1;

          x = x + self.particle_distance;
        }

        y = y + self.particle_distance;
      }

      z = z + self.particle_distance;
    }

    ParticleConfig::new(atoms)
  }
}
