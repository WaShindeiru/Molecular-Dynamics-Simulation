pub enum IntegrationAlgorithm {
  SemiImplicitEuler,
  VelocityVerlet,
  NoseHooverVerlet {
    desired_temperature: f64,
    q_effective_mass: f64
  },
}

