use crate::persistence::dto::atom::AtomDTO;
use crate::sim_core::world::boxed_world::box_container::box_container_config::BoxContainerConfig;

pub struct BoxContainerDTO {
  pub atoms: Vec<AtomDTO>,
  pub config: BoxContainerConfig,
}
