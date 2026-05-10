use crate::persistence::dto::world::boxed::BoxedWorldDTO;
use crate::persistence::dto::world::simple::SimpleWorldDTO;

pub mod boxed;
pub mod history;
pub mod simple;

pub enum WorldDTO {
  SimpleWorldDTO(SimpleWorldDTO),
  BoxedWorldDTO(BoxedWorldDTO),
}
