use crate::output::world::boxed::BoxedWorldDTO;
use crate::output::world::simple::SimpleWorldDTO;

pub mod simple;
pub mod boxed;
pub mod history;

pub enum WorldDTO {
  SimpleWorldDTO(SimpleWorldDTO),
  BoxedWorldDTO(BoxedWorldDTO),
}