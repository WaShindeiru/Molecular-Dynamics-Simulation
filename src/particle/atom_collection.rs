use crate::particle::atom::Atom;

pub trait AtomCollection {
    fn get_atom_by_id(&self, id: u64) -> Option<&Atom>;
    fn get_all_atoms(&self) -> &Vec<Atom>;
}