use crate::cell;
use crate::cell::HasProperties;
use debug_cell::Ref;

#[derive(Clone, Copy)]
pub struct BoundaryCondition<T: HasProperties>(pub fn(&mut T, &[cell::Face<T>], f64), pub f64);

#[derive(Clone, Copy)]
pub struct Equation<T: HasProperties>(pub fn(&mut T, f64), pub fn(&Ref<cell::Cell<T>>) -> f64);
