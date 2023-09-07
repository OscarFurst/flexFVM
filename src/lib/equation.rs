use crate::cell;
use crate::cell::HasProperties;
use debug_cell::Ref;

#[derive(Clone)]
pub struct BoundaryCondition<T: HasProperties>(
    pub fn(&mut T, &[cell::Face<T>], &Vec<f64>),
    pub Vec<f64>,
);

// An Equation is a combination of a
#[derive(Clone, Copy)]
pub struct Equation<T: HasProperties>(
    pub fn(&mut T, &Vec<f64>),
    pub fn(&Ref<cell::Cell<T>>) -> Vec<f64>,
);
