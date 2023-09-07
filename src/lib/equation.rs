use crate::cell;
use crate::cell::HasProperties;
use debug_cell::Ref;
use nalgebra::DVector;

#[derive(Clone)]
pub struct BoundaryCondition<T: HasProperties>(
    pub fn(&mut T, &[cell::Face<T>], &DVector<f64>),
    pub DVector<f64>,
);

// An Equation is a combination of a
#[derive(Clone, Copy)]
pub struct Equation<T: HasProperties>(
    pub fn(&mut T, &DVector<f64>),
    pub fn(&Ref<cell::Cell<T>>) -> DVector<f64>,
);
