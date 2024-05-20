use crate::computational_cell::HasProperties;
use nalgebra::DVector;

/// A BoundaryCondition is a combination of :
/// - a name
/// - a DVector
/// Basically, it only tells the equation which BC to apply, optionally including some data.
#[derive(Clone)]
pub struct BoundaryCondition(pub String, pub DVector<f64>);

// An Equation is a combination of :
// - a function that "gets" the derivated property
// - a function that "sets" the derivated property
// - a function that computes the derivative of the property
// - the dimension of the property, e.g. 1 for a scalar, 3 for a velocity vector, ...
// - a list (Vec) of possible boundary conditions
#[derive(Copy, Clone)]
pub struct Equation<T: HasProperties> {
    pub get: fn(&T, usize) -> DVector<f64>,
    pub lhs: fn(&mut T, usize, &DVector<f64>),
    pub rhs: fn(&crate::problem::Problem<T>, usize) -> DVector<f64>,
    pub dimension: usize,
}

// TODO : implement macros so that the user can simply provides a list of variable names and their dimensions
// and the macros will generate the property data structure, the get and the set procedures.
