// x,y,z coordinate system, x being left-to-right, y being back-to-front and z being bottom-to-top
use crate::computational_cell::*;
use crate::equation::{BoundaryCondition, Equation};

/// The layer owns the equations for a subset of cells
#[derive(Clone)]
pub struct EquationContainer<T: HasProperties> {
    pub cell_indices: Vec<usize>,
    pub equations: Vec<Equation<T>>,
    pub boundary_conditions: Vec<BoundaryCondition>,
}

#[allow(dead_code)]
impl<T: HasProperties> EquationContainer<T> {
    pub fn n_eq(&self) -> usize {
        // returns the total number of equations that need to be solved for the layer
        self.cell_indices.len() * self.equations.len()
    }
}
