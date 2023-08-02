use crate::cell;
use debug_cell::Ref;

#[derive(Clone, Copy)]
pub struct BoundaryCondition(pub fn(&mut cell::Properties, &[cell::Face], f64), pub f64);

#[derive(Clone, Copy)]
pub struct Equation(
    pub fn(&mut cell::Properties, f64),
    pub fn(&Ref<cell::Cell>) -> f64,
);
