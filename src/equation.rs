use crate::cell::Cell;
use debug_cell::{Ref, RefMut};
use std::borrow::BorrowMut;

///////////////////////////////
/// Heat diffusion equation ///
///////////////////////////////
pub fn heat_diffusion(p: &Ref<Cell>) -> f64 {
    let mut flux = 0.;
    for f in p.faces.iter() {
        let n = f.neighbor.as_ref().unwrap().cell.borrow();
        flux -= p.lambda * (p.T - n.T) / f.neighbor.as_ref().unwrap().distance * f.area;
    }
    flux / p.volume / p.rho / p.cp
}

pub fn heat_diffusion_dirichlet(p: &mut RefMut<Cell>, T: f64) {
    p.T = T;
}

pub fn heat_diffusion_neumann(p: &mut RefMut<Cell>, dTdx: f64) {
    let mut boundary_value = p.faces[0] // should not be mutably borrowed but somehow gives me an error if I don't...
        .neighbor
        .as_ref()
        .unwrap()
        .cell
        .borrow()
        .T;
    boundary_value += dTdx * p.borrow_mut().faces[0].neighbor.as_ref().unwrap().distance;
    p.borrow_mut().T = boundary_value;
}
