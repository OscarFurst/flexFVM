use crate::cell;
use debug_cell::Ref;

#[derive(Clone, Copy)]
pub struct BoundaryCondition(pub fn(&mut cell::Properties, &[cell::Face], f64), pub f64);

///////////////////////////////
/// Heat diffusion equation ///
///////////////////////////////

pub fn heat_convection_diffusion(p: &Ref<cell::Cell>) -> f64 {
    let mut diffusive_flux = 0.;
    let mut convective_flux = 0.;
    for f in p.faces.iter() {
        let n = f.neighbor.as_ref().unwrap().cell.borrow();
        diffusive_flux -= p.properties.D * (p.properties.T - n.properties.T)
            / f.neighbor.as_ref().unwrap().distance
            * f.area;
        convective_flux -= ((p.properties.T + n.properties.T) * (p.properties.u + n.properties.u)
            / 4.)
            .dot(&f.neighbor.as_ref().unwrap().n)
            * f.area;
    }
    (diffusive_flux + convective_flux) / p.volume / p.properties.rho / p.properties.cp
}

pub fn heat_diffusion(p: &Ref<cell::Cell>) -> f64 {
    let mut flux = 0.;
    for f in p.faces.iter() {
        let n = f.neighbor.as_ref().unwrap().cell.borrow();
        flux -= p.properties.lambda * (p.properties.T - n.properties.T)
            / f.neighbor.as_ref().unwrap().distance
            * f.area;
    }
    flux / p.volume / p.properties.rho / p.properties.cp
}

#[allow(non_snake_case, unused_variables)]
pub fn heat_dirichlet(p: &mut cell::Properties, faces: &[cell::Face], T: f64) {
    p.T = T;
}

#[allow(non_snake_case)]
pub fn heat_neumann(p: &mut cell::Properties, faces: &[cell::Face], dTdx: f64) {
    let onlyneighbor = faces[0].neighbor.as_ref().unwrap();
    let mut boundary_value = onlyneighbor.cell.borrow().properties.T;
    boundary_value += dTdx * onlyneighbor.distance;
    p.T = boundary_value;
}
