use crate::cell;
use debug_cell::Ref;

#[derive(Clone, Copy)]
pub struct BoundaryCondition(pub fn(&mut cell::Properties, &[cell::Face], f64), pub f64);

#[derive(Clone, Copy)]
pub struct Equation(
    pub fn(&mut cell::Properties, f64),
    pub fn(&Ref<cell::Cell>) -> f64,
);

///////////////////////////////
/// Heat diffusion equation ///
///////////////////////////////
#[allow(non_snake_case)]
pub fn update_T(p: &mut cell::Properties, T: f64) {
    p.T = T;
}

pub fn heat_convection_diffusion(p: &Ref<cell::Cell>) -> f64 {
    let mut diffusive_flux = 0.;
    let mut convective_flux = 0.;
    for f in p.faces.iter() {
        let n = f.neighbour();
        diffusive_flux -= p.properties.D * (p.properties.T - n.properties.T) / f.distance * f.area;
        convective_flux -=
            ((p.properties.h() + n.properties.h()) * (p.properties.u + n.properties.u) / 4.)
                .dot(&f.n)
                * f.area;
    }
    (diffusive_flux + convective_flux) / p.volume / p.properties.rho / p.properties.cp
}

pub fn heat_diffusion(p: &Ref<cell::Cell>) -> f64 {
    let mut flux = 0.;
    for f in p.faces.iter() {
        let n = f.neighbour();
        flux -= p.properties.lambda * (p.properties.T - n.properties.T) / f.distance * f.area;
    }
    flux / p.volume / p.properties.rho / p.properties.cp
}

#[allow(non_snake_case, unused_variables)]
pub fn heat_dirichlet(p: &mut cell::Properties, faces: &[cell::Face], T: f64) {
    p.T = T;
}

#[allow(non_snake_case)]
pub fn heat_neumann(p: &mut cell::Properties, faces: &[cell::Face], dTdx: f64) {
    let n = faces[0].neighbour();
    let mut boundary_value = n.properties.T;
    boundary_value += dTdx * faces[0].distance;
    p.T = boundary_value;
}
