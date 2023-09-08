use flex_fvm::{
    cell,
    equation::{BoundaryCondition, Equation},
    mesh, solver,
};
use nalgebra::{dvector, DVector, Matrix3, Vector3};
use std::cell::Ref;
use std::fmt;

// This example implements the equations found at https://levelup.gitconnected.com/create-your-own-finite-volume-fluid-simulation-with-python-8f9eab0b8305
// With some clarifications from https://www.clawpack.org/riemann_book/html/Euler.html
// The equations are:
// d(rho)/dt + div(rho * u) = 0
// d(rho*u) + div(rho * u x u + p * I) = 0
// d(E) + div(u * (E + p)) = 0
// with E = rho * e + 1/2 * rho * u^2
// and p = (gamma - 1) * rho * e

fn main() {
    let initial_mesh = create_mesh();
    let solution = solver::solve(&initial_mesh, 20., 0.0001);
    println!("{:?}", solution)
}

///////////////////////
// Cell Propertires ///
///////////////////////

#[allow(non_snake_case)]
#[derive(Default, Clone, Copy)]
pub struct Properties {
    // conservative properties
    pub rho: f64,            // kg/m3
    pub rho_u: Vector3<f64>, // kg/m2/s
    pub E: f64,              // J/m3
    // derived properties
    pub u: Vector3<f64>, // m/s
    pub e: f64,          // J/kg
    pub p: f64,          // Pa
    // material properties
    pub gamma: f64, // adiabatic index
}

impl fmt::Debug for Properties {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.u[0])
    }
}

impl cell::HasProperties for Properties {
    fn update(&mut self) {
        self.u = self.rho_u / self.rho;
        self.e = self.E / self.rho - 1. / 2. * self.u.dot(&self.u);
        self.p = (self.gamma - 1.) * self.rho * self.e;
    }
}

pub fn air() -> Properties {
    let rho = 1.3;
    let u = Vector3::new(1., 0., 0.);
    let p = 100_000.;
    let gamma = 1.4;
    Properties {
        rho,
        rho_u: rho * u,
        E: p / (gamma - 1.) + 1. / 2. * rho * u.dot(&u),
        u,
        e: p / (gamma - 1.) / rho + 1. / 2. * u.dot(&u),
        p,
        gamma,
    }
}

///////////////////
// Mesh creation //
///////////////////

pub fn create_mesh() -> mesh::Mesh<Properties> {
    let bc1 = BoundaryCondition(euler_inlet, dvector![1.3, 1., 0., 0., 100_000.]);
    let bc2 = BoundaryCondition(euler_wall, dvector![0.]);
    let bc3 = BoundaryCondition(euler_fully_developed, dvector![0.]);
    let boundary_conditions: [Vec<BoundaryCondition<Properties>>; 6] = [
        vec![bc1],
        vec![bc2.clone()],
        vec![bc2.clone()],
        vec![bc3],
        vec![bc2.clone()],
        vec![bc2],
    ];
    let mut air = air();
    air.u = Vector3::new(1., 0., 0.);
    let n_grid = (1000, 30, 1);
    let mut block = mesh::Mesh::new();
    block.add_layer(
        (0., 0., 0.),
        (10., 0.3, 0.01),
        n_grid,
        boundary_conditions,
        vec![
            Equation(update_rho, euler_mass),
            Equation(update_impulse, euler_impulse),
            Equation(update_energy, euler_energy),
        ],
        air,
        air,
    );
    block
}

////////////////////////
// Equation functions //
////////////////////////

// Reminder :
// This example implements the equations found at https://levelup.gitconnected.com/create-your-own-finite-volume-fluid-simulation-with-python-8f9eab0b8305
// With some clarifications from https://www.clawpack.org/riemann_book/html/Euler.html
// The equations are:
// d(rho)/dt + div(rho * u) = 0
// d(rho*u) + div(rho * u x u + p * I) = 0
// d(E) + div(u * (E + p)) = 0
// with E = rho * e + 1/2 * rho * u^2
// and p = (gamma - 1) * rho * e

#[allow(non_snake_case)]
pub fn update_rho(p: &mut Properties, rho: &DVector<f64>) {
    p.rho = rho[0];
}

pub fn euler_mass(p: &Ref<cell::Cell<Properties>>) -> DVector<f64> {
    // mass conservation equation
    // d(rho)/dt + div(rho * u) = 0
    let mut convective_flux = 0.;
    for f in p.faces.iter() {
        let n = f.neighbour();
        let rho_avg = (p.properties.rho + n.properties.rho) / 2.;
        let u_avg: Vector3<f64> = (p.properties.u + n.properties.u) / 2.;
        convective_flux -= (rho_avg * u_avg).dot(&f.n) * f.area;
    }
    dvector![convective_flux]
}

pub fn update_impulse(p: &mut Properties, rho_u: &DVector<f64>) {
    p.rho_u[0] = rho_u[0];
    p.rho_u[1] = rho_u[1];
    p.rho_u[2] = rho_u[2];
}

pub fn euler_impulse(p: &Ref<cell::Cell<Properties>>) -> DVector<f64> {
    // momentum conservation equation
    // d(rho*u)/dt + div(rho * u x u + p * I) = 0
    let mut convective_flux = dvector![0., 0., 0.];
    for f in p.faces.iter() {
        let n = f.neighbour();
        let rho_avg = (p.properties.rho + n.properties.rho) / 2.;
        let u_avg: Vector3<f64> = (p.properties.u + n.properties.u) / 2.;
        let p_avg = (p.properties.p + n.properties.p) / 2.;
        convective_flux -= (rho_avg * u_avg * u_avg.transpose() + p_avg * Matrix3::identity())
            .transpose()
            * f.n
            * f.area
        // TODO : check if this is correct...
    }
    convective_flux
}

#[allow(non_snake_case)]
pub fn update_energy(p: &mut Properties, E: &DVector<f64>) {
    p.E = E[0];
}

#[allow(non_snake_case)]
pub fn euler_energy(p: &Ref<cell::Cell<Properties>>) -> DVector<f64> {
    // energy conservation equation
    // d(E) + div(u * (E + p)) = 0
    let mut convective_flux = 0.;
    for f in p.faces.iter() {
        let n = f.neighbour();
        let E_avg = (p.properties.E + n.properties.E) / 2.;
        let p_avg = (p.properties.p + n.properties.p) / 2.;
        let u_avg: Vector3<f64> = (p.properties.u + n.properties.u) / 2.;
        convective_flux -= (u_avg * (E_avg + p_avg)).dot(&f.n) * f.area;
    }
    dvector![convective_flux]
}

pub fn euler_inlet(p: &mut Properties, _faces: &[cell::Face<Properties>], data: &DVector<f64>) {
    // this is the inlet boundary condition, where every quantity is known
    // data = [rho, u, v, w, p]
    p.rho = data[0];
    p.rho_u[0] = data[1] * p.rho;
    p.rho_u[1] = data[2] * p.rho;
    p.rho_u[2] = data[3] * p.rho;
    p.E = data[4] / (p.gamma - 1.)
        + 1. / 2. * p.rho * (data[1].powi(2) + data[2].powi(2) + data[3].powi(2));
}

pub fn euler_wall(p: &mut Properties, faces: &[cell::Face<Properties>], _data: &DVector<f64>) {
    // this is the impermeable wall boundary condition. The wall is normal to the face
    let n = faces[0].neighbour().properties;
    let normal = faces[0].n;
    p.rho = n.rho;
    p.rho_u = n.rho_u - 2. * n.rho_u.dot(&normal) * normal; // from https://math.stackexchange.com/questions/3654725/vector-equation-of-a-line-that-is-symmetrical-to-another-line-l-with-respect-to
    p.E = n.E;
}

pub fn euler_fully_developed(
    p: &mut Properties,
    faces: &[cell::Face<Properties>],
    _data: &DVector<f64>,
) {
    // this is the fully developed boundary condition, where the flow is assumed to not change in the normal direction
    // data = [rho, u, v, w, p]
    let n = faces[0].neighbour().properties;
    p.rho = n.rho;
    p.rho_u = n.rho_u;
    p.E = n.E;
}

// pub fn euler_periodic(p: &mut Properties, faces: &[cell::Face<Properties>], data: &DVector<f64>) {
//     // this one requaires much more work...
// }
