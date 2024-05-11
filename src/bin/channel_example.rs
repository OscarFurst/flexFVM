use flex_fvm::problem;
use flex_fvm::{
    computational_cell,
    equation::{BoundaryCondition, Equation},
    generate_equation_structure, solver,
};
use nalgebra::{dvector, DVector, Matrix3, Point3, Vector3};
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
    let output_args = solver::SolverOutputArgs {
        frequency: 0.02,
        filename: "channel_solution".to_string(),
    };
    let _solution = solver::solve(
        &initial_mesh,
        2.,
        0.001,
        solver::runge_kutta_4::<Properties>,
        output_args,
    );
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

impl computational_cell::HasProperties for Properties {
    fn update(&mut self) {
        self.u = self.rho_u / self.rho;
        self.e = self.E / self.rho - 1. / 2. * self.u.dot(&self.u);
        self.p = (self.gamma - 1.) * self.rho * self.e;
    }

    fn export(&self) -> Vec<(String, f64)> {
        vec![
            ("ux".to_string(), self.u[0]),
            ("uy".to_string(), self.u[1]),
            ("uz".to_string(), self.u[2]),
            ("p".to_string(), self.p),
            ("rho".to_string(), self.rho),
        ]
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
        E: p / (gamma - 1.) + 0.5 * rho * u.dot(&u),
        u,
        e: p / (gamma - 1.) / rho + 0.5 * u.dot(&u),
        p,
        gamma,
    }
}

fn properties_generator(position: &Point3<f64>) -> Properties {
    let w0 = 0.1;
    let sigma = 0.05f64.powi(2);
    let rho = 1. + (position.y - 0.5).abs().clamp(0., 0.25);
    let u = Vector3::new(
        -0.5 + (position.y - 0.5).abs().clamp(0., 0.25),
        w0 * f64::sin(4. * std::f64::consts::PI * position.x)
            * (f64::exp(-(position.y - 0.25).powi(2) / (2. * sigma))
                + f64::exp(-(position.y - 0.75).powi(2) / (2. * sigma))),
        0.,
    );
    let p = 2.5;
    let gamma = 1.4;
    Properties {
        rho,
        rho_u: rho * u,
        E: p / (gamma - 1.) + 0.5 * rho * u.dot(&u),
        u,
        e: p / (gamma - 1.) / rho + 0.5 * u.dot(&u),
        p,
        gamma,
    }
}

///////////////////
// Mesh creation //
///////////////////

pub fn create_mesh() -> problem::Problem<Properties> {
    let bc1 = BoundaryCondition("periodic".to_string(), dvector![1., 0., 0.]);
    let bc2 = BoundaryCondition("wall".to_string(), dvector![0.]);
    let bc3 = BoundaryCondition("periodic".to_string(), dvector![-1., 0., 0.]);
    let boundary_conditions: [Option<BoundaryCondition>; 6] = [
        Some(bc1),
        Some(bc2.clone()),
        None,
        Some(bc3),
        Some(bc2.clone()),
        None,
    ];
    let n_grid = (128, 128, 1);
    let mut block: problem::Problem<Properties> = problem::Problem::new();
    let mass_eq = generate_equation_structure!(get_rho, update_rho, euler_mass, 1);
    let impulse_eq = generate_equation_structure!(get_impulse, update_impulse, euler_impulse, 3);
    let energy_eq = generate_equation_structure!(get_energy, update_energy, euler_energy, 1);
    block.add_zone(
        (0., 0., 0.),
        (1., 1., 1.),
        n_grid,
        boundary_conditions,
        vec![mass_eq, impulse_eq, energy_eq],
        properties_generator,
    );
    assert!(block.check_periodic_bc());
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

pub fn get_rho(p: &Properties) -> f64 {
    p.rho
}

pub fn update_rho(p: &mut Properties, rho: f64) {
    p.rho = rho;
}

pub fn euler_mass(problem: &problem::Problem<Properties>, i_p: usize) -> f64 {
    // mass conservation equation
    // d(rho)/dt + div(rho * u) = 0
    let p = &problem.mesh.cells[i_p];
    let mut convective_flux = 0.;
    for neighbor in problem.neighbors(i_p) {
        match neighbor {
            problem::NeighborReturn::Internal(i_n, i_f, _) => {
                let n = &problem.mesh.cells[i_n];
                let f = &problem.mesh.faces[i_f];
                let rho_avg = (p.properties.rho + n.properties.rho) / 2.;
                let u_avg: Vector3<f64> = (p.properties.u + n.properties.u) / 2.;
                let invert_if_not_owner = if f.is_owner(i_p) { 1. } else { -1. };
                convective_flux -= (rho_avg * u_avg).dot(&f.normal) * f.area * invert_if_not_owner;
            }
            problem::NeighborReturn::Boundary(name, i_f, _, values) => {
                let f = &problem.mesh.faces[i_f];
                match name {
                    "inlet" => {
                        // flux with known properties rho (values[0]) and u (values[1..=3])
                        let rho_avg = (p.properties.rho + values[0]) / 2.;
                        let u_avg: Vector3<f64> =
                            (p.properties.u + Vector3::new(values[1], values[2], values[3])) / 2.;
                        convective_flux -= (rho_avg * u_avg).dot(&f.normal) * f.area;
                    }
                    "wall" => {
                        // impermeable wall
                        convective_flux -= 0.;
                    }
                    "fully_developed" => {
                        // fully developed flow which means the flow is assumed to not change in the normal direction
                        convective_flux -=
                            (p.properties.rho * p.properties.u).dot(&f.normal) * f.area;
                    }
                    _ => {
                        panic!("Unknown boundary condition : {}", name);
                    }
                }
            }
            problem::NeighborReturn::Nil => {
                // nothing to do
            }
        }
    }
    convective_flux
}

pub fn get_impulse(p: &Properties) -> Vector3<f64> {
    p.rho_u
}

pub fn update_impulse(p: &mut Properties, rho_u: &Vector3<f64>) {
    p.rho_u[0] = rho_u[0];
    p.rho_u[1] = rho_u[1];
    p.rho_u[2] = rho_u[2];
}

pub fn euler_impulse(problem: &problem::Problem<Properties>, i_p: usize) -> Vector3<f64> {
    // momentum conservation equation
    // d(rho*u)/dt + div(rho * u x u + p * I) = 0
    let p = &problem.mesh.cells[i_p];
    let mut convective_flux = Vector3::from_vec(vec![0., 0., 0.]);

    for neighbor in problem.neighbors(i_p) {
        match neighbor {
            problem::NeighborReturn::Internal(i_n, i_f, _) => {
                let n = &problem.mesh.cells[i_n];
                let f = &problem.mesh.faces[i_f];
                let rho_avg = (p.properties.rho + n.properties.rho) / 2.;
                let u_avg: Vector3<f64> = (p.properties.u + n.properties.u) / 2.;
                let p_avg = (p.properties.p + n.properties.p) / 2.;
                let invert_if_not_owner = if f.is_owner(i_p) { 1. } else { -1. };
                convective_flux -=
                    (rho_avg * u_avg * u_avg.transpose() + p_avg * Matrix3::identity()).transpose()
                        * f.normal
                        * f.area
                        * invert_if_not_owner;
                // TODO : check if this is correct...
            }
            problem::NeighborReturn::Boundary(name, i_f, _, values) => {
                let f = &problem.mesh.faces[i_f];
                match name {
                    "inlet" => {
                        // flux with known properties rho (values[0]) and u (values[1..=3])
                        let rho_avg = (p.properties.rho + values[0]) / 2.;
                        let u_avg: Vector3<f64> =
                            (p.properties.u + Vector3::new(values[1], values[2], values[3])) / 2.;
                        let p_avg = (p.properties.p + values[4]) / 2.;
                        convective_flux -= (rho_avg * u_avg * u_avg.transpose()
                            + p_avg * Matrix3::identity())
                        .transpose()
                            * f.normal
                            * f.area;
                    }
                    "wall" => {
                        // impermeable wall
                        convective_flux -= Vector3::zeros();
                    }
                    "fully_developed" => {
                        // fully developed flow which means the flow is assumed to not change in the normal direction
                        convective_flux -=
                            (p.properties.rho * p.properties.u * p.properties.u.transpose()
                                + p.properties.p * Matrix3::identity())
                            .transpose()
                                * f.normal
                                * f.area;
                    }
                    _ => {
                        panic!("Unknown boundary condition : {}", name);
                    }
                }
            }
            problem::NeighborReturn::Nil => {
                // nothing to do
            }
        }
    }
    convective_flux
}

#[allow(non_snake_case)]
pub fn get_energy(p: &Properties) -> f64 {
    p.E
}

#[allow(non_snake_case)]
pub fn update_energy(p: &mut Properties, E: f64) {
    p.E = E;
}

#[allow(non_snake_case)]
pub fn euler_energy(problem: &problem::Problem<Properties>, i_p: usize) -> f64 {
    // energy conservation equation
    // d(E) + div(u * (E + p)) = 0
    let p = &problem.mesh.cells[i_p];
    let mut convective_flux = 0.;

    for neighbors in problem.neighbors(i_p) {
        match neighbors {
            problem::NeighborReturn::Internal(i_n, i_f, _) => {
                let n = &problem.mesh.cells[i_n];
                let f = &problem.mesh.faces[i_f];
                let E_avg = (p.properties.E + n.properties.E) / 2.;
                let p_avg = (p.properties.p + n.properties.p) / 2.;
                let u_avg: Vector3<f64> = (p.properties.u + n.properties.u) / 2.;
                let invert_if_not_owner = if f.is_owner(i_p) { 1. } else { -1. };
                convective_flux -=
                    (u_avg * (E_avg + p_avg)).dot(&f.normal) * f.area * invert_if_not_owner;
            }
            problem::NeighborReturn::Boundary(name, i_f, _, values) => {
                let f = &problem.mesh.faces[i_f];

                match name {
                    "inlet" => {
                        let boundary_rho = values[0];
                        let boundary_u = Vector3::new(values[1], values[2], values[3]);
                        let boundary_p = values[4];
                        // WARNING : gamma should actually also be passed as part of the boundary condition data
                        let boundary_E = boundary_p / (p.properties.gamma - 1.)
                            + 1. / 2. * boundary_rho * boundary_u.dot(&boundary_u);
                        // flux with known properties rho (values[0]) and u (values[1..=3])
                        let E_avg = (p.properties.E + boundary_E) / 2.;
                        let p_avg = (p.properties.p + boundary_p) / 2.;
                        let u_avg: Vector3<f64> = (p.properties.u + boundary_u) / 2.;
                        convective_flux -= (u_avg * (E_avg + p_avg)).dot(&f.normal) * f.area;
                    }
                    "wall" => {
                        // impermeable wall, we assume it is also adiabatic
                        convective_flux -= 0.;
                    }
                    "fully_developed" => {
                        // fully developed flow which means the flow is assumed to not change in the normal direction
                        convective_flux -= (p.properties.u * (p.properties.E + p.properties.p))
                            .dot(&f.normal)
                            * f.area;
                    }
                    _ => {
                        panic!("Unknown boundary condition : {}", name);
                    }
                }
            }
            problem::NeighborReturn::Nil => {
                // nothing to do
            }
        }
    }
    convective_flux
}
