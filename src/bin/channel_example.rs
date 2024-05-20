use flex_fvm::computational_cell::HasProperties;
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
// Good explanation of numerical scheme can be found here : https://www2.math.ethz.ch/education/bachelor/lectures/fs2014/other/n_dgl/chap8.pdf
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
        print_csv: false,
        print_vtk: true,
    };
    let _solution = solver::solve(
        &initial_mesh,
        2.,
        0.0001,
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
    // conserved properties
    pub mass: f64,              // kg
    pub momentum: Vector3<f64>, // kg/m/s
    pub energy: f64,            // J
    // primitive properties
    pub rho: f64,        // kg/m3
    pub u: Vector3<f64>, // m/s
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
    fn update_primitive(&mut self, volume: f64) {
        // conserved properties : mass, momentum, energy
        self.rho = self.mass / volume;
        self.u = self.momentum / self.rho / volume;
        self.p = (self.energy / volume - 0.5 * self.rho * self.u.dot(&self.u)) * (self.gamma - 1.);
    }

    fn update_integral(&mut self, volume: f64) {
        // conserved properties : mass, momentum, energy
        self.mass = self.rho * volume;
        self.momentum = self.rho * self.u * volume;
        self.energy = (self.p / (self.gamma - 1.) + 0.5 * self.rho * self.u.dot(&self.u)) * volume;
    }

    fn export(&self) -> Vec<(String, f64)> {
        vec![
            ("rho".to_string(), self.rho),
            ("ux".to_string(), self.u[0]),
            ("uy".to_string(), self.u[1]),
            ("uz".to_string(), self.u[2]),
            ("p".to_string(), self.p),
        ]
    }
}

fn properties_generator(position: &Point3<f64>) -> Properties {
    let w0 = 0.1;
    let sigma = 0.05 / f64::sqrt(2.);
    let midsection = |x: f64| if (x - 0.5).abs() < 0.25 { 1. } else { 0. };
    let rho = 1. + midsection(position.y);
    let u = Vector3::new(
        -0.5 + midsection(position.y),
        // w0*np.sin(4*np.pi*X) * ( np.exp(-(Y-0.25)**2/(2 * sigma**2)) + np.exp(-(Y-0.75)**2/(2*sigma**2)) )
        w0 * f64::sin(4. * std::f64::consts::PI * position.x)
            * (f64::exp(-(position.y - 0.25).powf(2.) / (2. * sigma.powf(2.)))
                + f64::exp(-(position.y - 0.75).powf(2.) / (2. * sigma.powf(2.)))),
        0.,
    );
    let p = 2.5;
    let gamma = 5. / 3.;
    Properties {
        rho,
        u,
        p,
        gamma,
        ..Default::default()
    }
}

///////////////////
// Mesh creation //
///////////////////

pub fn create_mesh() -> problem::Problem<Properties> {
    let bc1 = BoundaryCondition("periodic".to_string(), dvector![1., 0., 0.]);
    let bc2 = BoundaryCondition("wall".to_string(), dvector![0.]);
    let bc3 = BoundaryCondition("periodic".to_string(), dvector![-1., 0., 0.]);
    let bc4 = BoundaryCondition("periodic".to_string(), dvector![0., 1., 0.]);
    let bc5 = BoundaryCondition("periodic".to_string(), dvector![0., -1., 0.]);
    let boundary_conditions: [Option<BoundaryCondition>; 6] =
        [Some(bc1), Some(bc4), None, Some(bc3), Some(bc5), None];
    // let boundary_conditions: [Option<BoundaryCondition>; 6] = [
    //     Some(bc1),
    //     Some(bc2.clone()),
    //     None,
    //     Some(bc3),
    //     Some(bc2.clone()),
    //     None,
    // ];
    let n_grid = (128, 128, 1);
    let mut block: problem::Problem<Properties> = problem::Problem::new();
    let mass_eq = generate_equation_structure!(get_mass, update_mass, euler_mass, 1);
    let impulse_eq = generate_equation_structure!(get_momentum, update_momentum, euler_impulse, 3);
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
    block.mesh.update_integral();
    block
}

/////////////////////////
/// Numerical methods ///
/////////////////////////

fn rusanov_scalar(
    flux: fn(&Properties) -> Vector3<f64>,
    var: fn(&Properties) -> f64,
    p1: &Properties,
    p2: &Properties,
    scaled_normal: &Vector3<f64>,
) -> f64 {
    let f1 = flux(p1);
    let f2 = flux(p2);
    let area = scaled_normal.norm();
    let lambda = f64::max(
        p1.u.dot(scaled_normal).abs() + f64::sqrt(p1.gamma * p1.p / p1.rho) * area,
        p2.u.dot(scaled_normal).abs() + f64::sqrt(p2.gamma * p2.p / p2.rho) * area,
    );
    let u1 = var(p1);
    let u2 = var(p2);
    return 0.5 * (f1 + f2).dot(&scaled_normal) - 0.5 * lambda * (u2 - u1);
}

fn rusanov_vector(
    flux: fn(&Properties) -> Matrix3<f64>,
    var: fn(&Properties) -> Vector3<f64>,
    p1: &Properties,
    p2: &Properties,
    scaled_normal: &Vector3<f64>,
) -> Vector3<f64> {
    let f1 = flux(p1);
    let f2 = flux(p2);
    let area = scaled_normal.norm();
    let lambda = f64::max(
        p1.u.dot(scaled_normal).abs() + f64::sqrt(p1.gamma * p1.p / p1.rho) * area,
        p2.u.dot(scaled_normal).abs() + f64::sqrt(p2.gamma * p2.p / p2.rho) * area,
    );
    let u1 = var(p1);
    let u2 = var(p2);
    return 0.5 * (f1 + f2) * scaled_normal - 0.5 * lambda * (u2 - u1);
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

pub fn get_mass(p: &Properties) -> f64 {
    p.mass
}

pub fn update_mass(p: &mut Properties, mass: f64) {
    p.mass = mass;
}

pub fn get_rho(p: &Properties) -> f64 {
    p.rho
}

fn euler_mass_flux(p: &Properties) -> Vector3<f64> {
    return p.rho * p.u;
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
                convective_flux -= rusanov_scalar(
                    euler_mass_flux,
                    get_rho,
                    &p.properties,
                    &n.properties,
                    &(f.out_normal(i_p) * f.area),
                );
            }
            problem::NeighborReturn::Boundary(name, i_f, _, _) => {
                let f = &problem.mesh.faces[i_f];
                match name {
                    "wall" => {
                        // impermeable wall
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

pub fn get_momentum(p: &Properties) -> Vector3<f64> {
    p.momentum
}

pub fn update_momentum(p: &mut Properties, momentum: &Vector3<f64>) {
    p.momentum[0] = momentum[0];
    p.momentum[1] = momentum[1];
    p.momentum[2] = momentum[2];
}

pub fn get_impulse(p: &Properties) -> Vector3<f64> {
    p.rho * p.u
}

fn euler_impulse_flux(p: &Properties) -> Matrix3<f64> {
    return p.rho * p.u * p.u.transpose() + p.p * Matrix3::identity();
}

pub fn euler_impulse(problem: &problem::Problem<Properties>, i_p: usize) -> Vector3<f64> {
    // momentum conservation equation
    // d(rho*u)/dt + div(1/rho * rho*u x rho*u + p * I) = 0
    let p = &problem.mesh.cells[i_p];
    let mut convective_flux = Vector3::new(0., 0., 0.);

    for neighbor in problem.neighbors(i_p) {
        match neighbor {
            problem::NeighborReturn::Internal(i_n, i_f, _) => {
                let n = &problem.mesh.cells[i_n];
                let f = &problem.mesh.faces[i_f];
                convective_flux -= rusanov_vector(
                    euler_impulse_flux,
                    get_impulse,
                    &p.properties,
                    &n.properties,
                    &(f.out_normal(i_p) * f.area),
                );
            }
            problem::NeighborReturn::Boundary(name, i_f, _, _) => {
                let f = &problem.mesh.faces[i_f];
                match name {
                    "wall" => {
                        convective_flux -=
                            (p.properties.p * Matrix3::identity()) * f.out_normal(i_p) * f.area;
                        // impermeable wall
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

pub fn get_energy(p: &Properties) -> f64 {
    p.energy
}

#[allow(non_snake_case)]
pub fn update_energy(p: &mut Properties, energy: f64) {
    p.energy = energy;
}

pub fn get_specific_energy(p: &Properties) -> f64 {
    p.p / (p.gamma - 1.) + 0.5 * p.rho * p.u.dot(&p.u)
}

fn euler_energy_flux(p: &Properties) -> Vector3<f64> {
    let en = get_specific_energy(p);
    return p.u * (en + p.p);
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
                convective_flux -= rusanov_scalar(
                    euler_energy_flux,
                    get_specific_energy,
                    &p.properties,
                    &n.properties,
                    &(f.out_normal(i_p) * f.area),
                );
            }
            problem::NeighborReturn::Boundary(name, i_f, _, _) => {
                let f = &problem.mesh.faces[i_f];

                match name {
                    "wall" => {
                        // impermeable wall, we assume it is also adiabatic
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
