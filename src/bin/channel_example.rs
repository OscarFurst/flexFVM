use flex_fvm::computational_cell::{HasProperties, NeighborType};
use flex_fvm::problem;
use flex_fvm::{
    equation::{BoundaryCondition, Equation},
    solver,
};
use nalgebra::{dvector, DVector, Matrix3, Matrix3xX, Point3, RowDVector, Vector3};
use nalgebra_sparse::csr::CsrMatrix;

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
        0.001,
        solver::runge_kutta_4::<Properties>,
        output_args,
    );
}

///////////////////////
// Cell Propertires ///
///////////////////////

#[allow(non_snake_case)]
#[derive(Clone)]
pub struct Properties {
    // conserved properties
    pub mass: RowDVector<f64>,    // kg
    pub momentum: Matrix3xX<f64>, // kg/m/s
    pub energy: RowDVector<f64>,  // J
    // primitive properties
    pub rho: RowDVector<f64>,             // kg/m3
    pub u: Matrix3xX<f64>,                // m/s
    pub p: RowDVector<f64>,               // Pa
    pub specific_energy: RowDVector<f64>, // J/kg
    // material properties
    pub gamma: RowDVector<f64>, // adiabatic index
}

impl HasProperties for Properties {
    fn new(size: usize) -> Self {
        if size == 0 {
            return Properties {
                mass: RowDVector::from(Vec::new()),
                momentum: Matrix3xX::from_vec(Vec::new()),
                energy: RowDVector::from(Vec::new()),
                rho: RowDVector::from(Vec::new()),
                u: Matrix3xX::from_vec(Vec::new()),
                p: RowDVector::from(Vec::new()),
                specific_energy: RowDVector::from(Vec::new()),
                gamma: RowDVector::from(Vec::new()),
            };
        }
        Properties {
            mass: RowDVector::from_element(size, 0.),
            momentum: Matrix3xX::from_element(size, 0.),
            energy: RowDVector::from_element(size, 0.),
            rho: RowDVector::from_element(size, 0.),
            u: Matrix3xX::from_element(size, 0.),
            p: RowDVector::from_element(size, 0.),
            specific_energy: RowDVector::from_element(size, 0.),
            gamma: RowDVector::from_element(size, 0.),
        }
    }

    fn update_primitive(&mut self, volume: &RowDVector<f64>) {
        // conserved properties : mass, momentum, energy
        self.rho = self.mass.component_div(volume);
        for i in 0..3 {
            self.u.set_row(
                i,
                &self
                    .momentum
                    .row(i)
                    .component_div(&self.rho)
                    .component_div(volume),
            );
        }
        let list_of_velocity_magnitude: RowDVector<f64> = self.u.component_mul(&self.u).row_sum();
        self.p = (self.energy.component_div(volume)
            - 0.5 * self.rho.component_mul(&list_of_velocity_magnitude))
        .component_mul(&self.gamma.add_scalar(-1.));
        self.specific_energy = self.p.component_div(&self.gamma.add_scalar(-1.))
            + 0.5 * self.rho.component_mul(&list_of_velocity_magnitude);
    }

    fn update_integral(&mut self, volume: &RowDVector<f64>) {
        // conserved properties : mass, momentum, energy
        self.mass = self.rho.component_mul(volume);
        for i in 0..3 {
            self.momentum
                .set_row(i, &self.mass.component_mul(&self.u.row(i)));
        }
        let list_of_velocity_magnitude: RowDVector<f64> = self.u.component_mul(&self.u).row_sum();
        self.energy = (self.p.component_div(&self.gamma.add_scalar(-1.))
            + 0.5 * self.rho.component_mul(&list_of_velocity_magnitude))
        .component_mul(volume);
    }

    fn interpolate_face_from_cells(
        &mut self,
        cell_properties: &Properties,
        interpolation_matrix: &CsrMatrix<f64>,
    ) {
        // interpolate only the primitive properties, as the integral properties depend on the cell volume!
        // TODO : this is probably terribly inefficient because of the transpose, which are only needed
        // because left multiplication is not implemented in nalgebra-sparse!
        self.rho = (interpolation_matrix * &cell_properties.rho.transpose()).transpose();
        self.u = (interpolation_matrix * &cell_properties.u.transpose()).transpose();
        self.p = (interpolation_matrix * &cell_properties.p.transpose()).transpose();
        self.gamma = (interpolation_matrix * &cell_properties.gamma.transpose()).transpose();
        // the following is a "derived property" and is therefore recomputed to avoid incosisntency
        let list_of_velocity_magnitude: RowDVector<f64> = self.u.component_mul(&self.u).row_sum();
        self.specific_energy = self.p.component_div(&self.gamma.add_scalar(-1.))
            + 0.5 * self.rho.component_mul(&list_of_velocity_magnitude);
    }

    fn export(&self) -> Vec<(String, Vec<f64>)> {
        vec![
            ("rho".to_string(), self.rho.data.as_slice().to_vec()),
            (
                "ux".to_string(),
                self.u.row(0).clone_owned().data.as_slice().to_vec(),
            ),
            (
                "uy".to_string(),
                self.u.row(1).clone_owned().data.as_slice().to_vec(),
            ),
            (
                "uz".to_string(),
                self.u.row(2).clone_owned().data.as_slice().to_vec(),
            ),
            ("p".to_string(), self.p.data.as_slice().to_vec()),
        ]
    }
}

pub fn push_to_matrix(matrix: &Matrix3xX<f64>, value: Vector3<f64>) -> Matrix3xX<f64> {
    Matrix3xX::from_vec([matrix.data.as_slice(), value.data.as_slice()].concat())
}

fn properties_generator(properties: &mut Properties, position: &Point3<f64>) {
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
    properties.rho = properties.rho.transpose().push(rho).transpose();
    properties.u = push_to_matrix(&properties.u, u);
    properties.p = properties.p.transpose().push(p).transpose();
    properties.gamma = properties.gamma.transpose().push(gamma).transpose();
    // will be updated by update_primitive:
    properties.mass = properties.mass.transpose().push(0.).transpose();
    properties.momentum = push_to_matrix(&properties.momentum, Vector3::new(0., 0., 0.));
    properties.energy = properties.energy.transpose().push(0.).transpose();
}

///////////////////
// Mesh creation //
///////////////////

pub fn create_mesh() -> problem::Problem<Properties> {
    let bc1 = BoundaryCondition("periodic".to_string(), dvector![1., 0., 0.]);
    let bc3 = BoundaryCondition("periodic".to_string(), dvector![-1., 0., 0.]);
    let bc4 = BoundaryCondition("periodic".to_string(), dvector![0., 1., 0.]);
    let bc5 = BoundaryCondition("periodic".to_string(), dvector![0., -1., 0.]);
    let boundary_conditions: [Option<BoundaryCondition>; 6] =
        [Some(bc1), Some(bc4), None, Some(bc3), Some(bc5), None];
    let n_grid = (128, 128, 1);
    let mut block: problem::Problem<Properties> = problem::Problem::new();
    let mass_eq = Equation {
        get: Properties::get_mass,
        lhs: Properties::update_mass,
        source: None,
        flux: Some(euler_mass),
        dimension: 1,
    };
    let impulse_eq = Equation {
        get: Properties::get_momentum,
        lhs: Properties::update_momentum,
        source: None,
        flux: Some(euler_impulse),
        dimension: 3,
    };
    let energy_eq = Equation {
        get: Properties::get_energy,
        lhs: Properties::update_energy,
        source: None,
        flux: Some(euler_energy),
        dimension: 1,
    };
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
        .cell_properties
        .update_integral(&block.mesh.cells.volume);
    block.prepare_MUSCL_interpolation();
    block
}

/////////////////////////
/// Numerical methods ///
/////////////////////////

fn rusanov_scalar(
    owner_properties: &Properties,
    neighbor_properties: &Properties,
    flux: fn(&Properties, &Properties, usize) -> Vector3<f64>,
    var: fn(&Properties, usize) -> f64,
    i_f: usize,
    scaled_normal: &Vector3<f64>,
) -> f64 {
    let f = flux(&owner_properties, &neighbor_properties, i_f);
    let area = scaled_normal.norm();
    let lambda = f64::max(
        // we take u_x to match the source but eventually we should take the normal velocity
        owner_properties.u[(0, i_f)].abs() * area
            + f64::sqrt(
                owner_properties.gamma[i_f] * owner_properties.p[i_f] / owner_properties.rho[i_f],
            ) * area,
        neighbor_properties.u[(0, i_f)].abs() * area
            + f64::sqrt(
                neighbor_properties.gamma[i_f] * neighbor_properties.p[i_f]
                    / neighbor_properties.rho[i_f],
            ) * area,
    );
    let u1 = var(owner_properties, i_f);
    let u2 = var(neighbor_properties, i_f);
    f.dot(scaled_normal) - 0.5 * lambda * (u2 - u1)
}

fn rusanov_vector(
    owner_properties: &Properties,
    neighbor_properties: &Properties,
    flux: fn(&Properties, &Properties, usize) -> Matrix3<f64>,
    var: fn(&Properties, usize) -> Vector3<f64>,
    i_f: usize,
    scaled_normal: &Vector3<f64>,
) -> Vector3<f64> {
    let f = flux(&owner_properties, &neighbor_properties, i_f);
    let area = scaled_normal.norm();
    let lambda = f64::max(
        owner_properties.u[(0, i_f)].abs() * area
            + f64::sqrt(
                owner_properties.gamma[i_f] * owner_properties.p[i_f] / owner_properties.rho[i_f],
            ) * area,
        neighbor_properties.u[(0, i_f)].abs() * area
            + f64::sqrt(
                neighbor_properties.gamma[i_f] * neighbor_properties.p[i_f]
                    / neighbor_properties.rho[i_f],
            ) * area,
    );
    let u1 = var(owner_properties, i_f);
    let u2 = var(neighbor_properties, i_f);
    f * scaled_normal - 0.5 * lambda * (u2 - u1)
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

impl Properties {
    pub fn get_mass(&self, i: usize) -> DVector<f64> {
        dvector![self.mass[i]]
    }

    pub fn update_mass(&mut self, i: usize, mass: &[f64]) {
        self.mass[i] = mass[0];
    }

    pub fn get_rho(&self, i: usize) -> f64 {
        self.rho[i]
    }

    fn euler_mass_flux(&self, other: &Properties, i: usize) -> Vector3<f64> {
        0.5 * (self.rho[i] * self.u.column(i) + other.rho[i] * other.u.column(i))
    }
}

pub fn euler_mass(problem: &problem::Problem<Properties>, i_f: usize) -> DVector<f64> {
    // mass conservation equation
    // d(rho)/dt + div(rho * u) = 0
    let mut convective_flux = 0.;
    let f = &problem.mesh.faces[i_f];
    match problem.neighbor(i_f) {
        problem::NeighborReturn::Internal(_i_n, _i_f, _d) => {
            convective_flux -= rusanov_scalar(
                &problem.face_owner_properties,
                &problem.face_neighbor_properties,
                Properties::euler_mass_flux,
                Properties::get_rho,
                i_f,
                &(f.normal * f.area),
            );
        }
        problem::NeighborReturn::Boundary(name, _i_f, _d, _val) => {
            match name {
                "wall" => {
                    // impermeable wall, nothing to do
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
    dvector![convective_flux]
}

impl Properties {
    pub fn get_momentum(&self, i: usize) -> DVector<f64> {
        // there has to be a better way to do this...
        dvector![
            self.momentum[(0, i)],
            self.momentum[(1, i)],
            self.momentum[(2, i)],
        ]
    }

    pub fn update_momentum(&mut self, i: usize, momentum: &[f64]) {
        self.momentum.column_mut(i).copy_from_slice(momentum);
    }

    pub fn get_impulse(&self, i: usize) -> Vector3<f64> {
        self.rho[i] * self.u.column(i)
    }

    fn euler_impulse_flux(&self, other: &Properties, i: usize) -> Matrix3<f64> {
        let mom_star = 0.5 * (self.rho[i] * self.u.column(i) + other.rho[i] * other.u.column(i));
        let rho_star = 0.5 * (self.rho[i] + other.rho[i]);
        let en_star = 0.5 * (self.specific_energy[i] + other.specific_energy[i]);
        let p_star = (self.gamma[i] - 1.0) * (en_star - 0.5 * mom_star.dot(&mom_star) / rho_star);
        mom_star * mom_star.transpose() / rho_star + p_star * Matrix3::identity()
    }
}

pub fn euler_impulse(problem: &problem::Problem<Properties>, i_f: usize) -> DVector<f64> {
    // momentum conservation equation
    // d(rho*u)/dt + div(1/rho * rho*u x rho*u + p * I) = 0
    let mut convective_flux = dvector![0., 0., 0.];
    let f = &problem.mesh.faces[i_f];
    match problem.neighbor(i_f) {
        problem::NeighborReturn::Internal(_i_n, i_f, _) => {
            convective_flux -= rusanov_vector(
                &problem.face_owner_properties,
                &problem.face_neighbor_properties,
                Properties::euler_impulse_flux,
                Properties::get_impulse,
                i_f,
                &(f.normal * f.area),
            );
        }
        problem::NeighborReturn::Boundary(name, i_f, _, _) => {
            match name {
                "wall" => {
                    convective_flux -= (problem.face_owner_properties.p[i_f] * Matrix3::identity())
                        * f.normal
                        * f.area;
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
    convective_flux
}

impl Properties {
    pub fn get_energy(&self, i: usize) -> DVector<f64> {
        dvector![self.energy[i]]
    }

    pub fn update_energy(&mut self, i: usize, energy: &[f64]) {
        self.energy[i] = energy[0];
    }

    pub fn get_specific_energy(&self, i: usize) -> f64 {
        self.specific_energy[i]
    }

    fn euler_energy_flux(&self, other: &Properties, i: usize) -> Vector3<f64> {
        let mom_star = 0.5 * (self.rho[i] * self.u.column(i) + other.rho[i] * other.u.column(i));
        let rho_star = 0.5 * (self.rho[i] + other.rho[i]);
        let en_star = 0.5 * (self.specific_energy[i] + other.specific_energy[i]);
        let p_star = (self.gamma[i] - 1.0) * (en_star - 0.5 * mom_star.dot(&mom_star) / rho_star);
        return mom_star / rho_star * (en_star + p_star);
    }
}

#[allow(non_snake_case)]
pub fn euler_energy(problem: &problem::Problem<Properties>, i_f: usize) -> DVector<f64> {
    // energy conservation equation
    // d(E) + div(u * (E + p)) = 0
    let mut convective_flux = 0.;

    let f = &problem.mesh.faces[i_f];
    match problem.neighbor(i_f) {
        problem::NeighborReturn::Internal(_i_n, i_f, _) => {
            convective_flux -= rusanov_scalar(
                &problem.face_owner_properties,
                &problem.face_neighbor_properties,
                Properties::euler_energy_flux,
                Properties::get_specific_energy,
                i_f,
                &(f.normal * f.area),
            );
        }
        problem::NeighborReturn::Boundary(name, _i_f, _, _) => {
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
    dvector![convective_flux]
}
