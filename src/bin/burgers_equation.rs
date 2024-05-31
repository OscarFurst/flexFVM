use flex_fvm::{computational_cell::HasProperties, equation::*, problem, solver};
use nalgebra::{dvector, DVector, Matrix3xX, Point3, RowDVector, Vector3};

fn main() {
    let initial_mesh = create_mesh();
    let output_args = solver::SolverOutputArgs {
        frequency: 0.1,
        filename: "burger_time_series".to_string(),
        print_csv: true,
        print_vtk: false,
    };
    let _solution = solver::solve(
        &initial_mesh,
        50.,
        0.0001,
        solver::euler_explicit::<Properties>,
        output_args,
    );
}

////////////////////////
/// Cell Propertires ///
////////////////////////

#[allow(non_snake_case)]
#[derive(Clone)]
pub struct Properties {
    // physical properties
    pub u: Matrix3xX<f64>,   // K
    pub nu: RowDVector<f64>, // kg/m3
}

impl HasProperties for Properties {
    fn new(size: usize) -> Self {
        if size == 0 {
            return Properties {
                u: Matrix3xX::from_vec(vec![]),
                nu: RowDVector::from_vec(vec![]),
            };
        }
        Properties {
            u: Matrix3xX::from_element(size, 0.),
            nu: RowDVector::from_element(size, 0.),
        }
    }
    fn export(&self) -> Vec<(String, Vec<f64>)> {
        vec![
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
        ]
    }
}

pub fn push_to_matrix(matrix: &Matrix3xX<f64>, value: Vector3<f64>) -> Matrix3xX<f64> {
    Matrix3xX::from_vec([matrix.data.as_slice(), value.data.as_slice()].concat())
}

pub fn material(properties: &mut Properties, position: &Point3<f64>) {
    properties.u = push_to_matrix(
        &properties.u,
        Vector3::new(f64::exp(-position.x.powi(2) / 2.), 0., 0.),
    );
    properties.nu = properties.nu.transpose().push(1.0).transpose();
}

///////////////////////
/// Mesh Generation ///
///////////////////////

#[allow(non_snake_case)]
pub fn create_mesh() -> problem::Problem<Properties> {
    // pack the equations
    let heat_diffusion_equations = Equation {
        get: Properties::get_u,
        lhs: Properties::update_u,
        source: burgers_equation,
        dimension: 3,
    };
    // define the boundary conditions
    // left
    let bc1 = BoundaryCondition("periodic".to_string(), dvector![1., 0., 0.]);
    // right
    let bc2 = BoundaryCondition("periodic".to_string(), dvector![-1., 0., 0.]);

    let n_grid = (500, 1, 1);

    let mut block: problem::Problem<Properties> = problem::Problem::new();
    // bottom left
    let boundary_conditions = [Some(bc1), None, None, Some(bc2), None, None];
    block.add_zone(
        (-7., 0., 0.),
        (7., 1., 1.),
        n_grid,
        boundary_conditions,
        vec![heat_diffusion_equations],
        material,
    );
    block
}

/////////////////////////
/// Burgers' equation ///
/////////////////////////
impl Properties {
    pub fn get_u(&self, i: usize) -> DVector<f64> {
        dvector![self.u[(0, i)], self.u[(1, i)], self.u[(2, i)]]
    }

    pub fn update_u(&mut self, i: usize, u: &[f64]) {
        self.u.column_mut(i).copy_from_slice(u);
    }

    pub fn upwind(&self, left: usize, right: usize, normal: Vector3<f64>) -> usize {
        if self.u.column(left).dot(&normal) >= 0. {
            left
        } else {
            right
        }
    }
}

pub fn burgers_equation(problem: &problem::Problem<Properties>, i_p: usize) -> DVector<f64> {
    // du/dt + div(u^2)/2 = 0
    let mut advective_flux = dvector![0., 0., 0.];
    for neighbor in problem.neighbors(i_p) {
        match neighbor {
            // flux term for a normal neighbor
            problem::NeighborReturn::Internal(i_n, i_f, _) => {
                let f = &problem.mesh.faces[i_f];
                let i_up = problem.cell_properties.upwind(i_p, i_n, f.out_normal(i_p));
                let u_upwind = problem.cell_properties.u.column(i_up);
                advective_flux -= 0.5 * u_upwind.dot(&u_upwind) * f.area * f.out_normal(i_p);
            }
            // flux term for a boundary condition
            problem::NeighborReturn::Boundary(name, _, _, _) => match name {
                // currently only periodic boundary condition is used in this case,
                // which does not require any action here
                _ => {
                    panic!("Unknown boundary condition : {}", name);
                }
            },
            problem::NeighborReturn::Nil => {
                // nothing to do
            }
        }
    }
    advective_flux / problem.mesh.cells.volume[i_p]
}
