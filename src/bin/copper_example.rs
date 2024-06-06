use flex_fvm::solver::SolverStrategy;
use flex_fvm::{computational_cell::HasProperties, equation::*, problem, solver};
use nalgebra::{dvector, DVector, Point3, RowDVector};

fn main() {
    let initial_mesh = create_mesh();
    let output_args = solver::SolverOutputArgs {
        frequency: 0.0001,
        filename: "copper_time_series".to_string(),
        print_csv: false,
        print_vtk: true,
    };
    let _solution = solver::solve(
        &initial_mesh,
        0.01,
        0.0000002,
        solver::EulerExplicitSolver::new(&initial_mesh),
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
    pub T: RowDVector<f64>,   // K
    pub rho: RowDVector<f64>, // kg/m3
    pub cp: RowDVector<f64>,  // J/k/kg
    // solid
    pub lambda: RowDVector<f64>, // W/m/k
}

impl HasProperties for Properties {
    fn new(size: usize) -> Self {
        if size == 0 {
            return Properties {
                T: RowDVector::from(Vec::new()),
                rho: RowDVector::from(Vec::new()),
                cp: RowDVector::from(Vec::new()),
                lambda: RowDVector::from(Vec::new()),
            };
        }
        Properties {
            T: RowDVector::from_element(size, 0.),
            rho: RowDVector::from_element(size, 0.),
            cp: RowDVector::from_element(size, 0.),
            lambda: RowDVector::from_element(size, 0.),
        }
    }

    fn export(&self) -> Vec<(String, Vec<f64>)> {
        vec![("T".to_string(), self.T.data.as_slice().to_vec())]
    }
}

pub fn copper(properties: &mut Properties, _position: &Point3<f64>) {
    properties.T = properties.T.transpose().push(0.).transpose();
    properties.rho = properties.rho.transpose().push(9000.).transpose();
    properties.cp = properties.cp.transpose().push(0.001).transpose();
    properties.lambda = properties.lambda.transpose().push(300.).transpose();
}

///////////////////////
/// Mesh Generation ///
///////////////////////

#[allow(non_snake_case)]
pub fn create_mesh() -> problem::Problem<Properties> {
    // pack the equations
    let heat_diffusion_equations = Equation {
        get: Properties::get_T,
        lhs: Properties::update_T,
        source: Some(heat_diffusion),
        flux: None,
        dimension: 1,
    };
    // define the boundary conditions
    // left
    let bc1 = BoundaryCondition("dirichlet".to_string(), dvector![300.]);
    // front, back, top, bottom
    let bc2 = BoundaryCondition("neumann".to_string(), dvector![0.]);
    // right
    let bc3 = BoundaryCondition("dirichlet".to_string(), dvector![0.]);
    // top
    let bc4 = BoundaryCondition("neumann".to_string(), dvector![100_000.]);

    let n_grid = (100, 1, 40);

    let mut block: problem::Problem<Properties> = problem::Problem::new();
    // bottom left
    let boundary_conditions = [
        Some(bc1.clone()),
        Some(bc2.clone()),
        Some(bc2.clone()),
        Some(bc3.clone()),
        Some(bc2.clone()),
        Some(bc4.clone()),
    ];
    block.add_zone(
        (0., 0., 0.),
        (1., 1., 1.),
        n_grid,
        boundary_conditions,
        vec![heat_diffusion_equations],
        copper,
    );
    block
}

///////////////////////////////
/// Heat diffusion equation ///
///////////////////////////////
#[allow(non_snake_case)]
impl Properties {
    pub fn get_T(&self, i: usize) -> DVector<f64> {
        dvector![self.T[i]]
    }

    pub fn update_T(&mut self, i: usize, T: &[f64]) {
        self.T[i] = T[0];
    }
}

pub fn heat_diffusion(
    problem: &problem::Problem<Properties>,
    i_p: usize,
    output: &mut DVector<f64>,
) {
    let flux = &mut output[0];
    for neighbor in problem.neighbors(i_p) {
        match neighbor {
            problem::NeighborReturn::Internal(i_n, i_f, distance) => {
                let f = &problem.mesh.faces.get_face(i_f);
                *flux -= problem.cell_properties.lambda[i_p]
                    * (problem.cell_properties.T[i_p] - problem.cell_properties.T[i_n])
                    / distance
                    * f.area;
            }
            problem::NeighborReturn::Boundary(name, i_f, distance, values) => {
                let f = &problem.mesh.faces.get_face(i_f);
                match name {
                    "dirichlet" => {
                        *flux -= problem.cell_properties.lambda[i_p]
                            * (problem.cell_properties.T[i_p] - values[0])
                            / distance
                            * f.area;
                    }
                    "neumann" => {
                        *flux -= values[0] * f.area;
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
    *flux /= problem.mesh.cells.volume[i_p]
        * problem.cell_properties.rho[i_p]
        * problem.cell_properties.cp[i_p]
}
