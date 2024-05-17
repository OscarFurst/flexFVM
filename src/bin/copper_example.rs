use flex_fvm::{
    computational_cell::HasProperties, equation::*, generate_equation_structure, problem, solver,
};
use nalgebra::{dvector, DVector, Point3};
use std::fmt;

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
        solver::euler_explicit::<Properties>,
        output_args,
    );
}

////////////////////////
/// Cell Propertires ///
////////////////////////

#[allow(non_snake_case)]
#[derive(Default, Clone, Copy)]
pub struct Properties {
    // physical properties
    pub T: f64,   // K
    pub rho: f64, // kg/m3
    pub cp: f64,  // J/k/kg
    // solid
    pub lambda: f64, // W/m/k
}

impl fmt::Debug for Properties {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.T)
    }
}

impl HasProperties for Properties {
    fn update_primitive(&mut self, volume: f64) {
        // here we can update derived variables if needed
    }

    fn export(&self) -> Vec<(String, f64)> {
        vec![("T".to_string(), self.T)]
    }
}

pub fn copper(_position: &Point3<f64>) -> Properties {
    Properties {
        T: 0.,
        rho: 9000.,
        cp: 0.001,
        lambda: 300.,
        //..Default::default()
    }
}

///////////////////////
/// Mesh Generation ///
///////////////////////

#[allow(non_snake_case)]
pub fn create_mesh() -> problem::Problem<Properties> {
    // pack the equations
    let heat_diffusion_equations =
        generate_equation_structure!(get_T, update_T_user, heat_diffusion_user, 1);
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
pub fn get_T(p: &Properties) -> f64 {
    p.T
}

#[allow(non_snake_case)]
pub fn update_T_user(p: &mut Properties, T: f64) {
    p.T = T;
}

pub fn heat_diffusion_user(problem: &problem::Problem<Properties>, i_p: usize) -> f64 {
    // TODO : deal with boundary conditions
    let p = &problem.mesh.cells[i_p];
    let mut flux = 0.;
    for neighbor in problem.neighbors(i_p) {
        match neighbor {
            // flux term for a normal neighbor
            problem::NeighborReturn::Internal(i_n, i_f, distance) => {
                let n = &problem.mesh.cells[i_n];
                let f = &problem.mesh.faces[i_f];
                flux -= p.properties.lambda * (p.properties.T - n.properties.T) / distance * f.area;
            }
            // flux term for a boundary condition
            problem::NeighborReturn::Boundary(name, i_f, distance, values) => {
                let f = &problem.mesh.faces[i_f];
                match name {
                    "dirichlet" => {
                        flux -=
                            p.properties.lambda * (p.properties.T - values[0]) / distance * f.area;
                    }
                    "neumann" => {
                        flux -= values[0] * f.area;
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
    flux / p.volume / p.properties.rho / p.properties.cp
}
