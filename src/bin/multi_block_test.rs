use flex_fvm::computational_cell::{PropertySlice, PropertySliceMut, PropertySliceOwned};
use flex_fvm::{
    computational_cell::HasProperties, equation::*, generate_equation_structure, problem, solver,
};
use nalgebra::{dvector, DVector, Point3};
use std::fmt;
use std::time::Instant;

fn main() {
    let start = Instant::now();
    let initial_mesh = create_mesh();
    let duration = start.elapsed();
    println!(
        "Meshing time : {} ms. Number of cells : {}",
        duration.as_millis(),
        initial_mesh.n_cells()
    );
    let output_args = solver::SolverOutputArgs {
        frequency: 0.,
        filename: "copper_solution".to_string(),
        print_csv: false,
        print_vtk: true,
    };
    let _solution = solver::solve(
        &initial_mesh,
        0.001,
        0.000001,
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
    pub T: DVector<f64>,   // K
    pub rho: DVector<f64>, // kg/m3
    pub cp: DVector<f64>,  // J/k/kg
    // solid
    pub lambda: DVector<f64>, // W/m/k
}

impl fmt::Debug for Properties {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.T)
    }
}

impl HasProperties for Properties {
    fn export(&self) -> Vec<(String, &DVector<f64>)> {
        vec![("T".to_string(), &self.T)]
    }

    fn slice(&self, index: usize) -> PropertySlice {
        PropertySlice {
            T: &self.T[index],
            rho: &self.rho[index],
            cp: &self.cp[index],
            lambda: &self.lambda[index],
        }
    }

    fn slice_mut(&mut self, index: usize) -> PropertySliceMut {
        PropertySliceMut {
            T: &mut self.T[index],
            rho: &mut self.rho[index],
            cp: &mut self.cp[index],
            lambda: &mut self.lambda[index],
        }
    }

    fn push(&mut self, slice: PropertySliceOwned) {
        self.T = self.T.push(slice.T);
        self.rho = self.rho.push(slice.rho);
        self.cp = self.cp.push(slice.cp);
        self.lambda = self.lambda.push(slice.lambda);
    }
}

impl Default for Properties {
    fn default() -> Self {
        Properties {
            T: DVector::from(Vec::new()),
            rho: DVector::from(Vec::new()),
            cp: DVector::from(Vec::new()),
            lambda: DVector::from(Vec::new()),
        }
    }
}
///////////////////
/// Copper data ///
///////////////////

pub fn copper(_position: &Point3<f64>) -> PropertySliceOwned {
    PropertySliceOwned {
        T: 0.,
        rho: 9000.,
        cp: 0.001,
        lambda: 300.,
    }
}

#[allow(non_snake_case)]
pub fn create_mesh() -> problem::Problem<Properties> {
    // pack the equations
    let heat_diffusion_equations =
        generate_equation_structure!(get_T, update_T_user, heat_diffusion_user, 1);
    // define the boundary conditions
    // left
    let bc1 = BoundaryCondition("dirichlet".to_string(), dvector![300.]);
    // right
    let bc2 = BoundaryCondition("neumann".to_string(), dvector![0.]);
    // front, back, top, bottom
    let bc3 = BoundaryCondition("dirichlet".to_string(), dvector![0.]);

    let n_grid = (50, 1, 20);

    let mut block: problem::Problem<Properties> = problem::Problem::new();
    // bottom left
    let boundary_conditions = [
        Some(bc1.clone()),
        Some(bc2.clone()),
        Some(bc2.clone()),
        None,
        Some(bc2.clone()),
        None,
    ];
    block.add_zone(
        (0., 0., 0.),
        (0.5, 1., 0.5),
        n_grid,
        boundary_conditions,
        vec![heat_diffusion_equations.clone()],
        copper,
    );
    // bottom right
    let boundary_conditions = [
        None,
        Some(bc2.clone()),
        Some(bc2.clone()),
        Some(bc3.clone()),
        Some(bc2.clone()),
        None,
    ];
    block.add_zone(
        (0.5, 0., 0.),
        (1., 1., 0.5),
        n_grid,
        boundary_conditions,
        vec![heat_diffusion_equations.clone()],
        copper,
    );
    // top left
    let boundary_conditions = [
        Some(bc1.clone()),
        Some(bc2.clone()),
        None,
        None,
        Some(bc2.clone()),
        Some(bc2.clone()),
    ];
    block.add_zone(
        (0., 0., 0.5),
        (0.5, 1., 1.),
        n_grid,
        boundary_conditions,
        vec![heat_diffusion_equations.clone()],
        copper,
    );
    // top right
    let boundary_conditions = [
        None,
        Some(bc2.clone()),
        None,
        Some(bc3.clone()),
        Some(bc2.clone()),
        Some(bc2.clone()),
    ];
    block.add_zone(
        (0.5, 0., 0.5),
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
pub fn get_T(p: PropertySlice) -> f64 {
    *p.T
}

#[allow(non_snake_case)]
pub fn update_T_user(p: PropertySliceMut, T: f64) {
    *p.T = T;
}

pub fn heat_diffusion_user(problem: &problem::Problem<Properties>, i_p: usize) -> f64 {
    // TODO : deal with boundary conditions
    let p = &problem
        .mesh
        .cells
        .get_cell(i_p, problem.cell_properties.slice(i_p));
    let mut flux = 0.;
    for neighbor in problem.neighbors(i_p) {
        match neighbor {
            problem::NeighborReturn::Internal(i_n, i_f, distance) => {
                let n = &problem
                    .mesh
                    .cells
                    .get_cell(i_n, problem.cell_properties.slice(i_n));
                let f = &problem.mesh.faces[i_f];
                flux -= p.properties.lambda * (p.properties.T - n.properties.T) / distance * f.area;
            }
            problem::NeighborReturn::Boundary(name, i_f, distance, values) => {
                let f = &problem.mesh.faces[i_f];
                match name {
                    "dirichlet" => {
                        flux -=
                            p.properties.lambda * (p.properties.T - values[0]) / distance * f.area;
                    }
                    "neumann" => {
                        flux -= values[0];
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
