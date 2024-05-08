use flex_fvm::{
    computational_cell::{self, ComputationalCell, HasProperties},
    equation::*,
    generate_equation_structure, mesh, problem, solver,
};
use nalgebra::{dvector, DVector};
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
    let solution = solver::solve(&initial_mesh, 0.001, 0.000001);
    println!("{:?}", solution)
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
    fn update(&mut self) {
        // here we can update derived variables if needed
    }
}

///////////////////
/// Copper data ///
///////////////////

pub fn copper() -> Properties {
    Properties {
        T: 0.,
        rho: 9000.,
        cp: 0.001,
        lambda: 300.,
        //..Default::default()
    }
}

#[allow(non_snake_case)]
pub fn create_mesh() -> problem::Problem<Properties> {
    // todo : add a macro that wraps all that into an "Equation" struct
    let heat_diffusion_equations =
        generate_equation_structure!(get_T, update_T_user, heat_diffusion_user, 1);
    // define the boundary conditions
    // left
    let bc1 = vec![BoundaryCondition(heat_dirichlet, dvector![300.])];
    // right
    let bc2 = vec![BoundaryCondition(heat_neumann, dvector![0.])];
    // front, back, top, bottom
    let bc3 = vec![BoundaryCondition(heat_dirichlet, dvector![0.])];

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
        vec![heat_diffusion_equations],
        copper(),
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
        boundary_conditions.clone(),
        vec![heat_diffusion_equations],
        copper(),
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
        boundary_conditions.clone(),
        vec![heat_diffusion_equations],
        copper(),
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
        copper(),
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

pub fn heat_diffusion_user(m: &Vec<ComputationalCell<Properties>>, i_p: &usize) -> f64 {
    // TODO : deal with boundary conditions
    let p = m[*i_p];
    let mut flux = 0.;
    for f in m[*i_p].faces.iter() {
        if let Some(n_n) = f.neighbour() {
            let n = m[n_n];
            flux -= p.properties.lambda * (p.properties.T - n.properties.T) / f.distance * f.area;
        }
    }
    flux / p.volume / p.properties.rho / p.properties.cp
}

#[allow(non_snake_case, unused_variables)]
pub fn heat_dirichlet(p: &mut Properties, face: &computational_cell::Face, T: &DVector<f64>) {
    // TODO
}

#[allow(non_snake_case)]
pub fn heat_neumann(p: &mut Properties, face: &computational_cell::Face, dTdx: &DVector<f64>) {
    // TODO
}
