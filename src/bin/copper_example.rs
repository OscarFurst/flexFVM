use flex_fvm::{
    cell::{self, HasProperties},
    equation::*,
    generate_equation_structure, mesh, solver,
};
use nalgebra::{dvector, DVector};
use std::cell::Ref;
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
        T: 300.,
        rho: 9000.,
        cp: 0.001,
        lambda: 300.,
        //..Default::default()
    }
}

#[allow(non_snake_case)]
pub fn create_mesh() -> mesh::Mesh<Properties> {
    // todo : add a macro that wraps all that into an "Equation" struct
    let heat_diffusion_equations =
        generate_equation_structure!(update_T_user, heat_diffusion_user, f64);
    // define the boundary conditions
    let bc1 = BoundaryCondition(heat_dirichlet, dvector![300.]);
    let bc2 = BoundaryCondition(heat_neumann, dvector![0.]);
    let bc3 = BoundaryCondition(heat_dirichlet, dvector![0.]);
    let boundary_conditions: [Vec<BoundaryCondition<Properties>>; 6] = [
        vec![bc1],
        vec![bc2.clone()],
        vec![bc2.clone()],
        vec![bc3],
        vec![bc2.clone()],
        vec![bc2],
    ];
    let n_grid = (50, 1, 20);
    let mut block = mesh::Mesh::new();
    block.add_layer(
        (0., 0., 0.),
        (0.5, 1., 0.5),
        n_grid,
        boundary_conditions.clone(),
        vec![heat_diffusion_equations],
        copper(),
        copper(),
    );
    block.add_layer(
        (0.5, 0., 0.),
        (1., 1., 0.5),
        n_grid,
        boundary_conditions.clone(),
        vec![heat_diffusion_equations],
        copper(),
        copper(),
    );
    block.add_layer(
        (0., 0., 0.5),
        (0.5, 1., 1.),
        n_grid,
        boundary_conditions.clone(),
        vec![heat_diffusion_equations],
        copper(),
        copper(),
    );
    block.add_layer(
        (0.5, 0., 0.5),
        (1., 1., 1.),
        n_grid,
        boundary_conditions,
        vec![heat_diffusion_equations],
        copper(),
        copper(),
    );
    block
}

///////////////////////////////
/// Heat diffusion equation ///
///////////////////////////////
#[allow(non_snake_case)]
pub fn update_T_user(p: &mut Properties, T: f64) {
    p.T = T;
}

pub fn heat_diffusion_user(p: &Ref<cell::Cell<Properties>>) -> f64 {
    let mut flux = 0.;
    for f in p.faces.iter() {
        let n = f.neighbour();
        flux -= p.properties.lambda * (p.properties.T - n.properties.T) / f.distance * f.area;
    }
    flux / p.volume / p.properties.rho / p.properties.cp
}

#[allow(non_snake_case, unused_variables)]
pub fn heat_dirichlet(p: &mut Properties, faces: &[cell::Face<Properties>], T: &DVector<f64>) {
    p.T = T[0];
}

#[allow(non_snake_case)]
pub fn heat_neumann(p: &mut Properties, faces: &[cell::Face<Properties>], dTdx: &DVector<f64>) {
    let n = faces[0].neighbour();
    let mut boundary_value = n.properties.T;
    boundary_value += dTdx[0] * faces[0].distance;
    p.T = boundary_value;
}
