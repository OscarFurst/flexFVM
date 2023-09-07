use debug_cell::Ref;
use flex_fvm::{
    cell::{self, HasProperties},
    equation::{BoundaryCondition, Equation},
    mesh, solver,
};
use std::fmt;

fn main() {
    let initial_mesh = create_mesh();
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
    const N_VARIABLES: usize = 1;
    fn update(&mut self) {
        // x contains the solution of the equations
        // the other variables may also be updated here
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
        ..Default::default()
    }
}

pub fn create_mesh() -> mesh::Mesh<Properties> {
    let bc1 = BoundaryCondition(heat_dirichlet, vec![300.]);
    let bc2 = BoundaryCondition(heat_neumann, vec![0.]);
    let bc3 = BoundaryCondition(heat_dirichlet, vec![0.]);
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
        vec![Equation(update_T, heat_diffusion)],
        copper(),
        copper(),
    );
    block.add_layer(
        (0.5, 0., 0.),
        (1., 1., 0.5),
        n_grid,
        boundary_conditions.clone(),
        vec![Equation(update_T, heat_diffusion)],
        copper(),
        copper(),
    );
    block.add_layer(
        (0., 0., 0.5),
        (0.5, 1., 1.),
        n_grid,
        boundary_conditions.clone(),
        vec![Equation(update_T, heat_diffusion)],
        copper(),
        copper(),
    );
    block.add_layer(
        (0.5, 0., 0.5),
        (1., 1., 1.),
        n_grid,
        boundary_conditions,
        vec![Equation(update_T, heat_diffusion)],
        copper(),
        copper(),
    );
    block
}

///////////////////////////////
/// Heat diffusion equation ///
///////////////////////////////
#[allow(non_snake_case)]
pub fn update_T(p: &mut Properties, T: &Vec<f64>) {
    p.T = T[0];
}

pub fn heat_diffusion(p: &Ref<cell::Cell<Properties>>) -> Vec<f64> {
    let mut flux = 0.;
    for f in p.faces.iter() {
        let n = f.neighbour();
        flux -= p.properties.lambda * (p.properties.T - n.properties.T) / f.distance * f.area;
    }
    vec![flux / p.volume / p.properties.rho / p.properties.cp]
}

#[allow(non_snake_case, unused_variables)]
pub fn heat_dirichlet(p: &mut Properties, faces: &[cell::Face<Properties>], T: &Vec<f64>) {
    p.T = T[0];
}

#[allow(non_snake_case)]
pub fn heat_neumann(p: &mut Properties, faces: &[cell::Face<Properties>], dTdx: &Vec<f64>) {
    let n = faces[0].neighbour();
    let mut boundary_value = n.properties.T;
    boundary_value += dTdx[0] * faces[0].distance;
    p.T = boundary_value;
}
