use debug_cell::Ref;
use flex_fvm::{
    cell,
    equation::{BoundaryCondition, Equation},
    materials, mesh, solver,
};

fn main() {
    let initial_mesh = create_mesh();
    let solution = solver::solve(&initial_mesh, 0.001, 0.000001);
    println!("{:?}", solution)
}

pub fn create_mesh() -> mesh::Mesh {
    let bc1 = BoundaryCondition(heat_dirichlet, 300.);
    let bc2 = BoundaryCondition(heat_neumann, 0.);
    let bc3 = BoundaryCondition(heat_dirichlet, 0.);
    let boundary_conditions: [Vec<BoundaryCondition>; 6] = [
        vec![bc1],
        vec![bc2],
        vec![bc2],
        vec![bc3],
        vec![bc2],
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
        materials::copper(),
        materials::copper(),
    );
    block.add_layer(
        (0.5, 0., 0.),
        (1., 1., 0.5),
        n_grid,
        boundary_conditions.clone(),
        vec![Equation(update_T, heat_diffusion)],
        materials::copper(),
        materials::copper(),
    );
    block.add_layer(
        (0., 0., 0.5),
        (0.5, 1., 1.),
        n_grid,
        boundary_conditions.clone(),
        vec![Equation(update_T, heat_diffusion)],
        materials::copper(),
        materials::copper(),
    );
    block.add_layer(
        (0.5, 0., 0.5),
        (1., 1., 1.),
        n_grid,
        boundary_conditions,
        vec![Equation(update_T, heat_diffusion)],
        materials::copper(),
        materials::copper(),
    );
    block
}

///////////////////////////////
/// Heat diffusion equation ///
///////////////////////////////
#[allow(non_snake_case)]
pub fn update_T(p: &mut cell::Properties, T: f64) {
    p.T = T;
}

pub fn heat_convection_diffusion(p: &Ref<cell::Cell>) -> f64 {
    let mut diffusive_flux = 0.;
    let mut convective_flux = 0.;
    for f in p.faces.iter() {
        let n = f.neighbour();
        diffusive_flux -= p.properties.D * (p.properties.T - n.properties.T) / f.distance * f.area;
        convective_flux -=
            ((p.properties.h() + n.properties.h()) * (p.properties.u + n.properties.u) / 4.)
                .dot(&f.n)
                * f.area;
    }
    (diffusive_flux + convective_flux) / p.volume / p.properties.rho / p.properties.cp
}

pub fn heat_diffusion(p: &Ref<cell::Cell>) -> f64 {
    let mut flux = 0.;
    for f in p.faces.iter() {
        let n = f.neighbour();
        flux -= p.properties.lambda * (p.properties.T - n.properties.T) / f.distance * f.area;
    }
    flux / p.volume / p.properties.rho / p.properties.cp
}

#[allow(non_snake_case, unused_variables)]
pub fn heat_dirichlet(p: &mut cell::Properties, faces: &[cell::Face], T: f64) {
    p.T = T;
}

#[allow(non_snake_case)]
pub fn heat_neumann(p: &mut cell::Properties, faces: &[cell::Face], dTdx: f64) {
    let n = faces[0].neighbour();
    let mut boundary_value = n.properties.T;
    boundary_value += dTdx * faces[0].distance;
    p.T = boundary_value;
}
