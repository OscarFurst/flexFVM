use crate::{equation::*, materials, mesh};
use nalgebra::Vector3;

pub fn water_channel_example() -> mesh::Mesh {
    let bc1 = BoundaryCondition(heat_dirichlet, 300.);
    let bc2 = BoundaryCondition(heat_neumann, 0.);
    let bc3 = BoundaryCondition(heat_neumann, 0.);
    let boundary_conditions: [Vec<BoundaryCondition>; 6] = [
        vec![bc1],
        vec![bc2],
        vec![bc2],
        vec![bc3],
        vec![bc2],
        vec![bc2],
    ];
    let mut water = materials::water();
    water.u = Vector3::new(1., 0., 0.);
    let n_grid = (9, 1, 1);
    let mut block = mesh::Mesh::new();
    block.add_layer(
        (0., 0., 0.),
        (5., 1., 1.),
        n_grid,
        boundary_conditions.clone(),
        vec![heat_convection_diffusion],
        water,
        water,
    );
    block
}

pub fn copper_example() -> mesh::Mesh {
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
        vec![heat_diffusion],
        materials::copper(),
        materials::copper(),
    );
    block.add_layer(
        (0.5, 0., 0.),
        (1., 1., 0.5),
        n_grid,
        boundary_conditions.clone(),
        vec![heat_diffusion],
        materials::copper(),
        materials::copper(),
    );
    block.add_layer(
        (0., 0., 0.5),
        (0.5, 1., 1.),
        n_grid,
        boundary_conditions.clone(),
        vec![heat_diffusion],
        materials::copper(),
        materials::copper(),
    );
    block.add_layer(
        (0.5, 0., 0.5),
        (1., 1., 1.),
        n_grid,
        boundary_conditions,
        vec![heat_diffusion],
        materials::copper(),
        materials::copper(),
    );
    block
}
