use flex_fvm::solver;
use flex_fvm::{equation, materials::copper, mesh};

fn main() {
    let bc1 = equation::BoundaryCondition(equation::heat_dirichlet, 300.);
    let bc2 = equation::BoundaryCondition(equation::heat_neumann, 0.);
    let bc3 = equation::BoundaryCondition(equation::heat_dirichlet, 0.);
    let boundary_conditions: [Vec<equation::BoundaryCondition>; 6] = [
        vec![bc1],
        vec![bc2],
        vec![bc2],
        vec![bc3],
        vec![bc2],
        vec![bc2],
    ];
    let (nx, ny, nz) = (50, 1, 20);
    let n_grid = (nx, ny, nz);
    let mut block: mesh::Mesh;
    block = mesh::Mesh::new();
    block.add_layer(
        (0., 0., 0.),
        (0.5, 1., 0.5),
        n_grid,
        boundary_conditions.clone(),
        vec![equation::heat_diffusion],
        copper(),
    );
    block.add_layer(
        (0.5, 0., 0.),
        (1., 1., 0.5),
        n_grid,
        boundary_conditions.clone(),
        vec![equation::heat_diffusion],
        copper(),
    );
    block.add_layer(
        (0., 0., 0.5),
        (0.5, 1., 1.),
        n_grid,
        boundary_conditions.clone(),
        vec![equation::heat_diffusion],
        copper(),
    );
    block.add_layer(
        (0.5, 0., 0.5),
        (1., 1., 1.),
        n_grid,
        boundary_conditions,
        vec![equation::heat_diffusion],
        copper(),
    );

    let solution = crate::solver::solve(&block, 0.001, 0.000001);
    println!("{:?}", solution)
}
