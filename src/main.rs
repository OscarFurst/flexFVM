use flex_fvm::{cell::BoundaryCondition, solver};
use nalgebra::DVector;
use std::time::Instant;

fn main() {
    let boundary_conditions = [
        BoundaryCondition::Dirichlet(300.),
        BoundaryCondition::Neumann(0.),
        BoundaryCondition::Neumann(0.),
        BoundaryCondition::Dirichlet(0.),
        BoundaryCondition::Neumann(0.),
        BoundaryCondition::Neumann(0.),
    ];
    let (nx, ny, nz) = (100, 1, 40);
    let mut block: solver::Domain;
    block = solver::Domain::new();
    block.add_layer((0., 0., 0.), (1., 1., 1.), nx, ny, nz, &boundary_conditions);

    let t_end: f64 = 0.001;
    let mut t = 0.;
    let mut it = 0;
    let h: f64 = 0.000001;
    let mut x = DVector::from_element(nx * ny * nz, 0.);
    let start = Instant::now();
    while t < t_end {
        //println!("Time = {} s", t);
        block.update(&x);
        let y = block.setup_equations(&t, &x);
        x += h * &y;
        if y.abs().max() < 0.0000001 {
            println!("Steady state reached at t = {} s", t);
            break;
        }
        t += h;
        it += 1;
        // println!("{:?}", y);
        if y.iter().any(|x| x.is_nan()) {
            break;
        }
    }
    let duration = start.elapsed();
    println!(
        "Elapsed time : {} ms. Iterations : {}",
        duration.as_millis(),
        it
    );
    println!("{:?}", block)
}
