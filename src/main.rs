use debug_cell::RefMut;
use flex_fvm::{cell::Cell, equation, solver};
use nalgebra::DVector;
use std::time::Instant;

fn main() {
    let bc1: (fn(&mut RefMut<Cell>, f64), f64) = (equation::heat_diffusion_dirichlet, 300.);
    let bc2: (fn(&mut RefMut<Cell>, f64), f64) = (equation::heat_diffusion_neumann, 0.);
    let bc3: (fn(&mut RefMut<Cell>, f64), f64) = (equation::heat_diffusion_dirichlet, 0.);
    let boundary_conditions: [Vec<(fn(&mut RefMut<Cell>, f64), f64)>; 6] = [
        vec![bc1],
        vec![bc2],
        vec![bc2],
        vec![bc3],
        vec![bc2],
        vec![bc2],
    ];
    let (nx, ny, nz) = (50, 1, 20);
    let mut block: solver::Domain;
    block = solver::Domain::new();
    block.add_layer(
        (0., 0., 0.),
        (0.5, 1., 0.5),
        nx,
        ny,
        nz,
        boundary_conditions.clone(),
        vec![equation::heat_diffusion],
    );
    block.add_layer(
        (0.5, 0., 0.),
        (1., 1., 0.5),
        nx,
        ny,
        nz,
        boundary_conditions.clone(),
        vec![equation::heat_diffusion],
    );
    block.add_layer(
        (0., 0., 0.5),
        (0.5, 1., 1.),
        nx,
        ny,
        nz,
        boundary_conditions.clone(),
        vec![equation::heat_diffusion],
    );
    block.add_layer(
        (0.5, 0., 0.5),
        (1., 1., 1.),
        nx,
        ny,
        nz,
        boundary_conditions.clone(),
        vec![equation::heat_diffusion],
    );

    let t_end: f64 = 0.001;
    let mut t = 0.;
    let mut it = 0;
    let h: f64 = 0.000001;
    let mut x = DVector::from_element(block.vec_size(), 0.);
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
