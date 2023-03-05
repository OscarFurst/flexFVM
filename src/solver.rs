use crate::mesh;
use nalgebra::DVector;
use std::time::Instant;

pub fn solve(initial_mesh: &mesh::Mesh, t_end: f64, h0: f64) -> mesh::Mesh {
    let mut mesh = initial_mesh.clone();
    let mut t = 0.;
    let mut it = 0;

    let mut x = DVector::from_element(mesh.vec_size(), 0.);
    let start = Instant::now();
    while t < t_end {
        mesh.update(&x);
        let y = mesh.setup_equations(&t, &x);
        x += h0 * &y;
        if y.abs().max() < 0.0000001 {
            println!("Steady state reached at t = {} s", t);
            break;
        }
        t += h0;
        it += 1;
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
    mesh
}
