use crate::cell::HasProperties;
use crate::mesh;
use nalgebra::{dvector, DVector};
use std::time::Instant;

pub fn solve<T: HasProperties>(initial_mesh: &mesh::Mesh<T>, t_end: f64, h0: f64) -> mesh::Mesh<T> {
    let mut mesh = initial_mesh.clone();
    let mut t = 0.;
    let mut it = 0;

    // x contains the initial solution
    // TODO : this is pretty hard to define once the set of equations is complex
    let mut x = vec![dvector![0.]; mesh.vec_size()];
    let start = Instant::now();
    while t < t_end {
        mesh.update(&x);
        let y = mesh.setup_equations(&t, &x);
        for i in 0..x.len() {
            for j in 0..x[i].len() {
                x[i][j] += h0 * y[i][j]; // loop is faster than DVector product
            }
        }
        if stop_criterial(&y, t) {
            break;
        };
        t += h0;
        it += 1;
    }
    let duration = start.elapsed();
    println!(
        "Elapsed time : {} ms. Iterations : {}",
        duration.as_millis(),
        it
    );
    mesh
}

fn stop_criterial(y: &[DVector<f64>], t: f64) -> bool {
    if y.iter().any(|x| x.iter().any(|x| x.is_nan())) {
        println!("NaN found in solution vector");
        return true;
    } else {
        let residuals = DVector::from_vec(y.iter().flatten().map(|x| x.abs()).collect::<Vec<_>>());
        if residuals.abs().max() < 0.0000001 {
            println!("Steady state reached at t = {} s", t);
            return true;
        }
    }

    false
}
