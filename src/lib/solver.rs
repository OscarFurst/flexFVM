use crate::computational_cell::HasProperties;
use crate::problem::Problem;
use indicatif::ProgressBar;
use nalgebra::DVector;
use std::time::Instant;

pub fn solve<T: HasProperties>(initial_problem: &Problem<T>, t_end: f64, h0: f64) -> Problem<T> {
    let mut problem = initial_problem.clone();
    let mut t = 0.;
    let mut it = 0;

    // x contains the initial solution
    let mut x = initial_problem.get_state();
    let start = Instant::now();
    println!("👩‍🔬 Solving ...");
    let pb = ProgressBar::new((t_end / h0) as u64);
    while t < t_end {
        pb.inc(1);
        problem.update(&x);
        let y = problem.setup_equations(&t, &x);
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
    pb.finish();
    let duration = start.elapsed();
    println!(
        "Elapsed time : {} ms. Iterations : {}",
        duration.as_millis(),
        it
    );
    problem
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
