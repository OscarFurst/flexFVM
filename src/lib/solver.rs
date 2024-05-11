use crate::computational_cell::HasProperties;
use crate::problem::Problem;
use indicatif::ProgressBar;
use nalgebra::DVector;
use std::time::Instant;

pub fn solve<T: HasProperties>(
    initial_problem: &Problem<T>,
    t_end: f64,
    h0: f64,
    solver: fn(&mut Problem<T>, f64, &mut Vec<DVector<f64>>, f64) -> Vec<DVector<f64>>,
) -> Problem<T> {
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
        let y = (solver)(&mut problem, t, &mut x, h0);
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
    println!("Exporting results ...");
    problem
}

pub fn euler_explicit<T: HasProperties>(
    problem: &mut Problem<T>,
    t: f64,
    x: &mut Vec<DVector<f64>>,
    h: f64,
) -> Vec<DVector<f64>> {
    problem.update(&x);
    let y = problem.setup_equations(&t, &x);
    for i in 0..x.len() {
        for j in 0..x[i].len() {
            x[i][j] += h * y[i][j]; // loop is faster than DVector product
        }
    }
    y
}

pub fn runge_kutta_4<T: HasProperties>(
    problem: &mut Problem<T>,
    t: f64,
    x: &mut Vec<DVector<f64>>,
    h: f64,
) -> Vec<DVector<f64>> {
    problem.update(&x);
    let k1 = problem.setup_equations(&t, &x);
    let mut x1 = x.clone();
    for i in 0..x.len() {
        for j in 0..x[i].len() {
            x1[i][j] += 0.5 * h * k1[i][j];
        }
    }
    problem.update(&x1);
    let k2 = problem.setup_equations(&(t + 0.5 * h), &x1);
    let mut x2 = x.clone();
    for i in 0..x.len() {
        for j in 0..x[i].len() {
            x2[i][j] += 0.5 * h * k2[i][j];
        }
    }
    problem.update(&x2);
    let k3 = problem.setup_equations(&(t + 0.5 * h), &x2);
    let mut x3 = x.clone();
    for i in 0..x.len() {
        for j in 0..x[i].len() {
            x3[i][j] += h * k3[i][j];
        }
    }
    problem.update(&x3);
    let k4 = problem.setup_equations(&(t + h), &x3);
    let mut y = x.clone();
    for i in 0..x.len() {
        for j in 0..x[i].len() {
            y[i][j] += h / 6. * (k1[i][j] + 2. * k2[i][j] + 2. * k3[i][j] + k4[i][j]);
        }
    }
    y
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
