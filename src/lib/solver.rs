use crate::computational_cell::HasProperties;
use crate::problem::Problem;
use indicatif::ProgressBar;
use nalgebra::DVector;
use std::time::Instant;

pub struct SolverOutputArgs {
    pub frequency: f64,
    pub filename: String,
    pub print_csv: bool,
    pub print_vtk: bool,
}

pub fn solve<T: HasProperties>(
    initial_problem: &Problem<T>,
    t_end: f64,
    h0: f64,
    solver: fn(&mut Problem<T>, f64, &mut Vec<DVector<f64>>, f64) -> Vec<DVector<f64>>,
    output_args: SolverOutputArgs,
) -> Problem<T> {
    let mut problem = initial_problem.clone();
    let mut t = 0.;
    let mut it = 0;

    // x contains the initial solution
    let mut x = initial_problem.get_state();
    let mut last_plotted_time = 0.;
    // we initialize the output directory if needed
    if output_args.frequency > 0. && (output_args.print_csv || output_args.print_vtk) {
        std::fs::remove_dir_all(format!("{}", &output_args.filename)).ok();
        std::fs::create_dir_all(format!("{}", &output_args.filename))
            .expect("Not able to create output directory");
        // export initial state
        if output_args.print_vtk {
            problem.export_to_vtk(format!("{}/step0.vtk", &output_args.filename).as_str());
        }
    }
    // make a nice progress bar
    println!("👩‍🔬 Solving ...");
    let pb = ProgressBar::new((t_end / h0) as u64);
    pb.set_style(
        indicatif::ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}")
            .unwrap()
            .progress_chars("##-"),
    );
    // main loop
    let start = Instant::now();
    while t < t_end {
        // update progress bar
        pb.inc(1);
        let mut cons_msg = String::new();
        for el in problem.sum_conservation_variables() {
            for val in el.iter() {
                cons_msg.push_str(&format!("  {:>+7.5} ", val));
            }
        }
        pb.set_message(format!("t = {:.2} s, convervation:", t) + cons_msg.as_str());
        // solve
        let y = (solver)(&mut problem, t, &mut x, h0);
        if stop_criterial(&y, t) {
            break;
        };
        t += h0;
        it += 1;
        // output
        if output_args.frequency > 0. && t > last_plotted_time + output_args.frequency {
            if output_args.print_vtk {
                problem.export_to_vtk(format!("{}/step{}.vtk", &output_args.filename, it).as_str());
            }
            if output_args.print_csv {
                let decimals = -f64::log10(h0).ceil() as usize;
                problem
                    .export_csv(
                        format!("{}/time{:.decimals$}.csv", &output_args.filename, t).as_str(),
                    )
                    .expect("Not able to export csv");
            }
            last_plotted_time = t;
        }
    }
    pb.finish();
    let duration = start.elapsed();
    println!(
        "Elapsed time : {} ms. Iterations : {}",
        duration.as_millis(),
        it
    );
    if output_args.frequency == 0. {
        println!("Exporting results ...");
        problem.export_to_vtk(format!("{}.vtk", &output_args.filename).as_str());
    }
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
            // final derivative
            y[i][j] = (k1[i][j] + 2. * k2[i][j] + 2. * k3[i][j] + k4[i][j]) / 6.;
            // final solution
            x[i][j] += h * y[i][j];
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
