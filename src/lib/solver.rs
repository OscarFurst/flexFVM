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

pub trait SolverStrategy {
    fn new<T: HasProperties>(problem: &Problem<T>) -> Self;

    fn solve<T: HasProperties>(
        &mut self,
        problem: &mut Problem<T>,
        t: f64,
        x: &mut DVector<f64>,
        h: f64,
    ) -> bool;
}

pub fn solve<T: HasProperties, S: SolverStrategy>(
    initial_problem: &Problem<T>,
    t_end: f64,
    h0: f64,
    mut solver: S,
    output_args: SolverOutputArgs,
) -> Problem<T> {
    let mut problem = initial_problem.clone();
    problem.set_equation_positions();
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
        pb.set_message(format!("t = {:.2} s", t));
        // solve
        if solver.solve(&mut problem, t, &mut x, h0) {
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

/// Euler explicit first order solver
pub struct EulerExplicitSolver {
    // the vector is allocated once and reused to improve performance
    y: DVector<f64>,
}

impl SolverStrategy for EulerExplicitSolver {
    fn new<T: HasProperties>(problem: &Problem<T>) -> Self {
        EulerExplicitSolver {
            y: DVector::zeros(problem.n_eq()),
        }
    }

    fn solve<T: HasProperties>(
        &mut self,
        problem: &mut Problem<T>,
        t: f64,
        x: &mut DVector<f64>,
        h: f64,
    ) -> bool {
        problem.update(x);
        problem.setup_equations_in_place(t, &mut self.y);
        for i in 0..x.len() {
            x[i] += h * self.y[i]; // loop is faster than DVector product last time I checked
        }
        stop_criterial(&self.y, t)
    }
}

/// Runge-Kutta fourth order explicit solver
pub struct RK4Solver {
    // all these vectors are allocated once and reused to improve performance
    k1: DVector<f64>,
    k2: DVector<f64>,
    k3: DVector<f64>,
    k4: DVector<f64>,
    x_alt: DVector<f64>,
    y: DVector<f64>,
}

impl SolverStrategy for RK4Solver {
    fn new<T: HasProperties>(problem: &Problem<T>) -> Self {
        RK4Solver {
            k1: DVector::zeros(problem.n_eq()),
            k2: DVector::zeros(problem.n_eq()),
            k3: DVector::zeros(problem.n_eq()),
            k4: DVector::zeros(problem.n_eq()),
            x_alt: DVector::zeros(problem.n_eq()),
            y: DVector::zeros(problem.n_eq()),
        }
    }

    fn solve<T: HasProperties>(
        &mut self,
        problem: &mut Problem<T>,
        t: f64,
        x: &mut DVector<f64>,
        h: f64,
    ) -> bool {
        problem.update(x);
        problem.setup_equations_in_place(t, &mut self.k1);
        self.x_alt.copy_from(x);
        for i in 0..x.len() {
            self.x_alt[i] += 0.5 * h * self.k1[i];
        }
        problem.update(&self.x_alt);
        problem.setup_equations_in_place(t + 0.5 * h, &mut self.k2);
        self.x_alt.copy_from(x);
        for i in 0..x.len() {
            self.x_alt[i] += 0.5 * h * self.k2[i];
        }
        problem.update(&self.x_alt);
        problem.setup_equations_in_place(t + 0.5 * h, &mut self.k3);
        self.x_alt.copy_from(x);
        for i in 0..x.len() {
            self.x_alt[i] += h * self.k3[i];
        }
        problem.update(&self.x_alt);
        problem.setup_equations_in_place(t + h, &mut self.k4);
        for i in 0..x.len() {
            // final derivative
            self.y[i] = (self.k1[i] + 2. * self.k2[i] + 2. * self.k3[i] + self.k4[i]) / 6.;
            // final solution
            x[i] += h * self.y[i];
        }
        stop_criterial(&self.y, t)
    }
}

fn stop_criterial(y: &DVector<f64>, t: f64) -> bool {
    if y.iter().any(|x| x.is_nan()) {
        println!("NaN found in solution vector");
        return true;
    } else {
        if y.abs().max() < 0.0000001 {
            println!("Steady state reached at t = {} s", t);
            return true;
        }
    }

    false
}
