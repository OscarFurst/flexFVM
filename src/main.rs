use nalgebra::DVector;
use std::io::{Write, Read};
use std::time::{Instant};
use std::fs;
use std::path::Path;
use flex_fvm::{solver,cell::BoundaryCondition};

fn main() {
    let boundary_conditions = [
        BoundaryCondition::Dirichlet(300.),
        BoundaryCondition::Neumann(0.),
        BoundaryCondition::Neumann(0.),
        BoundaryCondition::Dirichlet(0.),
        BoundaryCondition::Neumann(0.),
        BoundaryCondition::Neumann(0.)
        ];
    let (nx, ny, nz) = (100, 1, 40);
    let filename = "mesh.json";
    let mut block: solver::Domain;
    //fs::remove_file(filename).unwrap();
    if Path::new(filename).exists() {
        let mut file = fs::File::open(filename).unwrap();
        let mut serialized = String::new();
        file.read_to_string(&mut serialized).unwrap();
        block = serde_json::from_str(&serialized).unwrap();
    } else {
        block = solver::Domain::new();
        block.add_layer((0.,0.,0.), (1.,1.,1.), nx, ny, nz, boundary_conditions);
        // save mesh
        let serialized = serde_json::to_string(&block).unwrap();
        let mut file = fs::File::create("mesh.json").unwrap();
        file.write_all(serialized.as_bytes()).unwrap();
    }

    let t_end:f64 = 100.;
    let mut t = 0.;
    let mut it = 0;
    let h:f64 = 0.01;
    let mut x = DVector::from_element(nx*ny*nz, 0.);
    let start = Instant::now();
    while t < t_end {
        //println!("Time = {} s", t);
        block.update(&x);
        let y = block.setup_equations(&t, &x);
        x = x + h * &y;
        if y.abs().max() < 0.0000001 {
            println!("Steady state reached at t = {} s", t);
            break;
        }
        t = t+h;
        it += 1;
        // println!("{:?}", y);
        if y.iter().any(|x| x.is_nan()) {break}
    }
    let duration = start.elapsed();
    println!("Elapsed time : {} ms. Iterations : {}", duration.as_millis(), it);
    println!("{:?}", block)
}
