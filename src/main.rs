extern crate nalgebra as na;
use hello::mesh::*;
use na::DVector;
use std::time::{Instant};

fn func(_t:&f64, x:&DVector<f64>, mesh:&mut Mesh) -> DVector<f64> {
    mesh.update(x);
    let mut y = DVector::from_element(x.len(), 0.);
    for ii in 0..x.len() {
        let (i,_j,k) = mesh.ijk(&ii);
        let radial;
        let axial;
        if i == 0 {
            axial = mesh.p(&ii).rho * mesh.p(&ii).cp / mesh.p(&ii).dx * (
                (mesh.e(&ii).T - mesh.p(&ii).T) / (mesh.e(&ii).r - mesh.p(&ii).r)[0] * (mesh.e(&ii).lambda + mesh.p(&ii).lambda) / 2. - 
                (mesh.p(&ii).T - 100.         ) / (mesh.p(&ii).r - mesh.p(&ii).r1)[0]/2. * (mesh.p(&ii).lambda + mesh.p(&ii).lambda) / 2. );
        } else if i == mesh.nx-1 {
            axial = mesh.p(&ii).rho * mesh.p(&ii).cp / mesh.p(&ii).dx * (
                (0.            - mesh.p(&ii).T) / (mesh.p(&ii).r2 - mesh.p(&ii).r)[0]/2. * (mesh.p(&ii).lambda + mesh.p(&ii).lambda) / 2. - 
                (mesh.p(&ii).T - mesh.w(&ii).T) / (mesh.p(&ii).r - mesh.w(&ii).r)[0] * (mesh.p(&ii).lambda + mesh.w(&ii).lambda) / 2. );
        } else {
            axial = mesh.p(&ii).rho * mesh.p(&ii).cp / mesh.p(&ii).dx * (
                (mesh.e(&ii).T - mesh.p(&ii).T) / (mesh.e(&ii).r - mesh.p(&ii).r)[0] * (mesh.e(&ii).lambda + mesh.p(&ii).lambda) / 2. - 
                (mesh.p(&ii).T - mesh.w(&ii).T) / (mesh.p(&ii).r - mesh.w(&ii).r)[0] * (mesh.p(&ii).lambda + mesh.w(&ii).lambda) / 2. );
        }
        if k == 0 {
            radial = mesh.p(&ii).rho * mesh.p(&ii).cp / mesh.p(&ii).dy * (
                (mesh.n(&ii).T - mesh.p(&ii).T) / (mesh.n(&ii).r - mesh.p(&ii).r)[2] * (mesh.n(&ii).lambda + mesh.p(&ii).lambda) / 2. - 
                (0.) / 2. );
        } else if k == mesh.nz-1 {
            radial = mesh.p(&ii).rho * mesh.p(&ii).cp / mesh.p(&ii).dy * (
                (0.) / 2. - 
                (mesh.p(&ii).T - mesh.s(&ii).T) / (mesh.p(&ii).r - mesh.s(&ii).r)[2] * (mesh.p(&ii).lambda + mesh.s(&ii).lambda) / 2. );
        } else {
            radial = mesh.p(&ii).rho * mesh.p(&ii).cp / mesh.p(&ii).dy * (
                (mesh.n(&ii).T - mesh.p(&ii).T) / (mesh.n(&ii).r - mesh.p(&ii).r)[2] * (mesh.n(&ii).lambda + mesh.p(&ii).lambda) / 2. - 
                (mesh.p(&ii).T - mesh.s(&ii).T) / (mesh.p(&ii).r - mesh.s(&ii).r)[2] * (mesh.p(&ii).lambda + mesh.s(&ii).lambda) / 2. );
        }
        //println!("{} {}", radial, axial);
        y[ii] = radial + axial;
    }
    y
}

fn main() {
    let (nx, ny, nz) = (100, 1, 40);
    let mut mesh = Mesh::generate((0.,0.,0.), (1.,1.,1.), nx, ny, nz);

    let t_end:f64 = 10000.;
    let mut t = 0.;
    let mut it = 0;
    let h:f64 = 0.01;
    //let mut x = vec![0.;nx*ny*nz];
    let mut x = DVector::from_element(nx*ny*nz, 0.);
    let start = Instant::now();
    while t < t_end {
        println!("Time = {} s", t);
        mesh.update(&x);
        let y = func(&t, &x, &mut mesh);
        x = x + h * &y;
        if y.abs().max() < 0.0000001 {
            println!("Steady state reached at t = {} s", t);
            break;
        }
        t = t+h;
        it += 1;
        println!("{:?}", mesh);
    }
    let duration = start.elapsed();
    println!("Elapsed time : {} ms. Iterations : {}", duration.as_millis(), it);
    println!("{:?}", mesh.grid)
}
