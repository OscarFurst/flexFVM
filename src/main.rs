use flex_fvm::mesh::*;
use nalgebra::DVector;
use std::time::{Instant};

fn func(_t:&f64, x:&DVector<f64>, mesh:&mut Mesh) -> DVector<f64> {
    mesh.update(x);
    let mut y = DVector::from_element(x.len(), 0.);
    for ii in 0..x.len() {
        let (i,_j,k) = mesh.ijk(&ii);
        let radial;
        let axial;
        // TODO: improve behaviour in cases where some of the following cells is outside the mesh
        let p = mesh.p(&ii);
        let e = mesh.e(&ii);
        let w = mesh.w(&ii);
        let n = mesh.n(&ii);
        let s = mesh.s(&ii);
        if i == 0 {
            axial = p.rho * p.cp / p.dx * (
                (e.T - p.T) / p.de * (e.lambda + p.lambda) / 2. - 
                (p.T - 100.) / p.dw * (p.lambda + p.lambda) / 2. );
        } else if i == mesh.nx-1 {
            axial = p.rho * p.cp / p.dx * (
                (0. - p.T) / p.de * (p.lambda + p.lambda) / 2. - 
                (p.T - w.T) / p.dw * (p.lambda + w.lambda) / 2. );
        } else {
            axial = p.rho * p.cp / p.dx * (
                (e.T - p.T) / p.de * (e.lambda + p.lambda) / 2. - 
                (p.T - w.T) / p.dw * (p.lambda + w.lambda) / 2. );
        }
        if k == 0 {
            radial = p.rho * p.cp / p.dy * (
                (n.T - p.T) / p.dn * (n.lambda + p.lambda) / 2. - 
                (0.) / 2. );
        } else if k == mesh.nz-1 {
            radial = p.rho * p.cp / p.dy * (
                (0.) / 2. - 
                (p.T - s.T) / p.ds * (p.lambda + s.lambda) / 2. );
        } else {
            radial = p.rho * p.cp / p.dy * (
                (n.T - p.T) / p.dn * (n.lambda + p.lambda) / 2. - 
                (p.T - s.T) / p.ds * (p.lambda + s.lambda) / 2. );
        }
        y[ii] = radial + axial;
    }
    y
}

fn main() {
    let (nx, ny, nz) = (100, 1, 40);
    let mut mesh = Mesh::generate((0.,0.,0.), (1.,1.,1.), nx, ny, nz);

    let t_end:f64 = 100.;
    let mut t = 0.;
    let mut it = 0;
    let h:f64 = 0.01;
    let mut x = DVector::from_element(nx*ny*nz, 0.);
    let start = Instant::now();
    while t < t_end {
        //println!("Time = {} s", t);
        mesh.update(&x);
        let y = func(&t, &x, &mut mesh);
        x = x + h * &y;
        if y.abs().max() < 0.0000001 {
            println!("Steady state reached at t = {} s", t);
            break;
        }
        t = t+h;
        it += 1;
        //println!("{:?}", mesh);
    }
    let duration = start.elapsed();
    println!("Elapsed time : {} ms. Iterations : {}", duration.as_millis(), it);
    println!("{:?}", mesh.grid)
}
