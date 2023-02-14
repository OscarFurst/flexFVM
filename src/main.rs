use hello::mesh::*;

fn func(_t:&f64, x:&Vec<f64>, mesh:&mut Mesh) -> Vec<f64> {
    mesh.update(x);
    let mut y = vec![0.; x.len()];
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
        } else if k == mesh.ny-1 {
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
    let (nx, ny, nz) = (10, 4, 4);
    let mut mesh = Mesh::generate((0.,0.,0.), (1.,1.,1.), nx, ny, nz);

    let t_end:f64 = 10000.;
    let mut t = 0.;
    let h:f64 = 5.0;
    let mut x = vec![0.;nx*ny*nz];
    while t < t_end {
        println!("Time = {} s", t);
        mesh.update(&x);
        let y = func(&t, &x, &mut mesh);
        for ii in 0..x.len() {
            x[ii] = x[ii] + h * y[ii];
        }
        if y.iter().copied().fold(f64::NAN, f64::max) - y.iter().copied().fold(f64::NAN, f64::min) < 0.0000001 {
            print!("Steady state reached at t = {} s", t);
            break;
        }
        t = t+h;
        //println!("{:?}", mesh);
    }
    println!("{:?}", mesh.grid)
}
