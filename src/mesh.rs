// x,y,z coordinate system, x being left-to-right, y being back-to-front and z being bottom-to-top

extern crate nalgebra as na;
use na::{Vector3,DVector};
use std::fmt;

#[allow(non_snake_case)]
#[derive(Clone)]
pub struct Cell {
    pub r: Vector3<f64>,
    pub r1: Vector3<f64>,
    pub r2: Vector3<f64>,
    pub dx: f64,
    pub dy: f64,
    pub dz: f64,
    pub T: f64,
    pub rho: f64,
    pub cp: f64,
    pub lambda: f64,
}

impl fmt::Debug for Cell {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:.2}", self.T)
    }
}

impl Default for Cell {
    fn default() -> Self {
        Cell {
            r: Vector3::new(0.,0.,0.),
            r1: Vector3::new(0.,0.,0.),
            r2: Vector3::new(0.,0.,0.),
            dx: 0.,
            dy: 0.,
            dz: 0.,
            T: 0.0,
            rho: 1.0,
            cp: 0.1,
            lambda: 0.01
        }
    }
}

#[allow(dead_code)]
pub struct Mesh {
    pub grid: Vec<Cell>,
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
}

#[allow(dead_code)]
impl Mesh {
    pub fn ii(&self, (i,j,k):(&usize,&usize,&usize)) -> usize {
        self.nx*self.ny*k + self.nx*j + i
    }

    pub fn ijk(&self, ii:&usize) -> (usize, usize, usize) {
        let k = ii / (self.nx * self.ny);
        let j = ii % (self.nx * self.ny) / self.nx;
        let i = ii % (self.nx * self.ny) % self.nx;
        (i,j,k)
    }

    pub fn p(&mut self, ii:&usize) -> &mut Cell {
        &mut self.grid[*ii]
    }

    pub fn n(&mut self, ii:&usize) -> &mut Cell {
        let (i,j,k) = self.ijk(ii);
        let iii = self.ii((&i,&j,&(&k+1)));
        &mut self.grid[iii]
    }

    pub fn s(&mut self, ii:&usize) -> &mut Cell {
        let (i,j,k) = self.ijk(ii);
        let iii = self.ii((&i,&j,&(&k-1)));
        &mut self.grid[iii]
    }

    pub fn e(&mut self, ii:&usize) -> &mut Cell {
        let (i,j,k) = self.ijk(ii);
        let iii = self.ii((&(&i+1),&j,&k));
        &mut self.grid[iii]
    }

    pub fn w(&mut self, ii:&usize) -> &mut Cell {
        let (i,j,k) = self.ijk(ii);
        let iii = self.ii((&(&i-1),&j,&k));
        &mut self.grid[iii]
    }

    pub fn update(&mut self, x:&DVector<f64>) {
        for ii in 0..self.grid.len() {
            self.grid[ii].T = x[ii]; // TODO: needs some form of indexing once we have more equations
        }
    }

    pub fn generate(r1:(f64,f64,f64), r2:(f64,f64,f64), nx:usize, ny:usize, nz:usize) -> Mesh{
        let mut mesh = Mesh{grid: vec![Cell{..Default::default()}; nx*ny*nz], nx, ny, nz};
        let origin = Vector3::new(r1.0,r1.1,r1.2);
        let (dx, dy, dz) = ((r2.0-r1.0)/nx as f64, (r2.1-r1.1)/ny as f64, (r2.2-r1.2)/nz as f64);
        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    let index: usize = mesh.ii((&i,&j,&k));
                    mesh.grid[index].dx = dx;
                    mesh.grid[index].dy = dy;
                    mesh.grid[index].dz = dz;
                    mesh.grid[index].r1 = origin + Vector3::new(i as f64 * dx, j as f64 * dy, k as f64 * dz);
                    mesh.grid[index].r2 = origin + Vector3::new((i+1) as f64 * dx, (j+1) as f64 * dy, (k+1) as f64 * dz);

                    mesh.grid[index].r = (mesh.grid[index].r1 + mesh.grid[index].r2)/2.; // r meeds to be vector
                }
            }
        }
        mesh
    }
}

impl fmt::Debug for Mesh {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "[ ")?;
        for i in 0..self.nx {
            let ii = self.ii((&i,&0,&0));
            write!(f, "{:.2} ", self.grid[ii].T)?;
        }
        write!(f, "]")
    }
}

#[cfg(test)]
mod tests {
    use crate::mesh::*;
    #[test]
    fn test_ii() {
        let (nx, ny, nz) = (5, 4, 3);
        let mesh = Mesh::generate((0.,0.,0.), (1.,1.,1.), nx, ny, nz);
        assert_eq!(mesh.ii((&0,&0,&0)), 0);
        assert_eq!(mesh.ii((&1,&0,&0)), 1);
        assert_eq!(mesh.ii((&0,&1,&0)), 5);
        assert_eq!(mesh.ii((&0,&0,&1)), 20);
        assert_eq!(mesh.ii((&4,&3,&2)), 59);
    }
    #[test]
    fn test_ijk() {
        let (nx, ny, nz) = (5, 4, 3);
        let mesh = Mesh::generate((0.,0.,0.), (1.,1.,1.), nx, ny, nz);
        assert_eq!(mesh.ijk(&(nx*ny*nz-1)), (4,3,2));
        assert_eq!(mesh.ijk(&(1)), (1,0,0));
    }
    #[test]
    fn test_mesh() {
        let (nx, ny, nz) = (4, 4, 4);
        let mesh = Mesh::generate((0.,0.,0.), (1.,1.,1.), nx, ny, nz);
        assert_eq!(mesh.grid[63].r2, Vector3::new(1.,1.,1.));
        assert_eq!(mesh.grid[63].r1, Vector3::new(0.75,0.75,0.75));
    }
}