// x,y,z coordinate system, x being left-to-right, y being back-to-front and z being bottom-to-top
use nalgebra::{Vector3,DVector};
use crate::cell::*;
use std::{fmt, fs::DirBuilder};

#[allow(dead_code,non_snake_case)]
pub struct Mesh {
    pub grid: Vec<Cell>,
    // amount os simulated cells along each dimension
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    pub n: usize, // simulated cells
    pub N: usize, // all cells, incl. boundary conditions
    pub id: usize,
    pub n_eq: usize,
    //pub r1: Vector3<f64>,
    //pub r2: Vector3<f64>,
}

#[allow(dead_code)]
impl Mesh {
    pub fn get_cells(&self) -> &[Cell] {
        &self.grid[self.ii((1,1,1))..=self.ii((self.nx,self.ny,self.nz))]
    }
    pub fn get_mut_cells(&mut self) -> &mut[Cell] {
        let i0 = self.ii((1,1,1));
        let i1 = self.ii((self.nx,self.ny,self.nz));
        &mut self.grid[i0..=i1]
    }

    pub fn ii(&self, (i,j,k):(usize,usize,usize)) -> usize {
        self.nx*self.ny*k + self.nx*j + i
    }

    pub fn ijk(&self, ii:usize) -> (usize, usize, usize) {
        let k = ii / (self.nx * self.ny);
        let j = ii % (self.nx * self.ny) / self.nx;
        let i = ii % (self.nx * self.ny) % self.nx;
        (i,j,k)
    }

    pub fn update(&mut self, x:&DVector<f64>) {
        for ii in 0..self.grid.len() {
            self.grid[ii].T = x[ii]; // TODO: needs some form of indexing once we have more equations
        }
    }

    #[allow(non_snake_case)]
    pub fn generate(r1:(f64,f64,f64), r2:(f64,f64,f64), nx:usize, ny:usize, nz:usize, id: usize, boundary_conditions:Vec<BoundaryDirection>) -> Mesh{
        // nx, ny, nz is the size of the mesh without the boundary cells
        // NX, NY, NZ adds a layer of boundary cells around the domain
        let (NX,NY,NZ) = (nx+2,ny+2,nz+2);
        let (dx, dy, dz) = ((r2.0-r1.0)/nx as f64, (r2.1-r1.1)/ny as f64, (r2.2-r1.2)/nz as f64);
        let mut mesh = Mesh{grid: vec![Cell{dx, dy, dz, layer_id: id, ..Default::default()}; NX*NY*NZ], nx, ny, nz, n:nx*ny*nz, N:NX*NY*NZ, id, n_eq: 1};
        let origin = Vector3::new(r1.0,r1.1,r1.2);
        for i in 0..NX {
            for j in 0..NY {
                for k in 0..NZ {
                    let index: usize = mesh.ii((i,j,k));
                    mesh.grid[index].id = index;
                    mesh.grid[index].r1 = origin + Vector3::new((i as f64 -1.)* dx, (j as f64 -1.) * dy, (k as f64 -1.) * dz);
                    mesh.grid[index].r2 = origin + Vector3::new(i as f64 * dx, j as f64 * dy, k as f64 * dz);
                    mesh.grid[index].r = (mesh.grid[index].r1 + mesh.grid[index].r2)/2.; // r meeds to be vector
                    if i == 0 {
                        for boundary in &boundary_conditions {
                            if let BoundaryDirection::West(bc) = *boundary {
                                if let BoundaryCondition::Dirichlet(value, node) = bc {
                                    mesh.grid[index].variant = CellType::Boundary(BoundaryCondition::Dirichlet(value, Some((id, mesh.ii((i+1,j,k))))));
                                } else if let BoundaryCondition::Neumann(value, node) = bc {
                                    mesh.grid[index].variant = CellType::Boundary(BoundaryCondition::Neumann(value, Some((id, mesh.ii((i+1,j,k))))));
                                }
                            }
                        }
                    } else if i == NX-1 {
                        for boundary in &boundary_conditions {
                            if let BoundaryDirection::East(bc) = *boundary {
                                if let BoundaryCondition::Dirichlet(value, node) = bc {
                                    mesh.grid[index].variant = CellType::Boundary(BoundaryCondition::Dirichlet(value, Some((id, mesh.ii((i-1,j,k))))));
                                } else if let BoundaryCondition::Neumann(value, node) = bc {
                                    mesh.grid[index].variant = CellType::Boundary(BoundaryCondition::Neumann(value, Some((id, mesh.ii((i-1,j,k))))));
                                }
                            }
                        }
                    } else if j == 0 {
                        for boundary in &boundary_conditions {
                            if let BoundaryDirection::Close(bc) = *boundary {
                                if let BoundaryCondition::Dirichlet(value, node) = bc {
                                    mesh.grid[index].variant = CellType::Boundary(BoundaryCondition::Dirichlet(value, Some((id, mesh.ii((i,j+1,k))))));
                                } else if let BoundaryCondition::Neumann(value, node) = bc {
                                    mesh.grid[index].variant = CellType::Boundary(BoundaryCondition::Neumann(value, Some((id, mesh.ii((i,j+1,k))))));
                                }
                            }
                        }
                    } else if j == NY-1 {
                        for boundary in &boundary_conditions {
                            if let BoundaryDirection::Far(bc) = *boundary {
                                if let BoundaryCondition::Dirichlet(value, node) = bc {
                                    mesh.grid[index].variant = CellType::Boundary(BoundaryCondition::Dirichlet(value, Some((id, mesh.ii((i,j-1,k))))));
                                } else if let BoundaryCondition::Neumann(value, node) = bc {
                                    mesh.grid[index].variant = CellType::Boundary(BoundaryCondition::Neumann(value, Some((id, mesh.ii((i,j-1,k))))));
                                }
                            }
                        }
                    } else if k == 0 {
                        for boundary in &boundary_conditions {
                            if let BoundaryDirection::South(bc) = *boundary {
                                if let BoundaryCondition::Dirichlet(value, node) = bc {
                                    mesh.grid[index].variant = CellType::Boundary(BoundaryCondition::Dirichlet(value, Some((id, mesh.ii((i,j,k+1))))));
                                } else if let BoundaryCondition::Neumann(value, node) = bc {
                                    mesh.grid[index].variant = CellType::Boundary(BoundaryCondition::Neumann(value, Some((id, mesh.ii((i,j,k+1))))));
                                }
                            }
                        }
                    } else if k == NZ-1 {
                        for boundary in &boundary_conditions {
                            if let BoundaryDirection::North(bc) = *boundary {
                                if let BoundaryCondition::Dirichlet(value, node) = bc {
                                    mesh.grid[index].variant = CellType::Boundary(BoundaryCondition::Dirichlet(value, Some((id, mesh.ii((i,j,k-1))))));
                                } else if let BoundaryCondition::Neumann(value, node) = bc {
                                    mesh.grid[index].variant = CellType::Boundary(BoundaryCondition::Neumann(value, Some((id, mesh.ii((i,j,k-1))))));
                                }
                            }
                        }
                    } else {
                        mesh.grid[index].east = (id, mesh.ii((i+1,j,k)));
                        mesh.grid[index].west = (id, mesh.ii((i-1,j,k)));
                        mesh.grid[index].far = (id, mesh.ii((i,j+1,k)));
                        mesh.grid[index].close = (id, mesh.ii((i,j-1,k)));
                        mesh.grid[index].north = (id, mesh.ii((i,j,k+1)));
                        mesh.grid[index].south = (id, mesh.ii((i,j,k-1)));
                    }
                }
            }
        }
        mesh
    }
}

impl fmt::Debug for Mesh {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "[ ")?;
        for i in 0..self.nx+2 {
            let ii = self.ii((i,0,0));
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
        let mesh = Mesh::generate((0.,0.,0.), (1.,1.,1.), nx, ny, nz, 1, Vec::new());
        assert_eq!(mesh.ii((0,0,0)), 0);
        assert_eq!(mesh.ii((1,0,0)), 1);
        assert_eq!(mesh.ii((0,1,0)), 5);
        assert_eq!(mesh.ii((0,0,1)), 20);
        assert_eq!(mesh.ii((4,3,2)), 59);
    }
    #[test]
    fn test_ijk() {
        let (nx, ny, nz) = (5, 4, 3);
        let mesh = Mesh::generate((0.,0.,0.), (1.,1.,1.), nx, ny, nz, 1, Vec::new());
        assert_eq!(mesh.ijk(nx*ny*nz-1), (4,3,2));
        assert_eq!(mesh.ijk(1), (1,0,0));
    }
    #[test]
    fn test_mesh() {
        let (nx, ny, nz) = (4, 4, 4);
        let mesh = Mesh::generate((0.,0.,0.), (1.,1.,1.), nx, ny, nz, 1, Vec::new());
        assert_eq!(mesh.grid[63].r2, Vector3::new(1.,1.,1.));
        assert_eq!(mesh.grid[63].r1, Vector3::new(0.75,0.75,0.75));
    }
}