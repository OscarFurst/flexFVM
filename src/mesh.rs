// x,y,z coordinate system, x being left-to-right, y being back-to-front and z being bottom-to-top
use crate::cell::*;
use debug_cell::RefCell;
use nalgebra::{DVector, Point3, Vector3};
use std::{fmt, rc::Rc};

#[allow(non_snake_case)]
pub struct Mesh {
    pub cells: Vec<Rc<RefCell<Cell>>>,
    pub boundaries: [(Point3<f64>, Vector3<f64>, BoundaryCondition); 6],
    pub n_eq: usize,
}

#[allow(dead_code)]
impl Mesh {
    pub fn update(&mut self, x: &DVector<f64>) {
        for ii in 0..self.cells.len() {
            self.cells[ii].borrow_mut().T = x[ii]; // TODO: needs some form of indexing once we have more equations
        }
    }

    #[allow(non_snake_case)]
    pub fn generate(
        r1: (f64, f64, f64),
        r2: (f64, f64, f64),
        nx: usize,
        ny: usize,
        nz: usize,
        bc: &[BoundaryCondition; 6],
    ) -> Mesh {
        // generates a regular grid of rectangular cells
        // TODO: I could not think of an elegant way to pass boundary conditions
        // in the current state, the order in which they are passed are very important!
        let origin = Point3::new(r1.0, r1.1, r1.2);
        let ending = Point3::new(r2.0, r2.1, r2.2);
        let ex = (ending - origin)[0] / nx as f64 * Vector3::new(1., 0., 0.);
        let ey = (ending - origin)[1] / ny as f64 * Vector3::new(0., 1., 0.);
        let ez = (ending - origin)[2] / nz as f64 * Vector3::new(0., 0., 1.);
        let boundaries = [
            (origin, -ex, bc[0].clone()),
            (origin, -ey, bc[1].clone()),
            (origin, -ez, bc[2].clone()),
            (ending, ex, bc[3].clone()),
            (ending, ey, bc[4].clone()),
            (ending, ez, bc[5].clone()),
        ];
        let mut mesh = Mesh {
            cells: Vec::new(),
            boundaries: boundaries.clone(),
            n_eq: 1,
        };
        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    let mut cell = Cell {
                        variant: CellType::Internal,
                        ..Default::default()
                    };
                    cell.create_rectangle(
                        origin + ex * i as f64 + ey * j as f64 + ez * k as f64,
                        ex,
                        ey,
                        ez,
                    );
                    // check if at boundary. if so, create boundary condition
                    for (point, normal, bc) in &boundaries {
                        if let Some(face) = cell.face_on_plane(point, normal) {
                            // search for corresponding boundary condition
                            let boundary_cell =
                                Cell::boundary_from_cell(&cell, normal, face, bc.clone());
                            mesh.cells.push(Rc::new(RefCell::new(boundary_cell)));
                        }
                    }
                    // neighbors can be set from a specialized routine that uses the vertices of the cells. it is needed anyway to merge multiple layers
                    mesh.cells.push(Rc::new(RefCell::new(cell)));
                }
            }
        }
        mesh
    }

    pub fn link(self) -> Self {
        // compare cells two by two
        for cell in self.cells.iter() {
            for other_cell in self.cells.iter() {
                // exclude the possibility of them being the same cell
                let dr = other_cell.borrow().r - cell.borrow().r;
                if dr.norm() < f64::EPSILON {
                    continue;
                }
                // compare faces
                for face in cell
                    .borrow_mut()
                    .faces
                    .iter_mut()
                    .filter(|f| matches!(f.neighbor, None))
                {
                    for other_face in other_cell.borrow_mut().faces.iter_mut() {
                        // pair up faces that matched
                        if (face.center - other_face.center).norm() < f64::EPSILON {
                            face.neighbor = Some(Neighbor {
                                cell: Rc::clone(other_cell),
                                dr: dr,
                                distance: dr.norm(),
                            });
                            other_face.neighbor = Some(Neighbor {
                                cell: Rc::clone(cell),
                                dr: -dr,
                                distance: dr.norm(),
                            });
                        }
                    }
                }
            }
        }
        self
    }
}

impl fmt::Debug for Mesh {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "[ ")?;
        for cell in self
            .cells
            .iter()
            .map(|rc| rc.borrow())
            .filter(|c| matches!(c.variant, CellType::Internal))
        {
            write!(f, "{:.1} ", cell.T)?;
        }
        write!(f, "]")
    }
}

#[cfg(test)]
mod tests {
    use crate::{cell::BoundaryCondition, mesh::Mesh};
    #[test]
    fn test_generation() {
        let bc = [
            BoundaryCondition::Dirichlet(300.),
            BoundaryCondition::Neumann(0.),
            BoundaryCondition::Neumann(0.),
            BoundaryCondition::Dirichlet(0.),
            BoundaryCondition::Neumann(0.),
            BoundaryCondition::Neumann(0.),
        ];
        let (nx, ny, nz) = (2, 2, 2);
        let layer = Mesh::generate((0., 0., 0.), (1., 1., 1.), nx, ny, nz, &bc);
        for cell in layer.cells.iter().map(|rc| rc.borrow()) {
            println!("coord : {:?} ; type : {:?}", cell.r, cell.variant);
        }
        assert_eq!(layer.cells.len(), 32);
    }

    #[test]
    fn test_large_mesh() {
        let bc = [
            BoundaryCondition::Dirichlet(300.),
            BoundaryCondition::Neumann(0.),
            BoundaryCondition::Neumann(0.),
            BoundaryCondition::Dirichlet(0.),
            BoundaryCondition::Neumann(0.),
            BoundaryCondition::Neumann(0.),
        ];
        let (nx, ny, nz) = (100, 1, 1);
        let layer = Mesh::generate((0., 0., 0.), (1., 1., 1.), nx, ny, nz, &bc);
        for cell in layer.cells.iter().map(|rc| rc.borrow()) {
            println!("coord : {:?} ; type : {:?}", cell.r, cell.variant);
        }
        println!("Cell volume : {}", layer.cells[0].borrow().volume);
        assert_eq!(
            layer.cells.len(),
            (nx + 2) * (ny + 2) * (nz + 2) - 4 * nx - 4 * ny - 4 * nz - 8
        );
    }
}
