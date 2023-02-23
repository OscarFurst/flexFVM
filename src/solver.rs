use crate::cell;
use crate::cell::{BoundaryCondition::*, CellType::*};
use crate::mesh;
use nalgebra::DVector;
use std::fmt;
use std::rc::Rc;

pub struct Domain {
    layers: Vec<mesh::Mesh>,
}

impl Domain {
    pub fn new() -> Self {
        Domain { layers: Vec::new() }
    }

    pub fn add_layer(
        &mut self,
        r1: (f64, f64, f64),
        r2: (f64, f64, f64),
        nx: usize,
        ny: usize,
        nz: usize,
        bc: [cell::BoundaryCondition; 6],
    ) {
        self.layers
            .push(mesh::Mesh::generate(r1, r2, nx, ny, nz, bc));
        self.link_cells()
    }

    fn link_cells(&self) {
        // compare cells two by two
        for layer in self.layers.iter() {
            for cell in layer.cells.iter() {
                for other_layer in self.layers.iter() {
                    for other_cell in other_layer.cells.iter() {
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
                                    face.neighbor = Some(cell::Neighbor {
                                        cell: Rc::clone(other_cell),
                                        dr: dr,
                                        distance: dr.norm(),
                                    });
                                    other_face.neighbor = Some(cell::Neighbor {
                                        cell: Rc::clone(cell),
                                        dr: -dr,
                                        distance: dr.norm(),
                                    });
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    pub fn enforce_boundary_conditions(&mut self) {
        for layer in self.layers.iter() {
            for cell in layer.cells.iter() {
                let bc = cell.borrow().variant.clone();
                if let Boundary(Dirichlet(val)) = bc {
                    cell.borrow_mut().T = val;
                } else if let Boundary(Neumann(val)) = bc {
                    if cell.borrow().faces.len() != 1 {
                        panic!("boundary cells should only have a single neighbor");
                    }
                    let boundary_value = cell.borrow().faces[0]
                        .neighbor
                        .as_ref()
                        .unwrap()
                        .cell
                        .borrow()
                        .T; // gotta check if there is more than one face
                    cell.borrow_mut().T = boundary_value;
                }
            }
        }
    }

    pub fn setup_equations(&self, _t: &f64, x: &DVector<f64>) -> DVector<f64> {
        let mut y = DVector::from_element(x.len(), 0.);
        let mut i = 0;
        for layer in self.layers.iter() {
            for p in layer
                .cells
                .iter()
                .map(|c| c.borrow())
                .filter(|c| matches!(c.variant, Internal))
            {
                let mut flux = 0.;
                for f in p.faces.iter() {
                    let n = f.neighbor.as_ref().unwrap().cell.borrow();
                    flux -= p.lambda * (p.T - n.T) / f.neighbor.as_ref().unwrap().distance * f.area;
                }
                y[i] = flux / p.volume / p.rho / p.cp;
                i += 1;
            }
        }
        y
    }

    pub fn update(&mut self, y: &DVector<f64>) {
        // takes the solution and updates every cell with it
        let mut i = 0;
        for layer in self.layers.iter() {
            for mut cell in layer.cells.iter().map(|c| c.borrow_mut()) {
                if matches!(cell.variant, Internal) {
                    cell.T = y[i];
                    i += 1;
                }
            }
        }
        self.enforce_boundary_conditions();
    }
}

impl fmt::Debug for Domain {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}", self.layers[0])
    }
}

impl Default for Domain {
    fn default() -> Self {
        Domain::new()
    }
}

#[cfg(test)]
mod tests {
    use crate::cell::BoundaryCondition;
    use crate::solver;
    #[test]
    fn test_generation() {
        let boundary_conditions = [
            BoundaryCondition::Dirichlet(300.),
            BoundaryCondition::Neumann(0.),
            BoundaryCondition::Neumann(0.),
            BoundaryCondition::Dirichlet(0.),
            BoundaryCondition::Neumann(0.),
            BoundaryCondition::Neumann(0.),
        ];
        let (nx, ny, nz) = (3, 1, 1);
        let mut block: solver::Domain;
        block = solver::Domain::new();
        block.add_layer((0., 0., 0.), (1., 1., 1.), nx, ny, nz, boundary_conditions);
    }
}
