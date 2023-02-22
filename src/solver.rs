use crate::cell::{CellType,Neighbor};
use crate::mesh;
use crate::cell;
use nalgebra::DVector;
use std::rc::Rc;
use std::fmt;

pub struct Domain {
    layers: Vec<mesh::Mesh>,
}

impl Domain {
    pub fn new() -> Self {
        Domain{layers: Vec::new()}
    }

    pub fn add_layer(&mut self, r1:(f64,f64,f64), r2:(f64,f64,f64), nx:usize, ny:usize, nz:usize, bc:[cell::BoundaryCondition;6]) {
        self.layers.push(mesh::Mesh::generate(r1, r2, nx, ny, nz, bc));
        self.link_cells()
    }

    fn link_cells(&self) {
        // link cells inside layers
        for layer_index in 0..self.layers.len() {
            for cell_index in 0..self.layers[layer_index].grid.len() {
                self.assign_neighbors((layer_index,cell_index));
            }
        }
    }

    fn assign_neighbors(&self, (layer_index,cell_index):(usize,usize)) {
        // takes the index of a cell (layer,cell) and assigns the correct value to its cell.face.neighbor
        // the neighbor must first be found. the criteris is that its cell.face.r is the same

        // each face has a neighbor
        let n_faces = self.layers[layer_index].grid[cell_index].borrow().faces.len();
        for face_index in 0..n_faces {
            if matches!(self.layers[layer_index].grid[cell_index].borrow().faces[face_index].neighbor,Some(_)) {continue}
            // could be in any layer
            for neighbor_layer_index in 0..self.layers.len() {
                // could be any cell
                for neighbor_cell_index in 0..self.layers[neighbor_layer_index].grid.len() {
                    if (neighbor_layer_index,neighbor_cell_index) == (layer_index,cell_index) {continue}
                    // could be any of its faces
                    let n_faces_neighbor = self.layers[neighbor_layer_index].grid[neighbor_cell_index].borrow().faces.len();
                    for neighbor_face_index in 0..n_faces_neighbor {
                        if matches!(self.layers[neighbor_layer_index].grid[neighbor_cell_index].borrow().faces[neighbor_face_index].neighbor, Some(_)) {continue}
                        // check if the face centers coincide
                        if (self.layers[neighbor_layer_index].grid[neighbor_cell_index].borrow().faces[neighbor_face_index].center - self.layers[layer_index].grid[cell_index].borrow().faces[face_index].center).norm().abs() < f64::EPSILON {
                            let dr = self.layers[neighbor_layer_index].grid[neighbor_cell_index].borrow().r - self.layers[layer_index].grid[cell_index].borrow().r;
                            self.layers[layer_index].grid[cell_index].borrow_mut().faces[face_index].neighbor = Some(Neighbor{
                                cell: Rc::clone(&self.layers[neighbor_layer_index].grid[neighbor_cell_index]),
                                dr: dr,
                                distance: dr.norm(),
                            });
                            self.layers[neighbor_layer_index].grid[neighbor_cell_index].borrow_mut().faces[neighbor_face_index].neighbor = Some(Neighbor{
                                cell: Rc::clone(&self.layers[layer_index].grid[cell_index]),
                                dr: -dr,
                                distance: dr.norm(),
                            });
                        }
                    }
                }
            }
        }
    }

    pub fn enforce_boundary_conditions(&mut self) {
        for layer_id in 0..self.layers.len() {
            for cell_id in 0..self.layers[layer_id].grid.len() {
                let bc = self.layers[layer_id].grid[cell_id].borrow().variant.clone();
                if let cell::CellType::Boundary(cell::BoundaryCondition::Dirichlet(val)) = bc {
                    self.layers[layer_id].grid[cell_id].borrow_mut().T = val;
                } else if let cell::CellType::Boundary(cell::BoundaryCondition::Neumann(val)) = bc {
                    if self.layers[layer_id].grid[cell_id].borrow().faces.len() != 1 {
                        panic!("boundary cells should only have a single neighbor");
                    } else {
                        // TODO: only true for val = 0
                        if self.layers[layer_id].grid[cell_id].borrow().faces.len() != 1 {panic!("why does the boundary cell have more than one neighbor ?")}
                        let boundary_value = self.layers[layer_id].grid[cell_id].borrow().faces[0].neighbor.as_ref().unwrap().cell.borrow().T; // gotta check if there is more than one face
                        self.layers[layer_id].grid[cell_id].borrow_mut().T = boundary_value;
                    }
                    
                }
            }
        }
    }

    pub fn setup_equations(&self, _t:&f64, x:&DVector<f64>) -> DVector<f64> {
        let mut y = DVector::from_element(x.len(), 0.);
        let mut i = 0;
        for layer in self.layers.iter() {
            for p in layer.grid
                                            .iter()
                                            .map(|c| c.borrow())
                                            .filter(|c| matches!(c.variant, CellType::Internal)) {
                let mut flux = 0.;
                for f in p.faces.iter() {
                    let n = f.neighbor.as_ref().unwrap().cell.borrow();
                    flux -= p.lambda * (p.T - n.T) / f.neighbor.as_ref().unwrap().distance * f.area;
                }
                y[i] = flux/p.volume/p.rho/p.cp;
                i += 1;
            }
        }
        y
    }

    pub fn update(&mut self, y:&DVector<f64>) {
        // takes the solution and updates every cell with it
        let mut i = 0;
        for layer in self.layers.iter() {
            for mut cell in layer.grid.iter().map(|c| c.borrow_mut()) {
                if matches!(cell.variant,CellType::Internal) {
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