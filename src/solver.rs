use crate::cell::CellType;
use crate::cell::Neighbor;
use crate::mesh;
use crate::cell;
use nalgebra::DVector;
use std::fmt;

#[derive(Default)]
pub struct Domain {
    layers: Vec<mesh::Mesh>,
}

impl Domain {
    pub fn new() -> Self {
        Domain{layers: Vec::new()}
    }

    pub fn add_layer(&mut self, r1:(f64,f64,f64), r2:(f64,f64,f64), nx:usize, ny:usize, nz:usize, bc:[cell::BoundaryCondition;6]) {
        self.layers.push(mesh::Mesh::generate(r1, r2, nx, ny, nz, self.layers.len(),bc));
        self.link_cells()
    }

    fn get_cell(&self, (layer_id,cell_id):(usize,usize)) -> &cell::Cell {
        &self.layers[layer_id].grid[cell_id]
    }

    fn link_cells(&mut self) {
        // link cells inside layers
        for layer_index in 0..self.layers.len() {
            for cell_index in 0..self.layers[layer_index].grid.len() {
                self.assign_neighbors((layer_index,cell_index));
            }
        }
    }

    fn assign_neighbors(&mut self, (layer_index,cell_index):(usize,usize)) {
        // takes the index of a cell (layer,cell) and assigns the correct value to its cell.face.neighbor
        // the neighbor must first be found. the criteris is that its cell.face.r is the same

        // each face has a neighbor
        for face_index in 0..self.layers[layer_index].grid[cell_index].faces.len() {
            if matches!(self.layers[layer_index].grid[cell_index].faces[face_index].neighbor,Some(_)) {continue;}
            // could be in any layer
            for neighbor_layer_index in 0..self.layers.len() {
                // could be any cell
                for neighbor_cell_index in 0..self.layers[neighbor_layer_index].grid.len() {
                    if (neighbor_layer_index,neighbor_cell_index) == (layer_index,cell_index) {continue}
                    // could be any of its faces
                    for neighbor_face_index in 0..self.layers[neighbor_layer_index].grid[neighbor_cell_index].faces.len() {
                        if matches!(self.layers[neighbor_layer_index].grid[neighbor_cell_index].faces[neighbor_face_index].neighbor, Some(_)) {continue;}
                        // check if the face centers coincide
                        if (self.layers[neighbor_layer_index].grid[neighbor_cell_index].faces[neighbor_face_index].center - self.layers[layer_index].grid[cell_index].faces[face_index].center).norm().abs() < f64::EPSILON {
                            let dr = self.layers[neighbor_layer_index].grid[neighbor_cell_index].r - self.layers[layer_index].grid[cell_index].r;
                            self.layers[layer_index].grid[cell_index].faces[face_index].neighbor = Some(Neighbor{
                                index: (neighbor_layer_index,neighbor_cell_index),
                                dr: dr,
                                distance: dr.norm(),
                            });
                            self.layers[neighbor_layer_index].grid[neighbor_cell_index].faces[neighbor_face_index].neighbor = Some(Neighbor{
                                index: (layer_index,cell_index),
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
                if let cell::CellType::Boundary(cell::BoundaryCondition::Dirichlet(val)) = self.layers[layer_id].grid[cell_id].variant {
                    self.layers[layer_id].grid[cell_id].T = val;
                } else if let cell::CellType::Boundary(cell::BoundaryCondition::Neumann(val)) = self.layers[layer_id].grid[cell_id].variant {
                    if self.layers[layer_id].grid[cell_id].faces.len() != 1 {
                        panic!("boundary cells should only have a single neighbor");
                    } else {
                        // TODO: only true for val = 0
                        self.layers[layer_id].grid[cell_id].T = self.get_cell(self.layers[layer_id].grid[cell_id].faces[0].neighbor.unwrap().index).T;
                    }
                    
                }
            }
        }
    }

    pub fn setup_equations(&self, _t:&f64, x:&DVector<f64>) -> DVector<f64> {
        let mut y = DVector::from_element(x.len(), 0.);
        let mut i = 0;
        for layer in &self.layers {
            for p in layer.get_cells()
                                 .iter()
                                 .filter(|c| matches!(c.variant,CellType::Internal)) {
                let mut flux = 0.;
                for f in &p.faces {
                    let n = self.get_cell(f.neighbor.as_ref().unwrap().index);
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
        for layer in &mut self.layers {
            for cell in layer.get_mut_cells() {
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