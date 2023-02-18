use crate::mesh;
use crate::cell;
use nalgebra::DVector;
use std::fmt;

#[derive(Default)]
pub struct Domain {
    layers: Vec<mesh::Mesh>,
    n_layers: usize,
}

impl Domain {
    pub fn new() -> Self {
        Domain{layers: Vec::new(), n_layers: 0}
    }
    pub fn add_layer(&mut self, r1:(f64,f64,f64), r2:(f64,f64,f64), nx:usize, ny:usize, nz:usize, boundary_conditions:Vec<cell::BoundaryDirection>) {
        self.layers.push(mesh::Mesh::generate(r1, r2, nx, ny, nz, self.n_layers,boundary_conditions));
        self.n_layers += 1;
    }

    fn link_layers(&mut self) {
        for layer_index in 0..self.layers.len() {
            for cell_index in 0..layers[layer_index].grid.len() {
                let neighbor = self.find_neighbor((layer_index,cell_index));
            }
        }

        for layer_index in 0..self.layers.len() {
            for cell_index in 0..layers[layer_index].grid.len() {
                self.layers[layer_index].grid[cell_index].de = self.get_cell(self.layers[layer_index].grid[cell_index].east).r - self.layers[layer_index].grid[cell_index].r;
    }

    fn get_cell(&self, (layer_id,cell_id):(usize,usize)) -> &cell::Cell {
        &self.layers[layer_id].grid[cell_id]
    }

    #[allow(dead_code)]
    fn find_neighbor(&self, (id,index):(usize,usize)) -> (usize,usize) {
        (0,0)// TBD
    }

    pub fn enforce_boundary_conditions(&mut self) {
        for layer_id in 0..self.layers.len() {
            for cell_id in 0..self.layers[layer_id].grid.len() {
                if let cell::CellType::Boundary(cell::BoundaryCondition::Dirichlet(val, _)) = self.layers[layer_id].grid[cell_id].variant {
                    self.layers[layer_id].grid[cell_id].T = val;
                } else if let cell::CellType::Boundary(cell::BoundaryCondition::Neumann(val, neighbor)) = self.layers[layer_id].grid[cell_id].variant {
                    self.layers[layer_id].grid[cell_id].T = self.get_cell(neighbor.unwrap()).T; // TODO: only true for val = 0
                }

            }
        }
    }

    pub fn setup_equations(&self, _t:&f64, x:&DVector<f64>) -> DVector<f64> {
        let mut y = DVector::from_element(x.len(), 0.);
        let mut i = 0;
        for layer in &self.layers {
            for p in layer.get_cells() {
                let radial;
                let axial;
                let n = self.get_cell(p.north);
                let s = self.get_cell(p.south);
                let w = self.get_cell(p.west);
                let e = self.get_cell(p.east);

                axial = p.rho * p.cp / p.dx * (
                    (e.T - p.T) / p.de * (e.lambda + p.lambda) / 2. - 
                    (p.T - w.T) / p.dw * (p.lambda + w.lambda) / 2. );
                radial = p.rho * p.cp / p.dy * (
                    (n.T - p.T) / p.dn * (n.lambda + p.lambda) / 2. - 
                    (p.T - s.T) / p.ds * (p.lambda + s.lambda) / 2. );
                y[i] = radial + axial;
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
                cell.T = y[i];
                i += 1;
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