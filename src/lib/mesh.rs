use crate::computational_cell::{ComputationalCell, Face, HasProperties};
use nalgebra::{DVector, Point3, Vector3};

/// The mesh owns all geometrical features:
/// - vertices: points in space at the corners of computational cells
/// - faces: the faces of the computational cells
/// - computational cells
#[derive(Clone)]
pub struct Mesh<T: HasProperties> {
    pub vertices: Vec<Point3<f64>>,
    pub faces: Vec<Face>,
    pub cells: Vec<ComputationalCell<T>>,
}

impl<T: HasProperties> Mesh<T> {
    pub fn new() -> Self {
        Mesh {
            vertices: Vec::new(),
            faces: Vec::new(),
            cells: Vec::new(),
        }
    }

    pub fn n_cells(&self) -> usize {
        // returns the total number of cells in the domain
        self.cells.len()
    }

    pub fn add_hex_cell(&mut self, r1: Point3<f64>, r2: Point3<f64>, properties: T) -> usize {
        // create central node if it does not exist, or get the index of the existing one
        let node_index = match self
            .nodes
            .iter()
            .position(|n| (r1 + (r2 - r1) / 2. - n).norm() < f64::EPSILON)
        {
            Some(index) => index,
            None => {
                let node_index = self.nodes.len();
                self.nodes.push(r1 + (r2 - r1) / 2.);
                node_index
            }
        };
        // create vertices if they do not exist, else retrieve their index
        let mut current_vertices = vec![];
        let mut current_vertices_index = vec![];
        for i in 0..2 {
            for j in 0..2 {
                for k in 0..2 {
                    let r = Point3::new(
                        r1[0] + i as f64 * (r2[0] - r1[0]),
                        r1[1] + j as f64 * (r2[1] - r1[1]),
                        r1[2] + k as f64 * (r2[2] - r1[2]),
                    );
                    current_vertices.push(r);
                    // it this vertex exists, find its index
                    let vertex_index = match self
                        .vertices
                        .iter()
                        .position(|v| (r - v).norm() < f64::EPSILON)
                    {
                        Some(index) => index,
                        None => {
                            let vertex_index = self.vertices.len();
                            self.vertices.push(r);
                            vertex_index
                        }
                    };
                    current_vertices_index.push(vertex_index);
                }
            }
        }

        let cell_index = self.cells.len();
        self.cells.push(ComputationalCell::create_rectangle(
            node_index,
            current_vertices_index,
            current_vertices,
            properties,
        ));
        cell_index
    }

    pub fn get_face_on_plane(
        &mut self,
        point: Point3<f64>,
        normal: Vector3<f64>,
        cell_index: usize,
    ) -> Option<&mut Face> {
        // find the face that is on the plane defined by the point and normal
        let cell = &mut self.cells[cell_index];
        cell.faces.iter_mut().find(|f| {
            let mut scalar_vec = Vec::new();
            for id in f.vertices.iter() {
                scalar_vec.push(normal.dot(&(self.vertices[*id] - point)))
            }
            DVector::from_vec(scalar_vec).norm() < f64::EPSILON
        })
    }
}

impl<T: HasProperties> Default for Mesh<T> {
    fn default() -> Self {
        Mesh::new()
    }
}
