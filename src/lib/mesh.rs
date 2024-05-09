use crate::computational_cell::{ComputationalCell, Face, HasProperties, NeighborType};
use crate::geometry_utils::polyhedron_center;
use itertools::Itertools;
use nalgebra::{Point3, Vector3};

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
        // assumes that the cell does not exist yet
        let cell_index = self.cells.len();
        // create vertices if they do not exist, else retrieve their index
        let mut current_vertices = [[[Point3::origin(); 2]; 2]; 2];
        let mut current_vertices_index = [[[usize::MAX; 2]; 2]; 2];
        for i in 0..2 {
            for j in 0..2 {
                for k in 0..2 {
                    let r = Point3::new(
                        r1[0] + i as f64 * (r2[0] - r1[0]),
                        r1[1] + j as f64 * (r2[1] - r1[1]),
                        r1[2] + k as f64 * (r2[2] - r1[2]),
                    );
                    current_vertices[i][j][k] = r;
                    // it this vertex exists, find its index
                    match self
                        .vertices
                        .iter()
                        .position(|v| (v - r).norm() < f64::EPSILON)
                    {
                        Some(index) => current_vertices_index[i][j][k] = index,
                        None => {
                            let index = self.vertices.len();
                            self.vertices.push(r);
                            current_vertices_index[i][j][k] = index;
                        }
                    }
                }
            }
        }
        let cell_center = polyhedron_center(
            current_vertices
                .iter()
                .flatten()
                .flatten()
                .cloned()
                .collect::<Vec<_>>()
                .as_slice(),
        );
        // create faces if they do not exist, else retrieve their index
        // 6 faces, 4 vertices each
        let mut current_faces = [[usize::MAX; 4]; 6];
        // left face : i = 0
        current_faces[0] = [
            current_vertices_index[0][0][0],
            current_vertices_index[0][1][0],
            current_vertices_index[0][1][1],
            current_vertices_index[0][0][1],
        ];
        // right face : i = 1
        current_faces[1] = [
            current_vertices_index[1][0][0],
            current_vertices_index[1][0][1],
            current_vertices_index[1][1][1],
            current_vertices_index[1][1][0],
        ];
        // front face : j = 0
        current_faces[2] = [
            current_vertices_index[0][0][0],
            current_vertices_index[1][0][0],
            current_vertices_index[1][0][1],
            current_vertices_index[0][0][1],
        ];
        // back face : j = 1
        current_faces[3] = [
            current_vertices_index[0][1][0],
            current_vertices_index[0][1][1],
            current_vertices_index[1][1][1],
            current_vertices_index[1][1][0],
        ];
        // bottom face : k = 0
        current_faces[4] = [
            current_vertices_index[0][0][0],
            current_vertices_index[1][0][0],
            current_vertices_index[1][1][0],
            current_vertices_index[0][1][0],
        ];
        // top face : k = 1
        current_faces[5] = [
            current_vertices_index[0][0][1],
            current_vertices_index[1][0][1],
            current_vertices_index[1][1][1],
            current_vertices_index[0][1][1],
        ];
        // for each of these, check if there is a face with the same vertices
        // if not, create a new face
        let mut cell_faces = Vec::new();
        for vertices in current_faces.iter() {
            let mut found = false;
            for (index, face) in self.faces.iter_mut().enumerate() {
                if face.vertices.iter().sorted().eq(vertices.iter().sorted()) {
                    cell_faces.push(index);
                    let face_distance = (face.center - cell_center).norm();
                    face.neighbor = NeighborType::Internal(cell_index, face_distance);
                    found = true;
                    break;
                }
            }
            if !found {
                let mut face = Face::new(&self, vertices);
                let face_distance = (face.center - cell_center).norm();
                // face normal points from the owner to the neighbour
                face.owner = (cell_index, face_distance);
                if (face.center - cell_center).dot(&face.normal) < 0. {
                    face.normal *= -1.;
                }
                cell_faces.push(self.faces.len());
                self.faces.push(face);
            }
        }

        self.cells
            .push(ComputationalCell::new(&self, cell_faces, properties));
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
        for face_index in cell.faces.iter() {
            let face = &self.faces[*face_index];
            if (face.center - point).dot(&normal).abs() < f64::EPSILON {
                return Some(&mut self.faces[*face_index]);
            }
        }
        None // Return None if no face is found
    }

    pub fn faces_to_points(&self, faces: &[usize]) -> Vec<Vec<Point3<f64>>> {
        // returns the points of the faces
        faces
            .iter()
            .map(|i_face| {
                self.faces[*i_face]
                    .vertices
                    .iter()
                    .map(|i_vertex| self.vertices[*i_vertex])
                    .collect::<Vec<Point3<f64>>>()
            })
            .collect::<Vec<Vec<Point3<f64>>>>()
    }
}

impl<T: HasProperties> Default for Mesh<T> {
    fn default() -> Self {
        Mesh::new()
    }
}

#[cfg(test)]
mod tests {
    use crate::computational_cell::HasProperties;

    #[derive(Default, Copy, Clone, Debug)]
    pub struct Properties {}
    impl HasProperties for Properties {
        fn update(&mut self) {}
    }

    #[test]
    fn single_cell_mesh() {
        let mut mesh: crate::mesh::Mesh<Properties> = crate::mesh::Mesh::new();
        mesh.add_hex_cell(
            nalgebra::Point3::new(0., 0., 0.),
            nalgebra::Point3::new(1., 1., 1.),
            Properties::default(),
        );
        assert_eq!(mesh.n_cells(), 1);
        assert_eq!(mesh.vertices.len(), 8);
        assert_eq!(mesh.faces.len(), 6);
        assert_eq!(mesh.cells[0].faces.len(), 6);
        mesh.cells[0].faces.iter().for_each(|i_face| {
            let face = &mesh.faces[*i_face];
            assert_eq!(face.owner.0, 0);
            assert_eq!(face.owner.1, 0.5);
            assert_eq!(face.area, 1.0);
            // check that the face is a convex polygon with sorted vertices
            assert!(crate::geometry_utils::check_polygon_convexity(
                face.vertices
                    .iter()
                    .map(|i| mesh.vertices[*i])
                    .collect::<Vec<_>>()
                    .as_slice()
            ));
            // check that any vector from the center to the vertices is normal to the face normal
            for i in 0..face.vertices.len() {
                let v = &mesh.vertices[face.vertices[i]] - &face.center;
                assert!(v.dot(&face.normal).abs() < f64::EPSILON);
            }
            // check that the normal is pointing outwards of the cell
            let r = face.center - mesh.cells[0].r;
            assert!(face.normal.dot(&r) > 0.);
        });

        assert!(mesh.cells[0].volume - 1.0 < f64::EPSILON);
        assert_eq!(mesh.cells[0].r, nalgebra::Point3::new(0.5, 0.5, 0.5));
    }
}
