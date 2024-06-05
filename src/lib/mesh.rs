use crate::computational_cell::{CellData, FaceMut, Faces, NeighborType};
use crate::geometry_utils::{
    barycentric_coordinates, line_plane_intersection, polyhedron_center, polyhedron_volume,
};
use itertools::Itertools;
use nalgebra::{Point3, Vector3};
use std::collections::HashSet;

/// The mesh owns all geometrical features:
/// - vertices: points in space at the corners of computational cells
/// - faces: the faces of the computational cells
/// - computational cells
#[derive(Clone)]
pub struct Mesh {
    pub vertices: Vec<Point3<f64>>,
    pub faces: Faces,
    pub cells: CellData,
    pub periodic_pairs: Vec<(usize, Option<usize>)>,
    // optional, needs to be initialized with "collect_vertex_cells" before use
    pub vertex_cells: Vec<Vec<usize>>,
}

impl Mesh {
    pub fn new() -> Self {
        Mesh {
            vertices: Vec::new(),
            faces: Faces::new(),
            cells: CellData::default(),
            periodic_pairs: Vec::new(),
            vertex_cells: Vec::new(),
        }
    }

    pub fn n_cells(&self) -> usize {
        // returns the total number of cells in the domain
        self.cells.r.len()
    }

    fn vertex_index(&self, r: Point3<f64>) -> Option<usize> {
        // returns the index of the vertex if it exists
        // we parse the vertices in reverse order because what we are looking for is likely to be at the end of the list
        let return_value = self
            .vertices
            .iter()
            .rev()
            .position(|v| (v - r).norm() < f64::EPSILON)?;
        Some(self.vertices.len() - 1 - return_value)
    }
    /// Creates a new face from a list of point indices and initialize its basic geometrical properties.
    /// The properties relating to the neighbour need to be initialized later.
    pub fn add_face_from_vertices(&mut self, vertex_indices: &[usize]) -> usize {
        let vertex_coorinates: &[Point3<f64>] = &vertex_indices
            .iter()
            .map(|i_vertex| self.vertices[*i_vertex])
            .collect::<Vec<Point3<f64>>>();
        self.faces
            .add_from_vertices(vertex_coorinates, vertex_indices)
    }

    pub fn add_hex_cell(&mut self, r1: Point3<f64>, r2: Point3<f64>, zone: usize) -> usize {
        // assumes that the cell does not exist yet
        let cell_index = self.cells.volume.len();
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
                    match self.vertex_index(r) {
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
        // Create faces if they do not exist, else retrieve their index.
        // The faces need to follow the right-hand rule, meaning a counter-clockwise order of vertices
        // which leads to the normal pointing outwards of the cell.
        // There must be a smart way to achieve this, but for now I will hardcode the stupid way.
        // 6 faces, 4 vertices each
        let mut current_faces = [[usize::MAX; 4]; 6];
        // left face : i = 0
        current_faces[0] = [
            current_vertices_index[0][0][0],
            current_vertices_index[0][0][1],
            current_vertices_index[0][1][1],
            current_vertices_index[0][1][0],
        ];
        // right face : i = 1
        current_faces[1] = [
            current_vertices_index[1][0][0],
            current_vertices_index[1][1][0],
            current_vertices_index[1][1][1],
            current_vertices_index[1][0][1],
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
            current_vertices_index[0][1][0],
            current_vertices_index[1][1][0],
            current_vertices_index[1][0][0],
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
        let mut cell_faces = Vec::with_capacity(6);
        for vertices in current_faces.iter() {
            let mut found = false;
            for i_f in 0..self.faces.len() {
                let face = self.faces.get_face_mut(i_f);
                // WARNING : this assumes that no crazy face with a different number of vertices exists
                if face.vertices.iter().all(|p| vertices.contains(p)) {
                    cell_faces.push(i_f);
                    let face_distance = (*face.center - cell_center).norm();
                    *face.neighbor = NeighborType::Internal(cell_index, face_distance);
                    found = true;
                    break;
                }
            }
            if !found {
                let i_f = self.add_face_from_vertices(vertices);
                let face_distance = (self.faces.center[i_f] - cell_center).norm();
                // face normal points from the owner to the neighbour
                self.faces.owner[i_f] = (cell_index, face_distance);
                cell_faces.push(i_f);
            }
        }
        // check if the faces are correctly oriented, meaning that the face normals point outwards of the owner and inwards of the neighbour
        for &i_face in cell_faces.iter() {
            let face = self.faces.get_face(i_face);
            let to_face = face.center - cell_center;
            if cell_index == face.owner.0 {
                if to_face.dot(&face.normal) < 0. {
                    panic!("Face normal points into the owner cell")
                }
            } else if to_face.dot(&face.normal) > 0. {
                panic!("Face normal points into the neighbor cell")
            }
        }
        // create the cell
        // this is slow considering we actually already have all the necessary information,
        // but I want to test/showcase that function for now
        self.new_cell(cell_faces, zone);
        cell_index
    }

    /// Creates a new computational cell from a list of faces.
    /// The cell is not stored in the mesh, only its properties are initialized.
    pub fn new_cell(&mut self, faces: Vec<usize>, zone: usize) {
        let mut vertex_indices: Vec<usize> = faces
            .iter()
            .flat_map(|&i_face| self.faces.vertices[i_face].clone())
            .collect();
        vertex_indices.sort();
        vertex_indices.dedup();
        let vertex_coordinates: Vec<Point3<f64>> = vertex_indices
            .iter()
            .map(|i_vertex| self.vertices[*i_vertex])
            .collect();
        let r = polyhedron_center(&vertex_coordinates);
        let volume = polyhedron_volume(&self.faces_to_points(&faces));
        self.cells.push(faces, r, volume, zone);
    }

    pub fn get_face_on_plane(
        &mut self,
        point: &Point3<f64>,
        normal: &Vector3<f64>,
        cell_index: usize,
    ) -> Option<FaceMut> {
        match self.get_face_id_on_plane(point, normal, cell_index) {
            Some(i_f) => Some(self.faces.get_face_mut(i_f)),
            None => None,
        }
    }

    pub fn get_face_id_on_plane(
        &mut self,
        point: &Point3<f64>,
        normal: &Vector3<f64>,
        cell_index: usize,
    ) -> Option<usize> {
        // find the face that is on the plane defined by the point and normal
        for &i_f in self.cells.faces[cell_index].iter() {
            if (self.faces.center[i_f] - point).dot(normal).abs() < f64::EPSILON {
                return Some(i_f);
            }
        }
        None // Return None if no face is found
    }

    pub fn faces_to_points(&self, faces: &[usize]) -> Vec<Vec<Point3<f64>>> {
        // returns the points of the faces
        faces
            .iter()
            .map(|&i_face| {
                self.faces.vertices[i_face]
                    .iter()
                    .map(|i_vertex| self.vertices[*i_vertex])
                    .collect::<Vec<Point3<f64>>>()
            })
            .collect::<Vec<Vec<Point3<f64>>>>()
    }

    /// initialized the "vertex_cells" field of the mesh with the list of vertices that are part of each cell
    pub fn collect_cell_vertices(&mut self) {
        self.cells.vertices.reserve(self.n_cells());
        for list_of_faces in self.cells.faces.iter() {
            self.cells.vertices.push(
                list_of_faces
                    .iter()
                    .flat_map(|&i_face| self.faces.vertices[i_face].clone())
                    .collect::<HashSet<usize>>(),
            )
        }
    }

    /// initialized the "vertex_cells" field of the mesh with the list of cells that contain each vertex
    pub fn collect_vertex_cells(&mut self) {
        if self.cells.vertices.is_empty() {
            self.collect_cell_vertices();
        }
        self.vertex_cells = vec![Vec::new(); self.vertices.len()];
        for (i_cell, vertices) in self.cells.vertices.iter().enumerate() {
            for &i_vertex in vertices.iter() {
                self.vertex_cells[i_vertex].push(i_cell);
            }
        }
    }

    /// returns a list of cells that have at least one vertex in common with the cell i_cell
    /// the "vertex_cells" field of the mesh must be initialized before calling this function
    pub fn extended_neighborhood(&self, i_cell: usize) -> Vec<usize> {
        let mut neighbor_set = HashSet::new();
        let cell_vertices = &self.cells.vertices[i_cell];
        for &i_vertex in cell_vertices.iter() {
            for &i_neighbor in self.vertex_cells[i_vertex].iter() {
                neighbor_set.insert(i_neighbor);
            }
        }
        neighbor_set.remove(&i_cell);
        neighbor_set.into_iter().collect()
    }

    /// interpolates a point "P" on the line between "Mij" and "Bi" using barycentric interpolation
    /// of points in the sorted list of neighbors. The first points in the list are prioritized,
    /// and if more than 3 points are provided, the subsequent points are only used if the line does not
    /// intersect the triangle defined by the first 3 points.
    /// The result is a list of tuples with the index of the points used for interpolation and their corresponding
    /// barycentric coordinates.
    #[allow(non_snake_case)]
    pub fn best_barycentric_interpolation(
        &self,
        Mij: &Point3<f64>,
        Bi: &Point3<f64>,
        sorted_list_of_neighbors: &[usize],
    ) -> Option<[(usize, f64); 3]> {
        for list in sorted_list_of_neighbors
            .iter()
            .enumerate()
            // we consider all combinations of 3 points and sort them by the sum of their indices
            .combinations(3)
            .sorted_by(|a, b| {
                Ord::cmp(
                    &a.iter().map(|a| a.0).sum::<usize>(),
                    &b.iter().map(|b| b.0).sum::<usize>(),
                )
            })
        {
            let B1 = &self.cells.r[*list[0].1];
            let B2 = &self.cells.r[*list[1].1];
            let B3 = &self.cells.r[*list[2].1];
            let plane_normal = (B2 - B1).cross(&(B3 - B1));
            // compute the intersection of BiMij with the plane defined by the three points
            if let Some(H) = line_plane_intersection(Mij, Bi, B1, &plane_normal) {
                // we compute the barycentric coordinates of H in the triangle B1B2B3
                let (alpha, beta, gamma) = barycentric_coordinates(B1, B2, B3, &H);
                // check if the point is in the triangle
                if alpha >= 0.0 && beta >= 0.0 && gamma >= 0.0 {
                    // we have found the point H, now we can interpolate the properties of the point
                    return Some([(*list[0].1, alpha), (*list[1].1, beta), (*list[2].1, gamma)]);
                }
            }
        }
        None
    }
}

impl Default for Mesh {
    fn default() -> Self {
        Mesh::new()
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn single_cell_mesh() {
        let mut mesh: crate::mesh::Mesh = crate::mesh::Mesh::new();
        mesh.add_hex_cell(
            nalgebra::Point3::new(0., 0., 0.),
            nalgebra::Point3::new(1., 1., 1.),
            0,
        );
        assert_eq!(mesh.n_cells(), 1);
        assert_eq!(mesh.vertices.len(), 8);
        assert_eq!(mesh.faces.len(), 6);
        assert_eq!(mesh.cells.faces[0].len(), 6);
        mesh.cells.faces[0].iter().for_each(|i_face| {
            let face = mesh.faces.get_face(*i_face);
            assert_eq!(face.owner.0, 0);
            assert_eq!(face.owner.1, 0.5);
            assert_eq!(*face.area, 1.0);
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
                let v = &mesh.vertices[face.vertices[i]] - face.center;
                assert!(v.dot(&face.normal).abs() < f64::EPSILON);
            }
            // check that the normal is pointing outwards of the cell
            let r = face.center - mesh.cells.r[0];
            assert!(face.normal.dot(&r) > 0.);
        });

        assert!(mesh.cells.volume[0] - 1.0 < f64::EPSILON);
        assert_eq!(mesh.cells.r[0], nalgebra::Point3::new(0.5, 0.5, 0.5));
    }
}
