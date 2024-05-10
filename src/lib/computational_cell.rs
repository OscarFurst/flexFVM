// for now all cells are rectangular parallelepipeds
use crate::geometry_utils::*;
use crate::mesh::Mesh;
use nalgebra::{DVector, Point3, Vector3};
use std::fmt;

///////////////////////////
/// CELL and PROPERTIES ///
///////////////////////////

/// ComputationalCell is the main data structure of the mesh.
/// It does not directly store its node and vertices, but instead
/// stores their indices in the mesh.
/// The face itself is stored in the cell because it is unique to
/// the cell: the neighboring cell's face has a normal vector pointing
/// in the opposite direction, for example.
#[derive(Clone)]
pub struct ComputationalCell<T: HasProperties> {
    // faces of the volume (defines it)
    pub faces: Vec<usize>,
    // geometric properties
    pub r: Point3<f64>, // m
    pub volume: f64,    // m3
    pub properties: T,
    pub zone: usize,
}

pub trait HasProperties: Default + Clone + Copy + fmt::Debug {
    fn update(&mut self); // this function only updates derived variables, computed variables are updated by the equations
    fn export(&self) -> Vec<(String, f64)>;
}

impl<T: HasProperties> ComputationalCell<T> {
    /// Creates a new computational cell from a list of faces.
    /// The cell is not stored in the mesh, only its properties are initialized.
    pub fn new(mesh: &Mesh<T>, faces: Vec<usize>, properties: T, zone: usize) -> Self {
        let mut vertex_indices: Vec<usize> = faces
            .iter()
            .flat_map(|i_face| mesh.faces[*i_face].vertices.clone())
            .collect();
        vertex_indices.sort();
        vertex_indices.dedup();
        let vertex_coordinates: Vec<Point3<f64>> = vertex_indices
            .iter()
            .map(|i_vertex| mesh.vertices[*i_vertex])
            .collect();
        let r = polyhedron_center(&vertex_coordinates);
        let volume = polyhedron_volume(&mesh.faces_to_points(&faces));
        let cell = ComputationalCell {
            faces,
            r,
            volume,
            properties,
            zone,
        };
        cell
    }

    /// Returns a list of vertex indices which is correctly sorted for VTK.
    /// Only works for hexahedrons.
    pub fn hexahedron_vertices(&self, mesh: &Mesh<T>) -> Vec<usize> {
        // check that it is a hexahedron
        if self.faces.len() != 6
            || self
                .faces
                .iter()
                .any(|i_face| mesh.faces[*i_face].vertices.len() != 4)
        {
            panic!("Only hexahedra are supported for VTK export");
        }
        // create a vec of unique vertices
        let mut vertex_list: Vec<usize> = Vec::new();
        // start by adding the vertices of the first face
        vertex_list.extend(mesh.faces[self.faces[0]].vertices.iter());
        // add the vertices of the opposite face
        // opposite face is the one that shares no vertex with the first face
        let opposite_face = self
            .faces
            .iter()
            .find(|i_face| {
                mesh.faces[**i_face]
                    .vertices
                    .iter()
                    .all(|i_vertex| !vertex_list.contains(i_vertex))
            })
            .unwrap();
        let opposite_face = &mesh.faces[*opposite_face];
        let mut opposite_vertices = opposite_face.vertices.clone();
        // if the normal of the opposite face is pointing in the other direction, we need to reverse the order of the vertices
        if opposite_face.normal.dot(&mesh.faces[self.faces[0]].normal) < 0. {
            opposite_vertices.reverse();
        }
        // we shift the vertices of the "opposite_vertices" as long as the first vertex is not directly opposite to the first vertex of the "vertex_list"
        // being opposite means that the vector between the two vertices is parallel to the normal of the face
        let mut i = 0;
        while !self.are_connected(mesh, vertex_list[0], opposite_vertices[0]) {
            opposite_vertices.rotate_left(1);
            i += 1;
            if i > 4 {
                panic!("Could not find the opposite face");
            }
        }
        vertex_list.extend(opposite_vertices.iter());
        vertex_list
    }

    /// Returns true if the two vertices share two faces of the computational cell.
    fn are_connected(&self, mesh: &Mesh<T>, i_vertex1: usize, i_vertex2: usize) -> bool {
        let faces1: Vec<usize> = self
            .faces
            .iter()
            .map(|i_face| &mesh.faces[*i_face])
            .enumerate()
            .filter(|(i, face)| face.vertices.contains(&i_vertex1))
            .map(|(i, _)| i)
            .collect();
        let faces2: Vec<usize> = self
            .faces
            .iter()
            .map(|i_face| &mesh.faces[*i_face])
            .enumerate()
            .filter(|(i, face)| face.vertices.contains(&i_vertex2))
            .map(|(i, _)| i)
            .collect();
        let mut common_faces = 0;
        for i_face1 in faces1.iter() {
            for i_face2 in faces2.iter() {
                if *i_face1 == *i_face2 {
                    common_faces += 1;
                }
            }
        }
        common_faces == 2
    }
}

impl<T: HasProperties> fmt::Debug for ComputationalCell<T> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.properties.fmt(f)
    }
}

#[allow(dead_code)]
fn dedup(v: &mut Vec<usize>) {
    let mut set = std::collections::HashSet::new();

    v.retain(|x| set.insert(*x));
}

////////////
/// FACE ///
////////////
#[derive(Clone)]
pub enum NeighborType {
    Internal(usize, f64),
    Boundary(usize),
    Nil,
}

#[derive(Clone)]
pub struct Face {
    pub vertices: Vec<usize>,
    pub center: Point3<f64>,
    pub normal: Vector3<f64>,
    pub area: f64,
    pub owner: (usize, f64),
    pub neighbor: NeighborType,
}

impl Face {
    /// Creates a new face and initialize its basic geometrical properties.
    /// The properties relating to the neighbour need to be initialized later.
    pub fn new<T: HasProperties>(mesh: &Mesh<T>, vertices_indices: &[usize]) -> Self {
        let vertices_coordinates: &[Point3<f64>] = &vertices_indices
            .iter()
            .map(|i_vertex| mesh.vertices[*i_vertex])
            .collect::<Vec<Point3<f64>>>();
        // TODO: only implemented for rectangle
        // computing the center of the face
        let center: Point3<f64> = polyhedron_center(vertices_coordinates);
        let normal: Vector3<f64> = polygnon_normal(vertices_coordinates);
        let area: f64 = polyngon_area(vertices_coordinates);
        Face {
            vertices: vertices_indices.to_vec(),
            center,
            normal,
            area,
            owner: (usize::MAX, 0.),
            neighbor: NeighborType::Nil,
        }
    }
}
