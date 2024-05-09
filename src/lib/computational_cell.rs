// for now all cells are rectangular parallelepipeds
use crate::geometry_utils::*;
use crate::mesh::Mesh;
use nalgebra::{Point3, Vector3};
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
}

pub trait HasProperties: Default + Clone + Copy + fmt::Debug {
    fn update(&mut self); // this function only updates derived variables, computed variables are updated by the equations
}

impl<T: HasProperties> ComputationalCell<T> {
    /// Creates a new computational cell from a list of faces.
    /// The cell is not stored in the mesh, only its properties are initialized.
    pub fn new(mesh: &Mesh<T>, faces: Vec<usize>, properties: T) -> Self {
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
        };
        cell
    }
}

impl<T: HasProperties> fmt::Debug for ComputationalCell<T> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.properties.fmt(f)
    }
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

    /// Returns a tuple with the index of the neighbor to the cell passed as argument and the distance to it.
    pub fn neighbor(&self, i: usize) -> Option<(usize, f64)> {
        // check if i is the owner or the neighbour
        match self.owner.0 {
            i => match self.neighbor {
                NeighborType::Internal(j, distance) => Some((j, distance + self.owner.1)),
                _ => None,
            },
            _ => match self.neighbor {
                NeighborType::Internal(_, distance) => Some((i, distance + self.owner.1)),
                _ => None,
            },
        }
    }
}
