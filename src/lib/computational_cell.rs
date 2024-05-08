// for now all cells are rectangular parallelepipeds
use nalgebra::{Point3, Vector3};
use std::fmt;

use crate::equation;

///////////////////////////
/// CELL and PROPERTIES ///
///////////////////////////

/// ComputationalCell is the main data structure of the mesh.
/// It does not directly store its node and vertices, but instead
/// stores their indices in the mesh.
/// The face itself is stored in the cell because it is unique to
/// the cell: the neighboring cell's face has a normal vector pointing
/// in the opposite direction, for example.
#[allow(non_snake_case)]
#[derive(Clone)]
pub struct ComputationalCell<T: HasProperties> {
    // geometric properties
    pub r: usize,    // m
    pub volume: f64, // m3
    pub vertices: Vec<usize>,
    pub properties: T,
    // faces of the volume
    pub faces: Vec<Face>,
}

pub trait HasProperties: Default + Clone + Copy + fmt::Debug {
    fn update(&mut self); // this function only updates derived variables, computed variables are updated by the equations
}

impl<T: HasProperties> ComputationalCell<T> {
    pub fn create_rectangle(
        r: usize,
        vertex_indices: Vec<usize>,
        vertex_points: Vec<Point3<f64>>,
        properties: T,
    ) -> ComputationalCell<T> {
        // compute volume from vertex points
        let volume = {
            let v1 = vertex_points[1] - vertex_points[0];
            let v2 = vertex_points[2] - vertex_points[0];
            let v3 = vertex_points[4] - vertex_points[0];
            v1.dot(&v2.cross(&v3)).abs()
        };
        // create 6 faces of the volume from vertex points and indices
        let faces = vec![
            Face::new(
                vec![vertex_indices[0], vertex_indices[1], vertex_indices[2]],
                &vertex_points,
            ),
            Face::new(
                vec![vertex_indices[0], vertex_indices[1], vertex_indices[4]],
                &vertex_points,
            ),
            Face::new(
                vec![vertex_indices[0], vertex_indices[2], vertex_indices[4]],
                &vertex_points,
            ),
            Face::new(
                vec![vertex_indices[1], vertex_indices[2], vertex_indices[4]],
                &vertex_points,
            ),
            Face::new(
                vec![vertex_indices[0], vertex_indices[3], vertex_indices[5]],
                &vertex_points,
            ),
            Face::new(
                vec![vertex_indices[1], vertex_indices[3], vertex_indices[5]],
                &vertex_points,
            ),
        ];
        let mut cell = ComputationalCell {
            r,
            volume,
            vertices: vertex_indices,
            properties,
            faces,
        };
        cell
    }
}

impl<T: HasProperties> fmt::Debug for ComputationalCell<T> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.properties.fmt(f)
    }
}

impl<T: HasProperties> Default for ComputationalCell<T> {
    fn default() -> Self {
        ComputationalCell {
            r: usize::MAX,
            vertices: Vec::new(),
            volume: 0.,
            properties: T::default(),
            faces: Vec::new(),
        }
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
    fn new(vertices_indices: Vec<usize>, vertices_coordinates: &Vec<Point3<f64>>) -> Self {
        // TODO: only implemented for rectangle
        // computing the center of the face
        let vects: Vector3<f64> = vertices_coordinates.iter().map(|v| v.coords).sum(); // this is so ugly
        let center: Point3<f64> =
            Point3::new(0., 0., 0.) + vects / vertices_coordinates.len() as f64; // I just want to average...
        Face {
            vertices: vertices_indices,
            center,
            area: {
                let mut area = 0.;
                let v1 = vertices_coordinates[1] - vertices_coordinates[0];
                for i in 2..vertices_coordinates.len() {
                    if v1
                        .dot(&(vertices_coordinates[i] - vertices_coordinates[0]))
                        .abs()
                        < f64::EPSILON
                    {
                        area = v1
                            .cross(&(vertices_coordinates[i] - vertices_coordinates[0]))
                            .norm()
                            .abs()
                    }
                }
                area
            },
            neighbor: NeighborType::Nil,
            n: Vector3::zeros(),
            distance: 0.,
        }
    }

    pub fn neighbour(&self) -> Option<usize> {
        if let NeighborType::Internal(i) = self.neighbor {
            Some(i)
        } else {
            None
        }
    }
}
