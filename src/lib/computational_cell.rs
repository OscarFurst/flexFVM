// for now all cells are rectangular parallelepipeds
use crate::mesh::Mesh;
use crate::{geometry_utils::*, problem};
use nalgebra::{Point3, RowDVector, Vector3};
use nalgebra_sparse::csr::CsrMatrix;
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
pub struct CellData {
    pub faces: Vec<Vec<usize>>,
    pub r: Vec<Point3<f64>>,
    pub volume: RowDVector<f64>,
    pub zone: Vec<usize>,
}
impl Default for CellData {
    fn default() -> Self {
        CellData {
            faces: Vec::new(),
            r: Vec::new(),
            volume: RowDVector::from_vec(Vec::new()),
            zone: Vec::new(),
        }
    }
}
impl CellData {
    pub fn push(&mut self, faces: Vec<usize>, r: Point3<f64>, volume: f64, zone: usize) {
        self.faces.push(faces);
        self.r.push(r);
        self.volume = self.volume.transpose().push(volume).transpose();
        self.zone.push(zone);
    }

    pub fn get_cell(&self, i_cell: usize) -> ComputationalCell {
        ComputationalCell {
            faces: &self.faces[i_cell],
            r: &self.r[i_cell],
            volume: &self.volume[i_cell],
            zone: &self.zone[i_cell],
        }
    }

    pub fn get_cell_mut(&mut self, i_cell: usize) -> ComputationalCellMut {
        ComputationalCellMut {
            faces: &self.faces[i_cell],
            r: &mut self.r[i_cell],
            volume: &mut self.volume[i_cell],
            zone: &mut self.zone[i_cell],
        }
    }
}

pub struct ComputationalCell<'a> {
    // faces of the volume (defines it)
    pub faces: &'a [usize],
    // geometric properties
    pub r: &'a Point3<f64>, // m
    pub volume: &'a f64,    // m3
    pub zone: &'a usize,
}

pub struct ComputationalCellMut<'a> {
    // faces of the volume (defines it)
    pub faces: &'a [usize],
    // geometric properties
    pub r: &'a mut Point3<f64>, // m
    pub volume: &'a mut f64,    // m3
    pub zone: &'a mut usize,
}

pub trait HasProperties: Clone {
    /// new is called whenever a new set of properties is needed.
    /// It needs to be able to create an empty set of properties and a set of properties with a given size.
    fn new(size: usize) -> Self;
    /// update_primitive is called after the new solution variables are computed.
    /// It is meant to update the properties that are directly derived from the solution variables
    fn update_primitive(&mut self, _volume: &RowDVector<f64>) {}
    /// update_integral is called after the mesh is generated
    /// It is meant to update the properties that are derived from the geometry of the cell,
    /// so that densities can be used in the problem initialization while the mass is the actual
    /// variable that is solved for.
    fn update_integral(&mut self, _volume: &RowDVector<f64>) {}
    /// interpolate_face_from_cells is called after update_primitive.
    /// This is optional, as the interpolation can be done directly in the computatino of the derivatives,
    /// although that will be slower in most cases!
    /// It assigns to every side of every face the value determined by the selected interpolation scheme.
    /// Requires to select an interpolation scheme by calling an interpolation initialization function.
    fn interpolate_face_from_cells(
        &mut self,
        _cell_properties: &Self,
        _interpolation_matrix: &CsrMatrix<f64>,
    ) {
    }
    /// export is called when the solution is exported to a file.
    // I would love to use a data type that does not require a copy of the data, but vtkio requires a Vec \_(ツ)_/
    fn export(&self) -> Vec<(String, Vec<f64>)>;
}

impl<'a> ComputationalCell<'a> {
    /// Returns a list of vertex indices which is correctly sorted for VTK.
    /// Only works for hexahedrons.
    pub fn hexahedron_vertices(&self, mesh: &Mesh) -> Vec<usize> {
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
    fn are_connected(&self, mesh: &Mesh, i_vertex1: usize, i_vertex2: usize) -> bool {
        let faces1: Vec<usize> = self
            .faces
            .iter()
            .map(|i_face| &mesh.faces[*i_face])
            .enumerate()
            .filter(|(_, face)| face.vertices.contains(&i_vertex1))
            .map(|(i, _)| i)
            .collect();
        let faces2: Vec<usize> = self
            .faces
            .iter()
            .map(|i_face| &mesh.faces[*i_face])
            .enumerate()
            .filter(|(_, face)| face.vertices.contains(&i_vertex2))
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

impl<'a> fmt::Debug for ComputationalCell<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.r.fmt(f)
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
    /// Creates a new face from a list of point indices and initialize its basic geometrical properties.
    /// The properties relating to the neighbour need to be initialized later.
    pub fn new(mesh: &Mesh, vertices_indices: &[usize]) -> Self {
        let vertices_coordinates: &[Point3<f64>] = &vertices_indices
            .iter()
            .map(|i_vertex| mesh.vertices[*i_vertex])
            .collect::<Vec<Point3<f64>>>();
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

    pub fn is_owner(&self, i_cell: usize) -> bool {
        self.owner.0 == i_cell
    }

    pub fn out_normal(&self, i_cell: usize) -> Vector3<f64> {
        if self.is_owner(i_cell) {
            self.normal
        } else {
            -self.normal
        }
    }

    pub fn get_zone<T: HasProperties>(&self, problem: &problem::Problem<T>) -> usize {
        problem.mesh.cells.zone[self.owner.0]
    }
}
