use std::path::PathBuf;

use crate::computational_cell::{ComputationalCell, Face, HasProperties, NeighborType};
use crate::geometry_utils::polyhedron_center;
use nalgebra::{Point3, Vector3};
use vtkio::model::{Attributes, DataArray, ElementType};

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

    pub fn add_hex_cell(
        &mut self,
        r1: Point3<f64>,
        r2: Point3<f64>,
        properties: T,
        zone: usize,
    ) -> usize {
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
            for (index, face) in self.faces.iter_mut().enumerate() {
                // WARNING : this assumes that no crazy face with a different number of vertices exists
                if face.vertices.iter().all(|p| vertices.contains(p)) {
                    cell_faces.push(index);
                    let face_distance = (face.center - cell_center).norm();
                    face.neighbor = NeighborType::Internal(cell_index, face_distance);
                    found = true;
                    break;
                }
            }
            if !found {
                let mut face = Face::new(self, vertices);
                let face_distance = (face.center - cell_center).norm();
                // face normal points from the owner to the neighbour
                face.owner = (cell_index, face_distance);
                cell_faces.push(self.faces.len());
                self.faces.push(face);
            }
        }
        // check if the faces are correctly oriented, meaning that the face normals point outwards of the owner and inwards of the neighbour
        for i_face in cell_faces.iter() {
            let face = &self.faces[*i_face];
            let to_face = face.center - cell_center;
            if cell_index == face.owner.0 {
                if to_face.dot(&face.normal) < 0. {
                    panic!("Face normal points into the owner cell")
                }
            } else if to_face.dot(&face.normal) > 0. {
                panic!("Face normal points into the neighbor cell")
            }
        }

        self.cells
            .push(ComputationalCell::new(self, cell_faces, properties, zone));
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

    /// Exports the mesh to the VTK format.
    /// vtkio crate support for polyhedra is not yet implemented.
    /// This function will only work for hexahedra.
    pub fn export_to_vtk(&self, filename: &str) {
        // collect the lists of points that define the cells
        let cell_indices = self
            .cells
            .iter()
            .map(|cell| {
                cell.hexahedron_vertices(self)
                    .iter()
                    .map(|i| *i as u64)
                    .collect::<Vec<u64>>()
            })
            .collect::<Vec<Vec<u64>>>();
        let offsets: Vec<u64> = cell_indices
            .iter()
            .scan(0, |acc, x| {
                *acc += x.len() as u64;
                Some(*acc)
            })
            .collect::<Vec<u64>>();
        // collect the data of the cells
        let cell_data = self
            .cells
            .iter()
            .map(|cell| cell.properties.export())
            .collect::<Vec<_>>();
        // data is obtained as a list of lists, where the "parent" list is the cell index and the "child" list is the data values
        // I want to transpose this list of lists to have one "parent" list per variable with the "child" list being the values of the variable in all cells
        // TODO : here I assume that all cells posess the same properties.
        // I want this to be more generic at some point, although I don't know if it is compatible with the VTK format.
        let data_labels = cell_data[0]
            .iter()
            .map(|(label, _)| label.clone())
            .collect::<Vec<_>>();
        let data_values: Vec<Vec<f64>> = (0..data_labels.len())
            .map(|i| cell_data.iter().map(|data| data[i].1).collect::<Vec<f64>>())
            .collect::<Vec<Vec<f64>>>();

        // use the vtkio crate to write a VTK file
        let vtk_file = vtkio::Vtk {
            version: vtkio::model::Version { major: 4, minor: 2 },
            title: "Mesh".to_string(),
            byte_order: vtkio::model::ByteOrder::BigEndian,
            file_path: None,
            data: vtkio::model::DataSet::inline(vtkio::model::UnstructuredGridPiece {
                points: vtkio::IOBuffer::F64(
                    // vec that contains all point coordinates
                    self.vertices
                        .iter()
                        .flat_map(|p| vec![p.x, p.y, p.z])
                        .collect::<Vec<f64>>(),
                ),
                cells: vtkio::model::Cells {
                    cell_verts: vtkio::model::VertexNumbers::XML {
                        connectivity: cell_indices.concat(),
                        offsets,
                    },
                    types: cell_indices
                        .iter()
                        .map(|index_list| match index_list.len() {
                            8 => vtkio::model::CellType::Hexahedron,
                            _ => vtkio::model::CellType::Polyhedron,
                        })
                        .collect(),
                },
                data: Attributes {
                    cell: {
                        let mut attributes_items = Vec::new();
                        for (label, values) in data_labels.iter().zip(data_values.iter()) {
                            attributes_items.push(vtkio::model::Attribute::DataArray(DataArray {
                                name: label.clone(),
                                elem: ElementType::Scalars {
                                    num_comp: 1,
                                    lookup_table: None,
                                },
                                data: vtkio::IOBuffer::F64(values.clone()),
                            }));
                        }
                        attributes_items
                    },
                    ..Default::default()
                },
            }),
        };

        // export
        let file_path = PathBuf::from(filename);
        vtk_file
            .export_ascii(file_path)
            .expect("Failed to export mesh to VTK");
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

    #[allow(non_snake_case)]
    #[derive(Default, Copy, Clone, Debug)]
    pub struct Properties {
        T: f64,
        h: f64,
    }
    impl HasProperties for Properties {
        fn update(&mut self) {}
        fn export(&self) -> Vec<(String, f64)> {
            vec![("T".to_string(), self.T), ("h".to_string(), self.h)]
        }
    }

    #[test]
    fn single_cell_mesh() {
        let mut mesh: crate::mesh::Mesh<Properties> = crate::mesh::Mesh::new();
        mesh.add_hex_cell(
            nalgebra::Point3::new(0., 0., 0.),
            nalgebra::Point3::new(1., 1., 1.),
            Properties::default(),
            0,
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

    #[test]
    fn export_single_cell_mesh() {
        let mut mesh: crate::mesh::Mesh<Properties> = crate::mesh::Mesh::new();
        mesh.add_hex_cell(
            nalgebra::Point3::new(0., 0., 0.),
            nalgebra::Point3::new(1., 1., 1.),
            Properties::default(),
            0,
        );
        mesh.export_to_vtk("single_cell_mesh.vtk");
    }

    #[test]
    fn export_double_cell_mesh() {
        let mut mesh: crate::mesh::Mesh<Properties> = crate::mesh::Mesh::new();
        mesh.add_hex_cell(
            nalgebra::Point3::new(0., 0., 0.),
            nalgebra::Point3::new(1., 1., 1.),
            Properties::default(),
            0,
        );
        mesh.add_hex_cell(
            nalgebra::Point3::new(1., 0., 0.),
            nalgebra::Point3::new(2., 1., 1.),
            Properties::default(),
            0,
        );
        mesh.export_to_vtk("double_cell_mesh.vtk");
    }

    #[test]
    fn print_single_cell_vertices() {
        let mut mesh: crate::mesh::Mesh<Properties> = crate::mesh::Mesh::new();
        mesh.add_hex_cell(
            nalgebra::Point3::new(0., 0., 0.),
            nalgebra::Point3::new(1., 1., 1.),
            Properties::default(),
            0,
        );
        println!("Printing vertices of all faces of the single cell mesh:");
        for i_face in mesh.cells[0].faces.iter() {
            let face = &mesh.faces[*i_face];
            println!("Face vertices:");
            for i_vertex in face.vertices.iter() {
                println!("{:?}", mesh.vertices[*i_vertex]);
            }
        }
        println!("Printing the vertices given by the function:");
        for i_vertex in mesh.cells[0].hexahedron_vertices(&mesh) {
            println!("{:?}", mesh.vertices[i_vertex]);
        }
    }
}
