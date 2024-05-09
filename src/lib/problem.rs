use crate::computational_cell::{HasProperties, NeighborType};
use crate::equation;
use crate::equation_container::EquationContainer;
use crate::mesh::Mesh;
use nalgebra::{DVector, Point3, Vector3};

#[derive(Clone)]
pub struct Problem<T: HasProperties> {
    pub mesh: Mesh<T>,
    pub equation_containers: Vec<EquationContainer<T>>,
}

impl<T: HasProperties> Problem<T> {
    pub fn new() -> Self {
        Problem {
            mesh: Mesh::new(),
            equation_containers: Vec::new(),
        }
    }

    pub fn n_cells(&self) -> usize {
        self.mesh.cells.len()
    }

    pub fn n_eq(&self) -> usize {
        let mut n = 0;
        for layer in self.equation_containers.iter() {
            n += layer.n_eq();
        }
        n
    }

    /// Generates a rectangular parallelepipedal piece of a mesh that consists
    /// of a regular grid of rectangural parallelepipedal cells.
    pub fn add_zone(
        &mut self,
        // two points defining a rectangular parallelepiped
        r1: (f64, f64, f64),
        r2: (f64, f64, f64),
        // number of cells in each direction
        (nx, ny, nz): (usize, usize, usize),
        bc: [Option<Vec<equation::BoundaryCondition<T>>>; 6],
        equations: Vec<equation::Equation<T>>,
        internal_props: T,
    ) -> &mut Problem<T> {
        // TODO: I could not think of an elegant way to pass boundary conditions
        // in the current state, the order in which they are passed are very important!
        let origin = Point3::new(r1.0, r1.1, r1.2);
        let ending = Point3::new(r2.0, r2.1, r2.2);
        let ex = (ending - origin)[0] / nx as f64 * Vector3::new(1., 0., 0.);
        let ey = (ending - origin)[1] / ny as f64 * Vector3::new(0., 1., 0.);
        let ez = (ending - origin)[2] / nz as f64 * Vector3::new(0., 0., 1.);
        // make a list of the actually defined boundary conditions
        let boundary_order = [
            (origin, -ex),
            (origin, -ey),
            (origin, -ez),
            (ending, ex),
            (ending, ey),
            (ending, ez),
        ];
        let mut boundary_index = 0;
        let mut boundaries = Vec::new();
        for i in 0..6 {
            if let Some(bc) = &bc[i] {
                boundaries.push((
                    boundary_order[i].0,
                    boundary_order[i].1,
                    Some(bc.clone()),
                    boundary_index,
                ));
                boundary_index += 1;
            } else {
                boundaries.push((boundary_order[i].0, boundary_order[i].1, None, usize::MAX));
            }
        }
        // loop that creates one cell at a time
        let mut zone_cell_indices = Vec::new();
        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    let cell_index = self.mesh.add_hex_cell(
                        origin + ex * i as f64 + ey * j as f64 + ez * k as f64,
                        origin + ex * (i + 1) as f64 + ey * (j + 1) as f64 + ez * (k + 1) as f64,
                        internal_props.clone(),
                    );
                    zone_cell_indices.push(cell_index);
                    // check if at boundary. if so, create boundary condition
                    if i == 0 {
                        if let (point, normal, Some(_), index) = boundaries[0] {
                            // unwrap here because something is wrong if this fails
                            let mut boundary_face = self
                                .mesh
                                .get_face_on_plane(point, normal, cell_index)
                                .unwrap();
                            // we might need to do something more complex here in the future
                            boundary_face.neighbor = NeighborType::Boundary(index);
                        }
                    }
                    if j == 0 {
                        if let (point, normal, Some(_), index) = boundaries[1] {
                            // unwrap here because something is wrong if this fails
                            let mut boundary_face = self
                                .mesh
                                .get_face_on_plane(point, normal, cell_index)
                                .unwrap();
                            // we might need to do something more complex here in the future
                            boundary_face.neighbor = NeighborType::Boundary(index);
                        }
                    }
                    if k == 0 {
                        if let (point, normal, Some(_), index) = boundaries[2] {
                            // unwrap here because something is wrong if this fails
                            let mut boundary_face = self
                                .mesh
                                .get_face_on_plane(point, normal, cell_index)
                                .unwrap();
                            // we might need to do something more complex here in the future
                            boundary_face.neighbor = NeighborType::Boundary(index);
                        }
                    }
                    if i == nx - 1 {
                        if let (point, normal, Some(_), index) = boundaries[3] {
                            // unwrap here because something is wrong if this fails
                            let mut boundary_face = self
                                .mesh
                                .get_face_on_plane(point, normal, cell_index)
                                .unwrap();
                            // we might need to do something more complex here in the future
                            boundary_face.neighbor = NeighborType::Boundary(index);
                        }
                    }
                    if j == ny - 1 {
                        if let (point, normal, Some(_), index) = boundaries[4] {
                            // unwrap here because something is wrong if this fails
                            let mut boundary_face = self
                                .mesh
                                .get_face_on_plane(point, normal, cell_index)
                                .unwrap();
                            // we might need to do something more complex here in the future
                            boundary_face.neighbor = NeighborType::Boundary(index);
                        }
                    }
                    if k == nz - 1 {
                        if let (point, normal, Some(_), index) = boundaries[5] {
                            // unwrap here because something is wrong if this fails
                            let mut boundary_face = self
                                .mesh
                                .get_face_on_plane(point, normal, cell_index)
                                .unwrap();
                            // we might need to do something more complex here in the future
                            boundary_face.neighbor = NeighborType::Boundary(index);
                        }
                    }
                }
            }
        }
        let mut used_bc: Vec<Vec<equation::BoundaryCondition<T>>> = Vec::new();
        for b in boundaries {
            if let (_, _, Some(bb), _) = b {
                used_bc.push(bb);
            }
        }
        self.equation_containers.push(EquationContainer {
            cell_indices: zone_cell_indices,
            boundaries: used_bc,
            equations: equations,
        });
        self
    }

    pub fn get_state(&self) -> Vec<DVector<f64>> {
        // returns the current state of the mesh
        // TODO: Vec capacity should be calculated correctly for faster computation
        let mut x: Vec<DVector<f64>> = Vec::with_capacity(self.mesh.n_cells());
        for layer in self.equation_containers.iter() {
            for i_cell in layer.cell_indices.iter() {
                for eq in layer.equations.iter() {
                    x.push((eq.get)(&self.mesh.cells[*i_cell].properties));
                }
            }
        }
        x
    }

    pub fn setup_equations(&self, _t: &f64, _x: &[DVector<f64>]) -> Vec<DVector<f64>> {
        let mut y = vec![];
        for layer in self.equation_containers.iter() {
            for i_cell in layer.cell_indices.iter() {
                for eq in layer.equations.iter() {
                    y.push((eq.rhs)(&self.mesh.cells, i_cell))
                }
            }
        }
        y
    }

    pub fn update(&mut self, x: &[DVector<f64>]) {
        // takes the solution and updates every cell with it
        let mut i = 0;
        for layer in self.equation_containers.iter() {
            for i_cell in layer.cell_indices.iter() {
                for eq in layer.equations.iter() {
                    (eq.lhs)(&mut self.mesh.cells[*i_cell].properties, &x[i]);
                    i += 1;
                }
                self.mesh.cells[*i_cell].properties.update()
            }
        }
        self.enforce_boundary_conditions();
    }

    pub fn enforce_boundary_conditions(&mut self) {
        // for cell in self.cells.iter() {
        //     let cellref = &mut cell.borrow_mut();
        //     if let Boundary(bc_vector) = &cellref.variant {
        //         let mut newproperties = cellref.properties;
        //         for BoundaryCondition(function, value) in bc_vector {
        //             function(&mut newproperties, &cellref.faces, value);
        //         }
        //         cellref.properties = newproperties;
        //     }
        // }
    }
}
