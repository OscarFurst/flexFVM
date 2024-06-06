use crate::computational_cell::{HasProperties, NeighborType};
use crate::equation::{self, BoundaryCondition};
use crate::equation_container::EquationContainer;
use crate::geometry_utils::points_are_collinear;
use crate::mesh::Mesh;
use core::panic;
use indicatif::ProgressBar;
use nalgebra::{DVector, Point3, Vector3};
use nalgebra_sparse::{coo::CooMatrix, csr::CsrMatrix};
use std::error::Error;
use std::fs::File;
use std::path::PathBuf;
use vtkio::model::{Attributes, DataArray, ElementType};

pub enum NeighborReturn<'a> {
    Internal(usize, usize, f64),
    Boundary(&'a str, usize, f64, &'a DVector<f64>),
    Nil, // should only exist in 1D or 2D cases, where there is no neighbor in the missing dimension
}

#[derive(Clone)]
pub struct Problem<T: HasProperties> {
    pub mesh: Mesh,
    pub cell_properties: T,
    pub face_owner_properties: T,
    pub face_owner_interpolation: Option<CsrMatrix<f64>>,
    pub face_neighbor_properties: T,
    pub face_neighbor_interpolation: Option<CsrMatrix<f64>>,
    pub equation_containers: Vec<EquationContainer<T>>,
    // 2D (i_cell, i_eq) array which provides the starting index of equation i_eq of cell i_cell
    equation_positions: Vec<Vec<usize>>,
}

impl<T: HasProperties> Problem<T> {
    pub fn new() -> Self {
        Problem {
            mesh: Mesh::new(),
            cell_properties: T::new(0),
            face_owner_properties: T::new(0),
            face_owner_interpolation: Option::None,
            face_neighbor_properties: T::new(0),
            face_neighbor_interpolation: Option::None,
            equation_containers: Vec::new(),
            equation_positions: Vec::new(),
        }
    }

    pub fn n_cells(&self) -> usize {
        self.mesh.cells.volume.len()
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
        bc: [Option<BoundaryCondition>; 6],
        equations: Vec<equation::Equation<T>>,
        internal_props: fn(&mut T, &Point3<f64>),
    ) -> &mut Problem<T> {
        let zone_id = self.equation_containers.len();
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
        // equation container is already set here because it is needed to set up the boundary conditions
        let mut used_bc: Vec<BoundaryCondition> = Vec::new();
        for b in boundaries.clone() {
            if let (_, _, Some(bb), _) = b {
                used_bc.push(bb);
            }
        }
        self.equation_containers.push(EquationContainer {
            cell_indices: vec![],
            equations,
            boundary_conditions: used_bc,
        });
        // loop that creates one cell at a time
        println!("🔗 Meshing ...");
        let pb = ProgressBar::new((nx) as u64);
        pb.set_style(
            indicatif::ProgressStyle::default_bar()
                .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7}")
                .unwrap()
                .progress_chars("##-"),
        );
        let mut zone_cell_indices = Vec::new();
        for i in 0..nx {
            pb.inc(1);
            for j in 0..ny {
                for k in 0..nz {
                    let cell_origin = origin + ex * i as f64 + ey * j as f64 + ez * k as f64;
                    let cell_ending =
                        origin + ex * (i + 1) as f64 + ey * (j + 1) as f64 + ez * (k + 1) as f64;
                    let cell_center = cell_origin + (cell_ending - cell_origin) / 2.0;
                    let cell_index = self.mesh.add_hex_cell(cell_origin, cell_ending, zone_id);
                    internal_props(&mut self.cell_properties, &cell_center);
                    zone_cell_indices.push(cell_index);
                    // check if at boundary. if so, create boundary condition
                    if i == 0 {
                        self.setup_boundary_contidions(cell_index, &boundaries[0]);
                    }
                    if j == 0 {
                        self.setup_boundary_contidions(cell_index, &boundaries[1]);
                    }
                    if k == 0 {
                        self.setup_boundary_contidions(cell_index, &boundaries[2]);
                    }
                    if i == nx - 1 {
                        self.setup_boundary_contidions(cell_index, &boundaries[3]);
                    }
                    if j == ny - 1 {
                        self.setup_boundary_contidions(cell_index, &boundaries[4]);
                    }
                    if k == nz - 1 {
                        self.setup_boundary_contidions(cell_index, &boundaries[5]);
                    }
                }
            }
        }
        pb.finish();
        self.equation_containers[zone_id].cell_indices = zone_cell_indices;
        self
    }

    fn setup_boundary_contidions(
        &mut self,
        cell_index: usize,
        boundary: &(Point3<f64>, Vector3<f64>, Option<BoundaryCondition>, usize),
    ) {
        // TODO : this function is cursed, needs to be cleaned up.
        if let (point, normal, Some(BoundaryCondition(name, values)), index) = boundary {
            // unwrap here because something is wrong if this fails
            let face_id = self
                .mesh
                .get_face_id_on_plane(point, normal, cell_index)
                .unwrap();
            let face_center = self.mesh.faces.center[face_id];
            let mut neighb = NeighborType::Nil;
            match name as &str {
                "periodic" => {
                    // search for the other face with the same boundary condition
                    let mut found_twin = false;
                    if self.mesh.periodic_unpaired.len() > 0 {
                        let twin_direction = Vector3::new(values[0], values[1], values[2]);
                        for &i_f in self.mesh.periodic_unpaired.clone().iter() {
                            // we found another face with a periodic BC, now we need to check if it is the twin face
                            let face = self.mesh.faces.get_face_mut(i_f);
                            if (*face.center - face_center).cross(&twin_direction).norm()
                                < f64::EPSILON
                            {
                                found_twin = true;
                                self.mesh.periodic_unpaired.remove(&i_f);
                                self.mesh.periodic_pairs.insert(face_id, i_f);
                                self.mesh.periodic_pairs.insert(i_f, face_id);
                                // in this case we have to set up two faces instead of jsut one !
                                *face.neighbor = NeighborType::Internal(
                                    cell_index,
                                    (*face.center - self.mesh.cells.r[cell_index]).norm(),
                                );
                                neighb = NeighborType::Internal(
                                    face.owner.0,
                                    (self.mesh.cells.r[face.owner.0] - *face.center).norm(),
                                );
                                break;
                            }
                        }
                    }
                    if !found_twin {
                        self.mesh.periodic_unpaired.insert(face_id);
                        neighb = NeighborType::Boundary(*index);
                    }
                }
                _ => {
                    // for every other boundary condition, the BC index is enough
                    neighb = NeighborType::Boundary(*index);
                }
            }
            let boundary_face = self
                .mesh
                .get_face_on_plane(point, normal, cell_index)
                .unwrap();
            *boundary_face.neighbor = neighb;
            // make sure a boundary was assigned correctly
            if matches!(boundary_face.neighbor, NeighborType::Nil) {
                panic!("Boundary condition not assigned correctly");
            }
        }
    }

    /// Checks if the mesh boundary condition was properly set up.
    /// Since the boundary condition requires a 1-to-1 pairing of faces,
    /// no face should have a NeighborType::Boundary(index) with the index
    /// of a periodic boundary condition.
    pub fn check_periodic_bc(&self) -> bool {
        if self.mesh.periodic_unpaired.len() > 0 {
            println!("Some periodic boundary conditions were not paired.");
            return false;
        }
        for i_f in 0..self.mesh.faces.len() {
            if let NeighborType::Boundary(bc_index) = self.mesh.faces.neighbor[i_f] {
                if self.equation_containers[self.get_face_zone(i_f)].boundary_conditions[bc_index].0
                    == "periodic"
                {
                    return false;
                }
            }
        }
        true
    }

    /// Checks which zone of the mesh this face is in.
    pub fn get_face_zone(&self, i_face: usize) -> usize {
        self.mesh.cells.zone[self.mesh.faces.owner[i_face].0]
    }

    pub fn sum_conservation_variables(&self) -> Vec<DVector<f64>> {
        let mut sums: Vec<DVector<f64>> = Vec::new();
        for layer in self.equation_containers.iter() {
            for eq in layer.equations.iter() {
                sums.push(
                    layer
                        .cell_indices
                        .iter()
                        .map(|&i| (eq.get)(&self.cell_properties, i))
                        .sum(),
                );
            }
        }
        sums
    }

    pub fn get_state(&self) -> DVector<f64> {
        // returns the current state of the mesh
        let mut x: DVector<f64> = DVector::zeros(self.n_eq());
        let mut i = 0;
        for layer in self.equation_containers.iter() {
            for i_cell in layer.cell_indices.iter() {
                for eq in layer.equations.iter() {
                    x.rows_mut(i, eq.dimension)
                        .copy_from(&(eq.get)(&self.cell_properties, *i_cell));
                    i += eq.dimension;
                }
            }
        }
        x
    }

    /// Writes where to find the index of the first entry of equation i_eq for cell i_cell in the state vector.
    /// The result is written to "equation_positions".
    /// This is used to speed up the computation of the equations when using face fluxes.
    // TODO : this should not be public, it should be called by the problem itself at some point
    pub fn set_equation_positions(&mut self) {
        self.equation_positions.resize(self.n_cells(), Vec::new());
        let mut i = 0;
        for layer in self.equation_containers.iter() {
            for &i_cell in layer.cell_indices.iter() {
                for eq in layer.equations.iter() {
                    self.equation_positions[i_cell].push(i);
                    i += eq.dimension;
                }
            }
        }
    }

    pub fn setup_equations_in_place(&self, _t: f64, y: &mut DVector<f64>) {
        assert_eq!(
            y.len(),
            self.n_eq(),
            "Dimension mismatch in setup_equations_in_place"
        );
        if self.equation_positions.is_empty() {
            panic!("Equation positions not set up");
        }
        y.fill(0.0);
        // TODO : disable the "surface flux equations" completely if they are all None so that it's faster
        // TODO : the following lines are a test, to be improved!!
        let max_equation_dimension = self
            .equation_containers
            .iter()
            .map(|layer| layer.equations.iter().map(|eq| eq.dimension).sum::<usize>())
            .max()
            .unwrap();
        let mut result_buffer = DVector::zeros(max_equation_dimension);
        // end of test, to be improved
        for i_f in 0..self.mesh.faces.len() {
            for (i_eq, eq) in self.equation_containers[self.get_face_zone(i_f)]
                .equations
                .iter()
                .enumerate()
            {
                if let Some(eqn) = eq.flux {
                    result_buffer.fill(0.0);
                    (eqn)(self, i_f, &mut result_buffer);
                    // I really don't know what the most performant way to do this is.
                    // TODO : this looks terrible.
                    let i_o = self.mesh.faces.owner[i_f].0;
                    if let NeighborType::Internal(i_n, _) = self.mesh.faces.neighbor[i_f] {
                        // the "contains_key" is slower than I would like it to be...
                        if self.mesh.periodic_pairs.contains_key(&i_f) {
                            for i in 0..eq.dimension {
                                y[self.equation_positions[i_o][i_eq] + i] += result_buffer[i];
                            }
                        } else {
                            for i in 0..eq.dimension {
                                y[self.equation_positions[i_o][i_eq] + i] += result_buffer[i];
                                y[self.equation_positions[i_n][i_eq] + i] -= result_buffer[i];
                            }
                        }
                    } else {
                        for i in 0..eq.dimension {
                            y[self.equation_positions[i_o][i_eq] + i] += result_buffer[i];
                        }
                    }
                }
            }
        }
        let mut i = 0;
        for layer in self.equation_containers.iter() {
            for &i_cell in layer.cell_indices.iter() {
                for eq in layer.equations.iter() {
                    if let Some(eqn) = eq.source {
                        // get the source terms
                        result_buffer.fill(0.0);
                        (eqn)(self, i_cell, &mut result_buffer);
                        // put them in the right part of the derivative vector
                        for ii in 0..eq.dimension {
                            y[i + ii] += result_buffer[ii];
                        }
                    }
                    i += eq.dimension;
                }
            }
        }
    }

    pub fn setup_equations(&self, _t: f64) -> DVector<f64> {
        let mut y: DVector<f64> = DVector::zeros(self.n_eq());
        self.setup_equations_in_place(_t, &mut y);
        y
    }

    pub fn update(&mut self, x: &DVector<f64>) {
        // takes the solution and updates every cell with it
        let mut i = 0;
        for layer in self.equation_containers.iter() {
            for i_cell in layer.cell_indices.iter() {
                for eq in layer.equations.iter() {
                    (eq.lhs)(
                        &mut self.cell_properties,
                        *i_cell,
                        &x.as_slice()[i..i + eq.dimension],
                    );
                    i += eq.dimension;
                }
            }
        }
        self.cell_properties
            .update_primitive(&self.mesh.cells.volume);
        if let Some(interpolation) = &self.face_owner_interpolation {
            self.face_owner_properties
                .interpolate_face_from_cells(&self.cell_properties, interpolation);
        }
        if let Some(interpolation) = &self.face_neighbor_interpolation {
            self.face_neighbor_properties
                .interpolate_face_from_cells(&self.cell_properties, interpolation);
        }
    }

    pub fn update_integral(&mut self) {
        // takes the solution and updates every cell with it
        self.cell_properties
            .update_integral(&self.mesh.cells.volume);
    }

    pub fn initialize_face_properties(&mut self) {
        // initializes the face properties
        self.face_owner_properties = T::new(self.mesh.faces.len());
        self.face_neighbor_properties = T::new(self.mesh.faces.len());
    }

    /// Prepares the interpolation of the face properties using a constant interpolation scheme.
    /// A sparse matrix A is created so that X_face = A * X_cell.
    /// This has to be done only once, before the simulation starts.
    // TODO : needs its own module
    #[allow(non_snake_case)]
    pub fn prepare_constant_interpolation(&mut self) {
        let mut A_owner = CooMatrix::new(self.mesh.faces.len(), self.mesh.cells.volume.len());
        let mut A_neighbor = CooMatrix::new(self.mesh.faces.len(), self.mesh.cells.volume.len());
        for i_face in 0..self.mesh.faces.len() {
            let face = self.mesh.faces.get_face(i_face);
            let owner = face.owner.0;
            A_owner.push(i_face, owner, 1.0);
            // TODO : seems like a great place to put Dirichlet boundary conditions
            match *face.neighbor {
                NeighborType::Internal(neighbor, _) => A_neighbor.push(i_face, neighbor, 1.0),
                _ => {}
            };
        }
        self.face_owner_interpolation = Some(CsrMatrix::from(&A_owner));
        self.face_neighbor_interpolation = Some(CsrMatrix::from(&A_neighbor));
        self.initialize_face_properties();
    }

    /// Prepares the interpolation of the face properties using a linear interpolation scheme.
    /// A sparse matrix A is created so that X_face = A * X_cell.
    /// This has to be done only once, before the simulation starts.
    #[allow(non_snake_case)]
    pub fn prepare_linear_interpolation(&mut self) {
        let mut A = CooMatrix::new(self.mesh.faces.len(), self.mesh.cells.volume.len());
        for i_face in 0..self.mesh.faces.len() {
            let face = &self.mesh.faces.get_face(i_face);
            let owner = face.owner;
            let neighbor = match *face.neighbor {
                NeighborType::Internal(neighbor, distance) => (neighbor, distance),
                _ => (owner.0, owner.1),
            };
            A.push(i_face, owner.0, owner.1 / (owner.1 + neighbor.1));
            A.push(i_face, neighbor.0, neighbor.1 / (owner.1 + neighbor.1));
        }
        self.face_owner_interpolation = Some(CsrMatrix::from(&A));
        self.face_neighbor_interpolation = Some(CsrMatrix::from(&A));
        self.initialize_face_properties();
    }

    /// Prepares the interpolation of the face properties using a second order MUSCL interpolation scheme.
    /// A sparse matrix A is created so that X_face = A * X_cell.
    /// This has to be done only once, before the simulation starts.
    /// Method is from the following paper :
    /// Clément Le Touze, Angelo Murrone, Herve Guillard. Multislope MUSCL method for general unstructured
    /// meshes. Journal of Computational Physics, 2014, pp.44. 10.1016/j.jcp.2014.12.032 . hal-00939475v3
    ///
    /// Currently, 2D is only supported for structured meshes.
    #[allow(non_snake_case)]
    pub fn prepare_MUSCL_interpolation(&mut self) {
        let mut A_owner = CooMatrix::new(self.mesh.faces.len(), self.mesh.cells.volume.len());
        let mut A_neighbor = CooMatrix::new(self.mesh.faces.len(), self.mesh.cells.volume.len());
        // Allocate enough space. This is a bit more than required due to boundaries which take 1 entry instead of 7,
        // but the matrix is deallocated at the end of the function anyway.
        A_owner.reserve(self.mesh.faces.len() * 7);
        A_neighbor.reserve(self.mesh.faces.len() * 7);
        // make a nice progress bar
        println!("💪 Preparing MUSCL face-values interpolation ...");
        let pb: ProgressBar = ProgressBar::new(self.mesh.cells.r.len() as u64);
        pb.set_style(
            indicatif::ProgressStyle::default_bar()
                .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}")
                .unwrap()
                .progress_chars("##-"),
        );
        let mut good_faces = 0;
        let total_faces = self.mesh.faces.len();
        // this is a preparation for the algorithm that will follow
        self.mesh.collect_vertex_cells();
        for i_p in 0..self.mesh.cells.r.len() {
            // update progress bar
            pb.inc(1);
            pb.set_message(format!(
                "second order faces : {} / {}",
                good_faces, total_faces
            ));
            // first we look for all the cells sharing at least one vertex with the cell i_cell
            // TODO : this does not account for periodic boundary conditions in which the cells share neither vertices nor faces!!!
            let mut neighbors = self.mesh.extended_neighborhood(i_p);
            // at this point we have a list of neighbors.
            for &i_face in self.mesh.cells.faces[i_p].iter() {
                // now for each face of the cell, we want to find two points (forward and backward) so that we can compute a slope
                // the points are interpolated from a list of 3 cell centers and barycentric coordinates
                let mut point_forward = None;
                let mut point_backward = None;
                // vector pointing from the owner center to the face center
                let Bi = &self.mesh.cells.r[i_p];
                let Mij = &self.mesh.faces.center[i_face];
                let BiMij = Mij - Bi;
                // now we sort the list of neighbors in decreasing order of the value cos(BkBi, BiMij)
                neighbors.sort_by(|&i, &j| {
                    let BkMij = self.mesh.cells.r[i] - Bi;
                    let cos_i = BiMij.dot(&BkMij) / (BiMij.norm() * BkMij.norm());
                    let BkMij = self.mesh.cells.r[j] - Bi;
                    let cos_j = BiMij.dot(&BkMij) / (BiMij.norm() * BkMij.norm());
                    cos_j.partial_cmp(&cos_i).expect("Failed to compare floats when sorting cell neighbors during MUSCL interpolation")
                });
                // in this sorted list, the first elements are the most backward neighbors, and the last elements are the most forward neighbors
                // find the number of elemets for which cos(BkBi, BiMij) is positive
                let n_back = neighbors
                    .iter()
                    .take_while(|&&i| {
                        let BkMij = self.mesh.cells.r[i] - Bi;
                        BiMij.dot(&BkMij) > 0.0
                    })
                    .count();
                if points_are_collinear(Mij, Bi, &self.mesh.cells.r[neighbors[0]]) {
                    point_backward = Some([(neighbors[0], 1.0), (0, 0.0), (0, 0.0)]);
                } else if n_back >= 3 {
                    // we try combinations of 3 neighbors sorted after the sum of their indices
                    point_backward = self
                        .mesh
                        .best_barycentric_interpolation(Mij, Bi, &neighbors);
                }
                // find the number of elemets for which cos(BkBi, BiMij) is positive
                neighbors.reverse();
                let n_forward = neighbors.len() - n_back;
                if points_are_collinear(Mij, Bi, &self.mesh.cells.r[neighbors[0]]) {
                    point_forward = Some([(neighbors[0], 1.0), (0, 0.0), (0, 0.0)]);
                } else if n_forward >= 3 {
                    // we try combinations of 3 neighbors sorted after the sum of their indices
                    point_forward = self
                        .mesh
                        .best_barycentric_interpolation(Mij, Bi, &neighbors);
                }
                // check if we are setting the "owner" side of the face or the "neighbor" side.
                let matrix = if self.mesh.faces.owner[i_face].0 == i_p {
                    &mut A_owner
                } else {
                    &mut A_neighbor
                };
                if let (
                    Some([(F1, Falpha), (F2, Fbeta), (F3, Fgamma)]),
                    Some([(B1, Balpha), (B2, Bbeta), (B3, Bgamma)]),
                ) = (point_forward, point_backward)
                {
                    good_faces += 1;
                    let origin = Point3::origin();
                    // reconstruct the position of the forward and backward points
                    // this looks stupid because we need to transform points into vectors
                    let F = origin
                        + Falpha * (self.mesh.cells.r[F1] - origin)
                        + Fbeta * (self.mesh.cells.r[F2] - origin)
                        + Fgamma * (self.mesh.cells.r[F3] - origin);
                    let B = origin
                        + Balpha * (self.mesh.cells.r[B1] - origin)
                        + Bbeta * (self.mesh.cells.r[B2] - origin)
                        + Bgamma * (self.mesh.cells.r[B3] - origin);
                    // get the geometrical data needed for the interpolation
                    let distance = (F - B).norm();
                    let to_face = (Mij - Bi).norm();
                    // the value at the face is the value at the center times the slope times to_face.
                    // the slope is the difference between the forward and backward values divided by the distance.
                    matrix.push(i_face, i_p, 1.0);
                    matrix.push(i_face, F1, -to_face * Falpha / distance);
                    matrix.push(i_face, F2, -to_face * Fbeta / distance);
                    matrix.push(i_face, F3, -to_face * Fgamma / distance);
                    matrix.push(i_face, B1, to_face * Balpha / distance);
                    matrix.push(i_face, B2, to_face * Bbeta / distance);
                    matrix.push(i_face, B3, to_face * Bgamma / distance);
                } else {
                    // first order interpolation
                    matrix.push(i_face, i_p, 1.0);
                }
            }
        }
        pb.finish();
        // TODO : use a specialized enum instead of an Option
        self.face_owner_interpolation = Some(CsrMatrix::from(&A_owner));
        self.face_neighbor_interpolation = Some(CsrMatrix::from(&A_neighbor));
        // TODO : this is a workaround (first order) beacuse I don't want to implement periodic BC correctly yet
        for i_face in 0..self.mesh.faces.len() {
            if let NeighborType::Internal(i_n, _) = self.mesh.faces.neighbor[i_face] {
                if let Some(neighbor_matrix) = &mut self.face_neighbor_interpolation {
                    if neighbor_matrix.row(i_face).nnz() == 0 {
                        // this side of the face is currently not interpolated despite connecting to a neighbor!
                        // we perform a first order interpolation from the neighbor
                        A_neighbor.push(i_face, i_n, 1.0);
                    }
                }
            }
        }
        self.face_neighbor_interpolation = Some(CsrMatrix::from(&A_neighbor));
        self.initialize_face_properties();
    }

    /// Returns all neighbors of the specified cell
    /// A neighbor can be :
    /// - an internal cell (index, distance between centers)
    /// - a boundary condition (name, some values)
    /// - Nil (no neighbor, should only be encountered in 2D or 1D meshes)
    pub fn neighbors(&self, i_cell: usize) -> Vec<NeighborReturn> {
        self.mesh.cells.faces[i_cell]
            .iter()
            .map(|&i_face| match self.mesh.faces.neighbor[i_face] {
                NeighborType::Internal(i, n_distance) if i == i_cell => {
                    let (i, o_distance) = self.mesh.faces.owner[i_face];
                    NeighborReturn::Internal(i, i_face, o_distance + n_distance)
                }
                NeighborType::Internal(i, n_distance) => NeighborReturn::Internal(
                    i,
                    i_face,
                    self.mesh.faces.owner[i_face].1 + n_distance,
                ),
                NeighborType::Boundary(id) => {
                    let BoundaryCondition(name, values) = &self.equation_containers
                        [self.mesh.cells.zone[i_cell]]
                        .boundary_conditions[id];
                    NeighborReturn::Boundary(name, i_face, self.mesh.faces.owner[i_face].1, values)
                }
                _ => NeighborReturn::Nil,
            })
            .collect()
    }

    /// Retruns the "neighbor" side of a face.
    /// A neighbor can be :
    /// - an internal cell (index, distance between centers)
    /// - a boundary condition (name, some values)
    /// - Nil (no neighbor, should only be encountered in 2D or 1D meshes)
    pub fn neighbor(&self, i_face: usize) -> NeighborReturn {
        match self.mesh.faces.neighbor[i_face] {
            NeighborType::Internal(_, n_distance) => {
                let (i, o_distance) = self.mesh.faces.owner[i_face];
                NeighborReturn::Internal(i, i_face, o_distance + n_distance)
            }
            NeighborType::Boundary(id) => {
                let BoundaryCondition(name, values) =
                    &self.equation_containers[self.get_face_zone(i_face)].boundary_conditions[id];
                NeighborReturn::Boundary(name, i_face, self.mesh.faces.owner[i_face].1, values)
            }
            _ => NeighborReturn::Nil,
        }
    }

    /// Exports the mesh to the VTK format.
    /// vtkio crate support for polyhedra is not yet implemented.
    /// This function will only work for hexahedra.
    pub fn export_to_vtk(&self, filename: &str) {
        // collect the lists of points that define the cells
        // TODO : implement IntoIter for CellData and Properties to simplify this
        let mut cell_indices: Vec<Vec<u64>> = Vec::with_capacity(self.n_cells());
        for i in 0..self.n_cells() {
            cell_indices.push(
                self.mesh
                    .cells
                    .get_cell(i)
                    .hexahedron_vertices(&self.mesh)
                    .iter()
                    .map(|&x| x as u64)
                    .collect(),
            );
        }
        let offsets: Vec<u64> = cell_indices
            .iter()
            .scan(0, |acc, x| {
                *acc += x.len() as u64;
                Some(*acc)
            })
            .collect::<Vec<u64>>();
        // collect the data of the cells
        let mut data_labels = Vec::new();
        let mut data_values = Vec::new();
        for (label, values) in self.cell_properties.export() {
            data_labels.push(label);
            data_values.push(values);
        }

        // use the vtkio crate to write a VTK file
        let vtk_file = vtkio::Vtk {
            version: vtkio::model::Version { major: 4, minor: 2 },
            title: "Mesh".to_string(),
            byte_order: vtkio::model::ByteOrder::BigEndian,
            file_path: None,
            data: vtkio::model::DataSet::inline(vtkio::model::UnstructuredGridPiece {
                points: vtkio::IOBuffer::F64(
                    // vec that contains all point coordinates
                    self.mesh
                        .vertices
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
                                data: vtkio::IOBuffer::F64(values.to_vec()),
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

    pub fn export_csv(&self, filename: &str) -> Result<(), Box<dyn Error>> {
        // create the file
        File::create(filename)?;
        let mut wtr = csv::Writer::from_path(filename)?;

        // When writing records without Serde, the header record is written just
        // like any other record.
        let mut header = vec!["x", "y", "z"];
        let headerdata = self.cell_properties.export();
        for (label, _) in headerdata.iter() {
            header.push(label);
        }
        wtr.write_record(&header)?;
        for i in 0..self.n_cells() {
            let mut record = vec![
                self.mesh.cells.r[i].x.to_string(),
                self.mesh.cells.r[i].y.to_string(),
                self.mesh.cells.r[i].z.to_string(),
            ];
            for (_, value) in headerdata.iter() {
                record.push(value[i].to_string());
            }
            wtr.write_record(&record)?;
        }
        wtr.flush()?;
        Ok(())
    }
}

impl<T: HasProperties> Default for Problem<T> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use crate::computational_cell::HasProperties;
    use nalgebra::RowDVector;

    #[allow(non_snake_case)]
    #[derive(Clone)]
    pub struct Properties {
        // physical properties
        pub T: RowDVector<f64>,   // K
        pub rho: RowDVector<f64>, // kg/m3
        pub cp: RowDVector<f64>,  // J/k/kg
        // solid
        pub lambda: RowDVector<f64>, // W/m/k
    }

    impl HasProperties for Properties {
        fn new(_size: usize) -> Self {
            Properties {
                T: RowDVector::from(Vec::new()),
                rho: RowDVector::from(Vec::new()),
                cp: RowDVector::from(Vec::new()),
                lambda: RowDVector::from(Vec::new()),
            }
        }
        fn export(&self) -> Vec<(String, Vec<f64>)> {
            vec![("T".to_string(), self.T.data.as_slice().to_vec())]
        }
    }

    #[test]
    fn export_single_cell_mesh() {
        let mut problem: crate::problem::Problem<Properties> = crate::problem::Problem::new();
        problem.mesh.add_hex_cell(
            nalgebra::Point3::new(0., 0., 0.),
            nalgebra::Point3::new(1., 1., 1.),
            0,
        );
        problem.cell_properties = Properties {
            T: RowDVector::from_vec(vec![300.]),
            rho: RowDVector::from_vec(vec![1.]),
            cp: RowDVector::from_vec(vec![1000.]),
            lambda: RowDVector::from_vec(vec![1.]),
        };

        problem.export_to_vtk("single_cell_mesh.vtk");
    }

    #[test]
    fn export_double_cell_mesh() {
        let mut problem: crate::problem::Problem<Properties> = crate::problem::Problem::new();
        problem.mesh.add_hex_cell(
            nalgebra::Point3::new(0., 0., 0.),
            nalgebra::Point3::new(1., 1., 1.),
            0,
        );
        problem.mesh.add_hex_cell(
            nalgebra::Point3::new(1., 0., 0.),
            nalgebra::Point3::new(2., 1., 1.),
            0,
        );
        problem.cell_properties = Properties {
            T: RowDVector::from_vec(vec![300., 400.]),
            rho: RowDVector::from_vec(vec![1., 1.]),
            cp: RowDVector::from_vec(vec![1000., 1000.]),
            lambda: RowDVector::from_vec(vec![1., 1.]),
        };
        problem.export_to_vtk("double_cell_mesh.vtk");
    }
}
