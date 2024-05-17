use crate::computational_cell::{HasProperties, NeighborType};
use crate::equation::{self, BoundaryCondition};
use crate::equation_container::EquationContainer;
use crate::mesh::Mesh;
use indicatif::ProgressBar;
use nalgebra::{DVector, Point3, Vector3};

pub enum NeighborReturn<'a> {
    Internal(usize, usize, f64),
    Boundary(&'a str, usize, f64, &'a DVector<f64>),
    Nil, // should only exist in 1D or 2D cases, where there is no neighbor in the missing dimension
}

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
        bc: [Option<BoundaryCondition>; 6],
        equations: Vec<equation::Equation<T>>,
        internal_props: fn(&Point3<f64>) -> T,
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
        let mut zone_cell_indices = Vec::new();
        for i in 0..nx {
            pb.inc(1);
            for j in 0..ny {
                for k in 0..nz {
                    let cell_origin = origin + ex * i as f64 + ey * j as f64 + ez * k as f64;
                    let cell_ending =
                        origin + ex * (i + 1) as f64 + ey * (j + 1) as f64 + ez * (k + 1) as f64;
                    let cell_center = cell_origin + (cell_ending - cell_origin) / 2.0;
                    let cell_index = self.mesh.add_hex_cell(
                        cell_origin,
                        cell_ending,
                        internal_props(&cell_center),
                        zone_id,
                    );
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
            let face_center = self
                .mesh
                .get_face_on_plane(point, normal, cell_index)
                .unwrap()
                .center;
            let mut neighb = NeighborType::Nil;
            match name as &str {
                "periodic" => {
                    // search for the other face with the same boundary condition
                    let mut found_twin = false;
                    if self.mesh.faces.len() != 0 {
                        let twin_direction = Vector3::new(values[0], values[1], values[2]);
                        for i in 0..self.mesh.faces.len() {
                            if let NeighborType::Boundary(face_bc_index) =
                                self.mesh.faces[i].neighbor
                            {
                                if self.equation_containers[self.mesh.faces[i].get_zone(self)]
                                    .boundary_conditions[face_bc_index]
                                    .0
                                    == "periodic"
                                {
                                    // we found another face with a periodic BC, now we need to check if it is the twin face
                                    let face = &mut self.mesh.faces[i];
                                    if (face.center - face_center).cross(&twin_direction).norm()
                                        < f64::EPSILON
                                    {
                                        found_twin = true;
                                        // in this case we have to set up two faces instead of jsut one !
                                        face.neighbor = NeighborType::Internal(
                                            cell_index,
                                            (face.center - self.mesh.cells[cell_index].r).norm(),
                                        );
                                        neighb = NeighborType::Internal(
                                            face.owner.0,
                                            (self.mesh.cells[face.owner.0].r - face.center).norm(),
                                        );
                                        break;
                                    }
                                }
                            }
                        }
                    }
                    if !found_twin {
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
            boundary_face.neighbor = neighb;
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
        for face in self.mesh.faces.iter() {
            if let NeighborType::Boundary(bc_index) = face.neighbor {
                if self.equation_containers[face.get_zone(self)].boundary_conditions[bc_index].0
                    == "periodic"
                {
                    return false;
                }
            }
        }
        true
    }

    pub fn sum_conservation_variables(&self) -> Vec<DVector<f64>> {
        let mut sums: Vec<DVector<f64>> = Vec::new();
        for layer in self.equation_containers.iter() {
            for eq in layer.equations.iter() {
                sums.push(
                    layer
                        .cell_indices
                        .iter()
                        .map(|&i| &self.mesh.cells[i])
                        .map(|cell| (eq.get)(&cell.properties))
                        .sum(),
                );
            }
        }
        sums
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
                    y.push((eq.rhs)(self, *i_cell))
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
                self.mesh.cells[*i_cell].update_primitive()
            }
        }
    }

    /// Returns all neighbors of the specified cell
    /// A neighbor can be :
    /// - an internal cell (index, distance between centers)
    /// - a boundary condition (name, some values)
    /// - Nil (no neighbor, should only be encountered in 2D or 1D meshes)
    pub fn neighbors<'a>(&self, i_cell: usize) -> Vec<NeighborReturn> {
        self.mesh.cells[i_cell]
            .faces
            .iter()
            .map(|i_face| match self.mesh.faces[*i_face].neighbor {
                NeighborType::Internal(i, n_distance) if i == i_cell => {
                    let (i, o_distance) = self.mesh.faces[*i_face].owner;
                    NeighborReturn::Internal(i, *i_face, o_distance + n_distance)
                }
                NeighborType::Internal(i, n_distance) => NeighborReturn::Internal(
                    i,
                    *i_face,
                    self.mesh.faces[*i_face].owner.1 + n_distance,
                ),
                NeighborType::Boundary(id) => {
                    let BoundaryCondition(name, values) = &self.equation_containers
                        [self.mesh.cells[i_cell].zone]
                        .boundary_conditions[id];
                    NeighborReturn::Boundary(
                        name,
                        *i_face,
                        self.mesh.faces[*i_face].owner.1,
                        values,
                    )
                }
                _ => NeighborReturn::Nil,
            })
            .collect()
    }
}

impl<T: HasProperties> Default for Problem<T> {
    fn default() -> Self {
        Self::new()
    }
}
