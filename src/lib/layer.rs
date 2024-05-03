// x,y,z coordinate system, x being left-to-right, y being back-to-front and z being bottom-to-top
use crate::computational_cell::CellType::*;
use crate::computational_cell::*;
use crate::equation::{BoundaryCondition, Equation};
use indicatif::ProgressBar;
use nalgebra::{DVector, Point3, Vector3};
use std::{cell::RefCell, fmt, rc::Rc};

#[derive(Clone)]
pub struct Layer<T: HasProperties> {
    pub cells: Vec<Rc<RefCell<ComputationalCell<T>>>>,
    pub boundaries: [(Point3<f64>, Vector3<f64>, Vec<BoundaryCondition<T>>); 6],
    equations: Vec<Equation<T>>,
}

#[allow(dead_code)]
impl<T: HasProperties> Layer<T> {
    pub fn n_eq(&self) -> usize {
        // returns the total number of equations that need to be solved for the layer
        self.cells
            .iter()
            .filter(|c| matches!(c.borrow().variant, CellType::Internal))
            .count()
            * self.equations.len()
    }

    #[allow(non_snake_case)]
    pub fn generate(
        r1: (f64, f64, f64),
        r2: (f64, f64, f64),
        (nx, ny, nz): (usize, usize, usize),
        bc: [Vec<BoundaryCondition<T>>; 6],
        equations: Vec<Equation<T>>,
        internal_props: T,
        boundary_props: T,
    ) -> Layer<T> {
        // generates a regular grid of rectangular cells
        // TODO: I could not think of an elegant way to pass boundary conditions
        // in the current state, the order in which they are passed are very important!
        let origin = Point3::new(r1.0, r1.1, r1.2);
        let ending = Point3::new(r2.0, r2.1, r2.2);
        let ex = (ending - origin)[0] / nx as f64 * Vector3::new(1., 0., 0.);
        let ey = (ending - origin)[1] / ny as f64 * Vector3::new(0., 1., 0.);
        let ez = (ending - origin)[2] / nz as f64 * Vector3::new(0., 0., 1.);
        let boundaries = [
            (origin, -ex, bc[0].clone()),
            (origin, -ey, bc[1].clone()),
            (origin, -ez, bc[2].clone()),
            (ending, ex, bc[3].clone()),
            (ending, ey, bc[4].clone()),
            (ending, ez, bc[5].clone()),
        ];
        let mut mesh = Layer {
            cells: Vec::new(),
            boundaries: boundaries.clone(),
            equations,
        };
        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    let mut cell = ComputationalCell {
                        variant: CellType::Internal,
                        properties: internal_props,
                        ..Default::default()
                    };
                    cell.create_rectangle(
                        origin + ex * i as f64 + ey * j as f64 + ez * k as f64,
                        ex,
                        ey,
                        ez,
                    );
                    // check if at boundary. if so, create boundary condition
                    for (point, normal, bc) in &boundaries {
                        if let Some(face) = cell.face_on_plane(point, normal) {
                            // search for corresponding boundary condition
                            let mut boundary_cell = ComputationalCell::boundary_from_cell(
                                &cell,
                                normal,
                                face,
                                bc.clone(),
                            );
                            boundary_cell.properties = boundary_props;
                            mesh.cells.push(Rc::new(RefCell::new(boundary_cell)));
                        }
                    }
                    // neighbors can be set from a specialized routine that uses the vertices of the cells. it is needed anyway to merge multiple layers
                    mesh.cells.push(Rc::new(RefCell::new(cell)));
                }
            }
        }
        mesh
    }

    pub fn link(self) -> Self {
        // compare cells two by two
        let n_cells = self.cells.len();
        println!("🔗 Linking cells ...");
        let pb = ProgressBar::new(n_cells as u64);
        for i in 0..n_cells {
            pb.inc(1);
            let cell = &self.cells[i];
            for j in i + 1..n_cells {
                let other_cell = &self.cells[j];
                // exclude the possibility of them being the same cell
                let dr = other_cell.borrow().r - cell.borrow().r;
                if dr.norm() < f64::EPSILON {
                    continue;
                }
                // compare faces
                for face in cell
                    .borrow_mut()
                    .faces
                    .iter_mut()
                    .filter(|f| matches!(f.neighbour, None))
                {
                    for other_face in other_cell.borrow_mut().faces.iter_mut() {
                        // pair up faces that matched
                        if face.center.coords.relative_eq(
                            &other_face.center.coords,
                            f64::EPSILON,
                            f64::EPSILON,
                        ) {
                            face.neighbour = Some(Rc::clone(other_cell));
                            face.n = dr / dr.norm();
                            face.distance = dr.norm();

                            other_face.neighbour = Some(Rc::clone(cell));
                            other_face.n = -dr / dr.norm();
                            other_face.distance = dr.norm();
                        }
                    }
                }
            }
        }
        pb.finish();
        self
    }

    pub fn get_state(&self) -> Vec<DVector<f64>> {
        // returns the current state of the layer
        let mut x: Vec<DVector<f64>> = Vec::with_capacity(self.n_eq());
        for cell in self
            .cells
            .iter()
            .map(|rc| rc.borrow())
            .filter(|c| matches!(c.variant, CellType::Internal))
        {
            for eq in self.equations.iter() {
                x.push((eq.get)(&cell.properties));
            }
        }
        x
    }

    pub fn equations(&self) -> Vec<DVector<f64>> {
        let mut y = Vec::with_capacity(self.n_eq());
        for p in self
            .cells
            .iter()
            .map(|c| c.borrow())
            .filter(|c| matches!(c.variant, CellType::Internal))
        {
            for eq in self.equations.iter() {
                y.push((eq.rhs)(&p));
            }
        }
        y
    }

    pub fn update(&self, x: &[DVector<f64>]) {
        let mut i = 0;
        for mut cell in self.cells.iter().map(|c| c.borrow_mut()) {
            if matches!(cell.variant, CellType::Internal) {
                for eq in self.equations.iter() {
                    (eq.lhs)(&mut cell.properties, &x[i]);
                    i += 1;
                }
            }
            cell.properties.update();
        }
    }

    pub fn boundary_conditions(&self) {
        for cell in self.cells.iter() {
            let cellref = &mut cell.borrow_mut();
            if let Boundary(bc_vector) = &cellref.variant {
                let mut newproperties = cellref.properties;
                for BoundaryCondition(function, value) in bc_vector {
                    function(&mut newproperties, &cellref.faces, value);
                }
                cellref.properties = newproperties;
            }
        }
    }
}

impl<T: HasProperties> fmt::Debug for Layer<T> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "[ ")?;
        for cell in self
            .cells
            .iter()
            .map(|rc| rc.borrow())
            .filter(|c| matches!(c.variant, CellType::Internal))
        {
            write!(f, "{:?} ", cell.properties)?;
        }
        write!(f, "]")
    }
}

// #[cfg(test)]
// mod tests {
//     use crate::{equation::*, layer::Layer, materials};
//     #[test]
//     fn test_generation() {
//         let bc1 = BoundaryCondition(heat_dirichlet, 300.);
//         let bc2 = BoundaryCondition(heat_neumann, 0.);
//         let bc3 = BoundaryCondition(heat_dirichlet, 0.);
//         let bc: [Vec<BoundaryCondition>; 6] = [
//             vec![bc1],
//             vec![bc2],
//             vec![bc2],
//             vec![bc3],
//             vec![bc2],
//             vec![bc2],
//         ];
//         let (nx, ny, nz) = (2, 2, 2);
//         let layer = Layer::generate(
//             (0., 0., 0.),
//             (1., 1., 1.),
//             (nx, ny, nz),
//             bc,
//             vec![Equation(update_T, heat_diffusion)],
//             materials::copper(),
//             materials::copper(),
//         );
//         for cell in layer.cells.iter().map(|rc| rc.borrow()) {
//             println!("coord : {:?} ; type : {:?}", cell.r, cell.variant);
//         }
//         assert_eq!(layer.cells.len(), 32);
//     }

//     #[test]
//     fn test_large_mesh() {
//         let bc1 = BoundaryCondition(heat_dirichlet, 300.);
//         let bc2 = BoundaryCondition(heat_neumann, 0.);
//         let bc3 = BoundaryCondition(heat_dirichlet, 0.);
//         let bc: [Vec<BoundaryCondition>; 6] = [
//             vec![bc1],
//             vec![bc2],
//             vec![bc2],
//             vec![bc3],
//             vec![bc2],
//             vec![bc2],
//         ];
//         let (nx, ny, nz) = (9, 1, 1);
//         let layer = Layer::generate(
//             (0., 0., 0.),
//             (5., 1., 1.),
//             (nx, ny, nz),
//             bc,
//             vec![Equation(update_T, heat_diffusion)],
//             materials::copper(),
//             materials::copper(),
//         )
//         .link();
//         for cell in layer.cells.iter().map(|rc| rc.borrow()) {
//             println!("coord : {:?} ; type : {:?}", cell.r, cell.variant);
//         }
//         println!("Cell volume : {}", layer.cells[0].borrow().volume);
//         assert_eq!(
//             layer.cells.len(),
//             (nx + 2) * (ny + 2) * (nz + 2) - 4 * nx - 4 * ny - 4 * nz - 8
//         );
//         for cell in layer.cells.iter() {
//             for face in &cell.borrow().faces {
//                 match face.neighbour {
//                     None => {
//                         println!(
//                             "unbound face, coord : {:?} ; type : {:?}, face center : {:?}",
//                             cell.borrow().r,
//                             cell.borrow().variant,
//                             face.center
//                         )
//                     }
//                     _ => {}
//                 }
//             }
//         }
//     }
// }
