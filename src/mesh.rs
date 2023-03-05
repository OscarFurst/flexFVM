use crate::cell::{self, Cell, Neighbor, Properties};
use crate::equation;
use crate::layer;
use debug_cell::Ref;
use nalgebra::DVector;
use std::fmt;
use std::rc::Rc;

#[derive(Clone)]
pub struct Mesh {
    layers: Vec<layer::Layer>,
}

impl Mesh {
    pub fn new() -> Self {
        Mesh { layers: Vec::new() }
    }

    pub fn vec_size(&self) -> usize {
        // returns the total number of equations that need to be solved in the domain
        let mut n = 0;
        for layer in self.layers.iter() {
            n += layer.n_eq();
        }
        n
    }

    pub fn add_layer(
        &mut self,
        r1: (f64, f64, f64),
        r2: (f64, f64, f64),
        (nx, ny, nz): (usize, usize, usize),
        bc: [Vec<equation::BoundaryCondition>; 6],
        equations: Vec<fn(&Ref<Cell>) -> f64>,
        internal_props: Properties,
        boundary_props: Properties,
    ) {
        self.layers.push(
            layer::Layer::generate(
                r1,
                r2,
                (nx, ny, nz),
                bc,
                equations,
                internal_props,
                boundary_props,
            )
            .link(),
        );
        self.join_layers();
    }

    fn join_layers(&self) {
        // searches 'Boundary' cells inside of the domain, eliminates them and links their corresponding internal cells together
        let mut to_pop = vec![];
        for layer in self.layers.iter() {
            for other_layer in self.layers.iter() {
                if std::ptr::eq(layer, other_layer) {
                    continue;
                }
                for cell in layer.cells.iter() {
                    if matches!(cell.borrow().variant, cell::CellType::Internal) {
                        continue;
                    }
                    for other_cell in other_layer.cells.iter() {
                        if matches!(other_cell.borrow().variant, cell::CellType::Internal) {
                            continue;
                        }
                        let r = cell.borrow().faces[0].center;
                        if (r - other_cell.borrow().faces[0].center).norm() < f64::EPSILON {
                            // we found a pair
                            // make copies of the cells they are pointing to
                            let cell_ref =
                                Rc::clone(&cell.borrow().faces[0].neighbor.as_ref().unwrap().cell);
                            let other_cell_ref = Rc::clone(
                                &other_cell.borrow().faces[0].neighbor.as_ref().unwrap().cell,
                            );
                            let other_cell_ref_givaway = Rc::clone(
                                &other_cell.borrow().faces[0].neighbor.as_ref().unwrap().cell,
                            );
                            let dr = other_cell_ref.borrow().r - cell_ref.borrow().r;

                            cell_ref
                                .borrow_mut()
                                .faces
                                .iter_mut()
                                .find(|f| (r - f.center).norm() < f64::EPSILON)
                                .unwrap() // found the right face
                                .neighbor = Some(Neighbor {
                                cell: other_cell_ref_givaway,
                                n: dr,
                                distance: dr.norm(),
                            });
                            other_cell_ref
                                .borrow_mut()
                                .faces
                                .iter_mut()
                                .find(|f| (r - f.center).norm() < f64::EPSILON)
                                .unwrap() // found the right face
                                .neighbor = Some(Neighbor {
                                cell: cell_ref,
                                n: -dr,
                                distance: dr.norm(),
                            });
                        }
                        to_pop.push(Rc::clone(cell));
                        to_pop.push(Rc::clone(other_cell));
                    }
                }
            }
        }
    }

    pub fn enforce_boundary_conditions(&mut self) {
        for layer in self.layers.iter() {
            layer.boundary_conditions();
        }
    }

    pub fn setup_equations(&self, _t: &f64, _x: &DVector<f64>) -> DVector<f64> {
        let mut y = vec![];
        for layer in self.layers.iter() {
            y.append(&mut layer.equations());
        }
        DVector::from_vec(y)
    }

    pub fn update(&mut self, y: &DVector<f64>) {
        // takes the solution and updates every cell with it
        let mut i = 0;
        for layer in self.layers.iter() {
            for mut cell in layer.cells.iter().map(|c| c.borrow_mut()) {
                if matches!(cell.variant, cell::CellType::Internal) {
                    cell.properties.T = y[i];
                    i += 1;
                }
            }
        }
        self.enforce_boundary_conditions();
    }
}

impl fmt::Debug for Mesh {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for layer in self.layers.iter() {
            write!(f, "{:?}", layer)?
        }
        write!(f, " ")
    }
}

impl Default for Mesh {
    fn default() -> Self {
        Mesh::new()
    }
}

#[cfg(test)]
mod tests {
    use crate::{equation, materials, mesh};
    #[test]
    fn test_generation() {
        let bc1 = equation::BoundaryCondition(equation::heat_dirichlet, 300.);
        let bc2 = equation::BoundaryCondition(equation::heat_neumann, 0.);
        let bc3 = equation::BoundaryCondition(equation::heat_dirichlet, 0.);
        let boundary_conditions: [Vec<equation::BoundaryCondition>; 6] = [
            vec![bc1],
            vec![bc2],
            vec![bc2],
            vec![bc3],
            vec![bc2],
            vec![bc2],
        ];
        let (nx, ny, nz) = (3, 1, 1);
        let mut block: mesh::Mesh;
        block = mesh::Mesh::new();
        block.add_layer(
            (0., 0., 0.),
            (1., 1., 1.),
            (nx, ny, nz),
            boundary_conditions,
            vec![equation::heat_diffusion],
            materials::copper(),
            materials::copper(),
        );
    }
}
