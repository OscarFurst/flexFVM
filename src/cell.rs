#[allow(dead_code)]
use nalgebra::Vector3;
use std::fmt;

#[allow(non_snake_case)]
#[derive(Clone)]
pub struct Cell {
    // geometric properties
    pub r: Vector3<f64>,
    pub r1: Vector3<f64>,
    pub r2: Vector3<f64>,
    pub dx: f64,
    pub dy: f64,
    pub dz: f64,
    pub dn: f64, pub ds: f64, pub de: f64, pub dw: f64, pub dc: f64, pub df: f64,
    pub variant: CellType,
    // physical properties
    pub T: f64,
    pub rho: f64,
    pub cp: f64,
    pub lambda: f64,
    // this pair of values uniquely identifies the cell
    pub id: usize,
    pub layer_id: usize,
    // identifyers of neighbors
    pub north: (usize,usize),
    pub south: (usize,usize),
    pub east: (usize,usize),
    pub west: (usize,usize),
    pub close: (usize,usize),
    pub far: (usize,usize),
}

#[derive(Clone,Copy)]
pub enum CellType {
    Internal,
    Boundary(BoundaryCondition),
    None,
}

#[derive(Clone,Copy)]
pub enum BoundaryCondition {
    Dirichlet(f64, Option<(usize,usize)>),
    Neumann(f64, Option<(usize,usize)>),
}

pub enum BoundaryDirection {
    North(BoundaryCondition),
    South(BoundaryCondition),
    East(BoundaryCondition),
    West(BoundaryCondition),
    Close(BoundaryCondition),
    Far(BoundaryCondition),
}

impl fmt::Debug for Cell {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:.2}", self.T)
    }
}

impl Default for Cell {
    fn default() -> Self {
        Cell {
            r: Vector3::new(0.,0.,0.),
            r1: Vector3::new(0.,0.,0.),
            r2: Vector3::new(0.,0.,0.),
            dx: 0.,
            dy: 0.,
            dz: 0.,
            dn: 0., ds: 0., de: 0., dw: 0., dc: 0., df: 0.,
            T: 0.0,
            rho: 1.0,
            cp: 0.1,
            lambda: 0.01,
            variant: CellType::Internal,
            id: 0,
            layer_id: 0,
            north: (0,0),
            south: (0,0),
            east: (0,0),
            west: (0,0),
            close: (0,0),
            far: (0,0),
        }
    }
}