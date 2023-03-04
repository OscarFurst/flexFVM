// for now all cells are rectangular parallelepipeds
use debug_cell::RefCell;
use nalgebra::{Point3, Vector3};
use std::fmt;
use std::rc::Rc;

use crate::equation;

///////////////////////////
/// CELL and PROPERTIES ///
///////////////////////////

#[allow(non_snake_case)]
pub struct Cell {
    // geometric properties
    pub r: Point3<f64>, // m
    pub volume: f64,    // m3
    pub vertices: Vec<Point3<f64>>,
    pub variant: CellType,
    pub properties: Properties,
    // faces of the volume
    pub faces: Vec<Face>,
}

#[allow(non_snake_case)]
#[derive(Clone, Copy)]
pub struct Properties {
    // physical properties
    pub T: f64,      // K
    pub rho: f64,    // kg/m3
    pub cp: f64,     // J/k/kg
    pub lambda: f64, // W/m/k
}

impl Cell {
    pub fn create_rectangle(
        &mut self,
        origin: Point3<f64>,
        ex: Vector3<f64>,
        ey: Vector3<f64>,
        ez: Vector3<f64>,
    ) {
        // creates a rectangular parallelepiped
        self.r = origin + ex / 2. + ey / 2. + ez / 2.;
        self.volume = ex[0] * ey[1] * ez[2];
        // values are not initialized by this function!
        let vertices = vec![
            origin,
            origin + ex,
            origin + ey,
            origin + ez,
            origin + ex + ey,
            origin + ex + ez,
            origin + ey + ez,
            origin + ex + ey + ez,
        ];
        self.faces = vec![
            Face::new(&vec![vertices[0], vertices[1], vertices[2], vertices[4]]), // south
            Face::new(&vec![vertices[3], vertices[5], vertices[6], vertices[7]]), // north
            Face::new(&vec![vertices[0], vertices[1], vertices[3], vertices[5]]), // close
            Face::new(&vec![vertices[2], vertices[4], vertices[6], vertices[7]]), // far
            Face::new(&vec![vertices[1], vertices[4], vertices[5], vertices[7]]), // east
            Face::new(&vec![vertices[0], vertices[2], vertices[3], vertices[6]]), // west
        ];
        self.vertices = vertices;
    }

    pub fn face_on_plane(&self, point: &Point3<f64>, normal: &Vector3<f64>) -> Option<&Face> {
        // finds out if "self" has a face on a boudnary plane
        self.faces
            .iter()
            .find(|&face| (face.center - point).dot(normal).abs() <= f64::EPSILON)
    }

    pub fn boundary_from_cell(
        cell: &Cell,
        normal: &Vector3<f64>,
        face: &Face,
        bc: Vec<equation::BoundaryCondition>,
    ) -> Cell {
        // creates a new "boundary cell" from the provided internal cell
        // the boundary condition is not linked to the cell after the function, it needs to be done when linking the cells
        let boundary_cell = Cell {
            r: cell.r + normal,
            volume: cell.volume,
            vertices: Vec::from_iter(cell.vertices.iter().map(|v| v + normal)),
            variant: CellType::Boundary(bc),
            faces: vec![face.clone()],
            ..Default::default()
        };
        boundary_cell
    }
}

impl fmt::Debug for Cell {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:.1}", self.properties.T)
    }
}

impl Default for Cell {
    fn default() -> Self {
        Cell {
            r: Point3::new(0., 0., 0.),
            vertices: Vec::new(),
            volume: 0.,
            properties: Properties {
                ..Default::default()
            },
            variant: CellType::Internal,
            faces: Vec::new(),
        }
    }
}

impl Default for Properties {
    fn default() -> Self {
        Properties {
            T: 0.0,
            rho: 9000.0,
            cp: 0.001,
            lambda: 300.,
        }
    }
}

/////////////////////////
/// NEIGHBOR and FACE ///
/// /////////////////////

#[derive(Clone)]
pub struct Neighbor {
    pub cell: Rc<RefCell<Cell>>,
    pub dr: Vector3<f64>,
    pub distance: f64,
}

#[derive(Clone)]
pub struct Face {
    pub vertices: Vec<Point3<f64>>,
    pub center: Point3<f64>,
    pub area: f64,
    pub neighbor: Option<Neighbor>,
}

impl Face {
    fn new(vertices: &Vec<Point3<f64>>) -> Self {
        // TODO: only implemented for rectangle
        let vects: Vector3<f64> = vertices.iter().map(|v| v.coords).sum(); // this is so ugly
        let center: Point3<f64> = Point3::new(0., 0., 0.) + vects / vertices.len() as f64; // I just want to average...
        Face {
            vertices: vertices.clone(),
            center,
            area: {
                let mut area = 0.;
                let v1 = vertices[1] - vertices[0];
                for i in 2..vertices.len() {
                    if v1.dot(&(vertices[i] - vertices[0])).abs() < f64::EPSILON {
                        area = v1.cross(&(vertices[i] - vertices[0])).norm().abs()
                    }
                }
                area
            },
            neighbor: None,
        }
    }

    pub fn get_origin(&self) -> Point3<f64> {
        // find the point closest to origin
        // that does not look nice...
        let length_vector = self.vertices.iter().map(|v| v.coords.abs().sum());
        let index_of_min = length_vector
            .enumerate()
            .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .map(|(index, _)| index);
        self.vertices[index_of_min.unwrap()]
    }
}

//////////////////////////////////////
/// CELLTYPE and BOUNDARYCONDITION ///
//////////////////////////////////////
#[derive(Clone)]
pub enum CellType {
    Internal,
    Boundary(Vec<equation::BoundaryCondition>),
    Nil,
}

impl fmt::Debug for CellType {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "TODO")
    }
}
