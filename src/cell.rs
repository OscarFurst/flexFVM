// for now all cells are rectangular parallelepipeds
use nalgebra::{Vector3,Point3};
use std::{fmt, f64::EPSILON};

#[allow(non_snake_case)]
#[derive(Clone)]
pub struct Cell {
    // geometric properties
    pub r: Point3<f64>,
    pub volume: f64,
    pub vertices: Vec<Point3<f64>>,
    pub variant: CellType,
    // physical properties
    pub T: f64,
    pub rho: f64,
    pub cp: f64,
    pub lambda: f64,
    // this pair of values uniquely identifies the cell
    pub id: (usize,usize),
    // identifyers of neighbors
    pub faces: Vec<Face>,
}

impl Cell {
    pub fn create_rectangle(&mut self, origin:Point3<f64>, ex:Vector3<f64>, ey:Vector3<f64>, ez:Vector3<f64>) {
        // unit vectors ei, ej and ek in x, y and z
        self.r = origin + ex/2. + ey/2. + ez/2.;
        self.volume = ex[0] * ey[1] * ez[2];
        // dont touch variant and values and id
        let vertices = vec![
            origin,
            origin+ex,
            origin+ey,
            origin+ez,
            origin+ex+ey,
            origin+ex+ez,
            origin+ey+ez,
            origin+ex+ey+ez
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

    pub fn face_on_plane(&self, point:&Point3<f64>, normal:&Vector3<f64>) -> Option<&Face> {
        //println!("Testing cell {:?}", self.r);
        for face in &self.faces {
            //println!("   face {:?} ; res = {}", face.center, (face.center - point).dot(normal));
            if (face.center - point).dot(normal).abs() <= f64::EPSILON {
                return Some(face);
            }
        }
        return None;
    }

    pub fn boundary_from_cell(cell:&Cell, normal:&Vector3<f64>, face:&Face, bc:BoundaryCondition) -> Cell {
        // the boundary condition is not linked to the cell after the function, it needs to be done when linking the cells
        let boundary_cell = Cell{
            r: cell.r + normal,
            volume: cell.volume,
            vertices: Vec::from_iter(cell.vertices.iter().map(|v| v+normal)),
            variant: CellType::Boundary(bc),
            faces: vec![face.clone()],
            ..Default::default()
        };
        boundary_cell
    }
}

#[derive(Clone,Copy)]
pub struct Neighbor {
    pub index: (usize,usize),
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
    fn new(vertices:&Vec<Point3<f64>>) -> Self {
        // TODO: only implemented for rectangle
        let vects:Vector3<f64> = vertices.iter().map(|v| v.coords).sum(); // this is so ugly
        let center: Point3<f64> = Point3::new(0.,0.,0.) + vects / vertices.len() as f64; // I just want to average...
        let newface = Face{
            vertices: vertices.clone(),
            center: center,
            area: {
                let mut area = 0.;
                let v1 = vertices[1] - vertices[0];
                for i in 2..vertices.len() {
                    if v1.dot(&(vertices[i] - vertices[0])).abs() < f64::EPSILON {
                        area = v1.cross(&(vertices[i] - vertices[0])).norm().abs()
                    }
                }area
            },
            neighbor: None,
        };
        newface
    }

    pub fn get_origin(&self) -> Point3<f64>{
        // find the point closest to origin
        // that does not look nice...
        let length_vector = self.vertices.iter().map(|v| v.coords.abs().sum());
        let index_of_min = length_vector
            .enumerate()
            .min_by(|(_,a), (_,b)| a.partial_cmp(b).unwrap())
            .map(|(index,_)| index);
        self.vertices[index_of_min.unwrap()]
    }
}

#[derive(Clone,Debug)]
pub enum CellType {
    Internal,
    Boundary(BoundaryCondition),
    None,
}

#[derive(Clone,Debug)]
pub enum BoundaryCondition {
    Dirichlet(f64),
    Neumann(f64),
}

impl fmt::Debug for Cell {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:.1}", self.T)
    }
}

impl Default for Cell {
    fn default() -> Self {
        Cell {
            r: Point3::new(0.,0.,0.),
            vertices: Vec::new(),
            volume: 0.,
            T: 0.0,
            rho: 1.0,
            cp: 0.1,
            lambda: 0.01,
            variant: CellType::Internal,
            id: (0,0),
            faces: Vec::new(),
        }
    }
}