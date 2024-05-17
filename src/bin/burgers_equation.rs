use flex_fvm::{
    computational_cell::HasProperties, equation::*, generate_equation_structure, problem, solver,
};
use nalgebra::{dvector, DVector, Point3, Vector3};
use std::fmt;

fn main() {
    let initial_mesh = create_mesh();
    let output_args = solver::SolverOutputArgs {
        frequency: 0.1,
        filename: "burger_time_series".to_string(),
        print_csv: true,
        print_vtk: false,
    };
    let _solution = solver::solve(
        &initial_mesh,
        50.,
        0.0001,
        solver::euler_explicit::<Properties>,
        output_args,
    );
}

////////////////////////
/// Cell Propertires ///
////////////////////////

#[allow(non_snake_case)]
#[derive(Default, Clone, Copy)]
pub struct Properties {
    // physical properties
    pub u: Vector3<f64>, // K
    pub nu: f64,         // kg/m3
}

impl fmt::Debug for Properties {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.u)
    }
}

impl HasProperties for Properties {
    fn update_primitive(&mut self, volume: f64) {
        // here we can update derived variables if needed
    }

    fn export(&self) -> Vec<(String, f64)> {
        vec![
            ("ux".to_string(), self.u[0]),
            ("uy".to_string(), self.u[1]),
            ("uz".to_string(), self.u[2]),
        ]
    }
}

pub fn material(position: &Point3<f64>) -> Properties {
    Properties {
        u: Vector3::new(f64::exp(-position.x.powi(2) / 2.), 0., 0.),
        nu: 1.,
    }
}

///////////////////////
/// Mesh Generation ///
///////////////////////

#[allow(non_snake_case)]
pub fn create_mesh() -> problem::Problem<Properties> {
    // pack the equations
    let heat_diffusion_equations =
        generate_equation_structure!(get_u, update_u, burgers_equation, 3);
    // define the boundary conditions
    // left
    let bc1 = BoundaryCondition("periodic".to_string(), dvector![1., 0., 0.]);
    // right
    let bc2 = BoundaryCondition("periodic".to_string(), dvector![-1., 0., 0.]);

    let n_grid = (500, 1, 1);

    let mut block: problem::Problem<Properties> = problem::Problem::new();
    // bottom left
    let boundary_conditions = [Some(bc1), None, None, Some(bc2), None, None];
    block.add_zone(
        (-7., 0., 0.),
        (7., 1., 1.),
        n_grid,
        boundary_conditions,
        vec![heat_diffusion_equations],
        material,
    );
    block
}

/////////////////////////
/// Burgers' equation ///
/////////////////////////
#[allow(non_snake_case)]
pub fn get_u(p: &Properties) -> Vector3<f64> {
    p.u
}

#[allow(non_snake_case)]
pub fn update_u(p: &mut Properties, u: &Vector3<f64>) {
    p.u[0] = u[0];
    p.u[1] = u[1];
    p.u[2] = u[2];
}

fn upwind<'a>(left: &'a Properties, right: &'a Properties, normal: Vector3<f64>) -> &'a Properties {
    if left.u.dot(&normal) >= 0. {
        left
    } else {
        right
    }
}

pub fn burgers_equation(problem: &problem::Problem<Properties>, i_p: usize) -> Vector3<f64> {
    // du/dt + div(u^2)/2 = 0
    let p = &problem.mesh.cells[i_p];
    let mut advective_flux = Vector3::zeros();
    for neighbor in problem.neighbors(i_p) {
        match neighbor {
            // flux term for a normal neighbor
            problem::NeighborReturn::Internal(i_n, i_f, _) => {
                let n = &problem.mesh.cells[i_n];
                let f = &problem.mesh.faces[i_f];
                let upwind_properties = upwind(&p.properties, &n.properties, f.out_normal(i_p));
                let u_upwind = upwind_properties.u;
                advective_flux -= 0.5 * u_upwind.dot(&u_upwind) * f.area * f.out_normal(i_p);
            }
            // flux term for a boundary condition
            problem::NeighborReturn::Boundary(name, _, _, _) => match name {
                _ => {
                    panic!("Unknown boundary condition : {}", name);
                }
            },
            problem::NeighborReturn::Nil => {
                // nothing to do
            }
        }
    }
    advective_flux / p.volume
}
