use crate::computational_cell::HasProperties;
use nalgebra::DVector;

/// A BoundaryCondition is a combination of :
/// - a name
/// - a DVector
/// Basically, it only tells the equation which BC to apply, optionally including some data.
#[derive(Clone)]
pub struct BoundaryCondition(pub String, pub DVector<f64>);

// An Equation is a combination of :
// - a function that "gets" the derivated property
// - a function that "sets" the derivated property
// - a function that computes the derivative of the property
// - the dimension of the property, e.g. 1 for a scalar, 3 for a velocity vector, ...
// - a list (Vec) of possible boundary conditions
#[derive(Clone)]
pub struct Equation<T: HasProperties> {
    pub get: fn(&T) -> DVector<f64>,
    pub lhs: fn(&mut T, &DVector<f64>),
    pub rhs: fn(&crate::mesh::Mesh<T>, usize) -> DVector<f64>,
    pub dimension: usize,
}

#[macro_export]
macro_rules! generate_equation_structure {
    ($get_function:ident, $update_function:ident, $equation_function:ident, 1) => {
        Equation {
            get: flex_fvm::generate_get_function!($get_function, 1),
            lhs: flex_fvm::generate_update_function!($update_function, 1),
            rhs: flex_fvm::generate_equation_function!($equation_function, 1),
            dimension: 1,
        }
    };
    ($get_function:ident, $update_function:ident, $equation_function:ident, $SomeType:expr) => {
        Equation {
            get: flex_fvm::generate_get_function!($get_function, $SomeType),
            lhs: flex_fvm::generate_update_function!($update_function, $SomeType),
            rhs: flex_fvm::generate_equation_function!($equation_function, $SomeType),
            dimension: $SomeType,
        }
    };
}

#[macro_export]
macro_rules! generate_get_function {
    // this macro transforms a user defined function that operates on a "convenient" data type like f64 or Vector3<f64>
    // into a function that operates on a DVector<f64>, which is used by the solver
    ($user_get_function:ident, 1) => {
        |p: &Properties| -> DVector<f64> { dvector![$user_get_function(p)] }
            as fn(&Properties) -> DVector<f64>
    };
    ($user_get_function:ident, $SomeType:expr) => {
        |p: &Properties| -> DVector<f64> {
            nalgebra::DVector::from_vec($user_get_function(p).as_slice().to_vec()) // TODO : this might be slow as heck
        } as fn(&Properties) -> DVector<f64>
    };
}

#[macro_export]
macro_rules! generate_update_function {
    // this macro transforms a user defined function that operates on a "convenient" data type like f64 or Vector3<f64>
    // into a function that operates on a DVector<f64>, which is used by the solver
    ($user_update_function:ident, 1) => {
        |p: &mut Properties, v2: &DVector<f64>| -> () {
            let v1: f64 = v2[0];
            $user_update_function(p, v1);
        } as fn(&mut Properties, &DVector<f64>)
    };
    ($user_update_function:ident, $SomeType:expr) => {
        |p: &mut Properties, v2: &DVector<f64>| -> () {
            let v1: nalgebra::SVector<f64, $SomeType> = nalgebra::SVector::from_vec(v2.as_slice().to_vec()); // TODO : this might be slow as heck
            $user_update_function(p, &v1);
        } as fn(&mut Properties, &DVector<f64>);
    };
}

#[macro_export]
macro_rules! generate_equation_function {
    // this macro transforms a user defined function that operates on a "convenient" data type like f64 or Vector3<f64>
    // into a function that operates on a DVector<f64>, which is used by the solver
    ($user_equation_function:ident, 1) => {
        |m: &crate::mesh::Mesh<Properties>, p: usize| -> DVector<f64> {
            let v1: f64 = $user_equation_function(m, p);
            dvector![v1]
        } as fn(&crate::mesh::Mesh<Properties>, usize) -> DVector<f64>
    };
    ($user_equation_function:ident, $SomeType:expr) => {
        |m: &crate::mesh::Mesh<Properties>, p: usize| -> DVector<f64> {
            let v1: nalgebra::SVector<f64, $SomeType> = $user_equation_function(m, p);
            DVector::from_vec(v1.as_slice().to_vec())
        } as fn(&crate::mesh::Mesh<Properties>, usize) -> DVector<f64>;
    };
}
