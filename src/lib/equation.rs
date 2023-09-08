use crate::cell;
use crate::cell::HasProperties;
use debug_cell::Ref;
use nalgebra::DVector;

#[derive(Clone)]
pub struct BoundaryCondition<T: HasProperties>(
    pub fn(&mut T, &[cell::Face<T>], &DVector<f64>),
    pub DVector<f64>,
);

// An Equation is a combination of a
#[derive(Clone, Copy)]
pub struct Equation<T: HasProperties>(
    pub fn(&mut T, &DVector<f64>),
    pub fn(&Ref<cell::Cell<T>>) -> DVector<f64>,
);

#[macro_export]
macro_rules! generate_update_function {
    // this macro transforms a user defined function that operates on a "convenient" data type like f64 or Vector3<f64>
    // into a function that operates on a DVector<f64>, which is used by the solver
    ($user_update_function:ident, f64) => {
        |p: &mut Properties, v2: &DVector<f64>| -> () {
            let v1: f64 = v2[0];
            $user_update_function(p, v1);
        } as fn(&mut Properties, &DVector<f64>)
    };
    ($user_update_function:ident, $SomeType:ty) => {
        |p: &mut Properties, v2: &DVector<f64>| -> () {
            let v1: $SomeType = $SomeType::from_slice(v2.as_slice());
            $user_update_function(p, v1);
        } as fn(&mut Properties, &DVector<f64>);
    };
}

#[macro_export]
macro_rules! generate_equation_function {
    // this macro transforms a user defined function that operates on a "convenient" data type like f64 or Vector3<f64>
    // into a function that operates on a DVector<f64>, which is used by the solver
    ($user_update_function:ident, f64) => {
        |p: &Ref<cell::Cell<Properties>>| -> DVector<f64> {
            let v1: f64 = $user_update_function(p);
            dvector![v1]
        } as fn(&Ref<cell::Cell<Properties>>) -> DVector<f64>
    };
    ($user_update_function:ident, $SomeType:ty) => {
        |p: &Ref<cell::Cell<Properties>>| -> DVector<f64> {
            let v1: $SomeType = $user_update_function(p);
            DVector::from_vec(v1.as_slice().to_vec())
        } as fn(&Ref<cell::Cell<Properties>>) -> DVector<f64>;
    };
}
