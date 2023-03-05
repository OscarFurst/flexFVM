use crate::cell::Properties;
use nalgebra::Vector3;

pub fn air() -> Properties {
    Properties {
        T: 300.,
        rho: 1000.,
        cp: 4200.,
        D: 0.,
        u: Vector3::new(1., 0., 0.),
        ..Default::default()
    }
}

pub fn copper() -> Properties {
    Properties {
        T: 300.,
        rho: 9000.,
        cp: 0.001,
        lambda: 300.,
        ..Default::default()
    }
}
