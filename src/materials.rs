use crate::cell::Properties;

pub fn water() -> Properties {
    Properties {
        T: 300.,
        rho: 1000.,
        cp: 4200.,
        D: 0.,
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
