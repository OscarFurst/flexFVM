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
