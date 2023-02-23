use crate::cell::Cell;
use debug_cell::Ref;

pub fn heat_diffusion(p: &Ref<Cell>) -> f64 {
    let mut flux = 0.;
    for f in p.faces.iter() {
        let n = f.neighbor.as_ref().unwrap().cell.borrow();
        flux -= p.lambda * (p.T - n.T) / f.neighbor.as_ref().unwrap().distance * f.area;
    }
    flux / p.volume / p.rho / p.cp
}
