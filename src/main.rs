use flex_fvm::examples;
use flex_fvm::solver;

fn main() {
    let initial_mesh = examples::copper_example();
    let solution = solver::solve(&initial_mesh, 0.001, 0.000001);
    println!("{:?}", solution)
}
