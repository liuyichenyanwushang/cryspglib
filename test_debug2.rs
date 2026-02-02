use rspglib::cell::Cell;
use rspglib::mathfunc::mat_copy_matrix_d3;
use rspglib::primitive::Primitive;
use rspglib::spacegroup::spa_search_spacegroup;

fn main() {
    // Construct a simple P1 cell
    let lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let position = [[0.0, 0.0, 0.0]];
    let types = [1];
    let mut cell = Cell::new(1, rspglib::cell::TensorRank::NoSpin);
    cell.set_cell(&lattice, &position, &types);

    let mut primitive = Primitive::new(1);
    primitive.cell = Some(cell);
    mat_copy_matrix_d3(&mut primitive.orig_lattice, &lattice);

    println!("Testing spa_search_spacegroup...");
    let spg = spa_search_spacegroup(&primitive, 0, 1e-5, -1.0);

    match spg {
        Some(s) => {
            println!("Success! Found spacegroup: {}", s.number);
        }
        None => {
            println!("Failed! spa_search_spacegroup returned None");
        }
    }
}
