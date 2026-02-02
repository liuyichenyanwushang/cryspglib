extern crate rspglib;

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
            // Let's test sym_get_operation directly
            test_sym_get_operation();
        }
    }
}

fn test_sym_get_operation() {
    use rspglib::cell::Cell;
    use rspglib::symmetry::sym_get_operation;

    println!("\nTesting sym_get_operation directly...");
    let lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let position = [[0.0, 0.0, 0.0]];
    let types = [1];
    let mut cell = Cell::new(1, rspglib::cell::TensorRank::NoSpin);
    cell.set_cell(&lattice, &position, &types);

    let sym = sym_get_operation(&cell, 1e-5, -1.0);
    match sym {
        Some(s) => {
            println!("sym_get_operation succeeded with {} operations", s.size);
        }
        None => {
            println!("sym_get_operation returned None");
            // Test get_lattice_symmetry
            test_get_lattice_symmetry();
        }
    }
}

fn test_get_lattice_symmetry() {
    use rspglib::cell::Cell;
    use rspglib::symmetry::get_lattice_symmetry;

    println!("\nTesting get_lattice_symmetry directly...");
    let lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let position = [[0.0, 0.0, 0.0]];
    let types = [1];
    let mut cell = Cell::new(1, rspglib::cell::TensorRank::NoSpin);
    cell.set_cell(&lattice, &position, &types);

    let lat_sym = get_lattice_symmetry(&cell, 1e-5, -1.0);
    println!("get_lattice_symmetry returned size: {}", lat_sym.size);

    if lat_sym.size == 0 {
        println!("get_lattice_symmetry returned size 0!");
        // Lattice symmetry detection failed
    }
}

