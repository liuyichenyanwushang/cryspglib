extern crate rspglib;

use rspglib::cell::Cell;
use rspglib::symmetry::sym_get_operation;
use rspglib::delaunay::del_delaunay_reduce;

fn test_delaunay() {
    println!("=== 测试 Delaunay 约化 ===");

    let lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let mut min_lattice = [[0.0; 3]; 3];

    println!("输入晶格:");
    for row in &lattice {
        println!("  {:?}", row);
    }

    let success = del_delaunay_reduce(&mut min_lattice, &lattice, 1e-5);
    println!("del_delaunay_reduce 成功: {}", success);
    if success {
        println!("约化后晶格:");
        for row in &min_lattice {
            println!("  {:?}", row);
        }

        // 检查是否与原始相同
        let mut is_identical = true;
        for i in 0..3 {
            for j in 0..3 {
                if f64::abs(min_lattice[i][j] - lattice[i][j]) > 1e-10 {
                    is_identical = false;
                }
            }
        }
        println!("与原始晶格相同: {}", is_identical);
    }
}

fn test_get_lattice_symmetry() {
    println!("\n=== 测试 get_lattice_symmetry ===");

    let lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let position = [[0.0, 0.0, 0.0]];
    let types = [1];
    let mut cell = Cell::new(1, rspglib::cell::TensorRank::NoSpin);
    cell.set_cell(&lattice, &position, &types);

    println!("测试立方晶格 (单位矩阵)");
    println!("symprec = 1e-5, angle_tolerance = -1.0");

    let sym_opt = sym_get_operation(&cell, 1e-5, -1.0);
    if let Some(sym) = &sym_opt {
        println!("找到的对称操作数量: {}", sym.size);
    } else {
        println!("sym_get_operation 返回 None");
    }

    match sym_opt {
        Some(sym) => {
            println!("成功！应为48个操作（立方晶格）");
            if sym.size == 48 {
                println!("✅ 找到正确的48个对称操作");
            } else {
                println!("⚠ 找到 {} 个操作，期望48个", sym.size);
            }
        }
        None => {
            println!("警告: sym_get_operation 返回 None！");

            // 尝试更大的容差
            println!("\n尝试更大的容差: symprec = 1e-4");
            let sym2 = sym_get_operation(&cell, 1e-4, -1.0);
            match sym2 {
                Some(s) => println!("找到的对称操作数量: {}", s.size),
                None => println!("仍然返回 None"),
            }

            println!("\n尝试更大的容差: symprec = 1e-3");
            let sym3 = sym_get_operation(&cell, 1e-3, -1.0);
            match sym3 {
                Some(s) => println!("找到的对称操作数量: {}", s.size),
                None => println!("仍然返回 None"),
            }

            println!("\n尝试更大的容差: symprec = 1e-2");
            let sym4 = sym_get_operation(&cell, 1e-2, -1.0);
            match sym4 {
                Some(s) => println!("找到的对称操作数量: {}", s.size),
                None => println!("仍然返回 None"),
            }

            println!("\n尝试更大的容差: symprec = 1e-1");
            let sym5 = sym_get_operation(&cell, 1e-1, -1.0);
            match sym5 {
                Some(s) => println!("找到的对称操作数量: {}", s.size),
                None => println!("仍然返回 None"),
            }
        }
    }
}

fn main() {
    println!("调试晶格对称性检测");
    println!("====================\n");

    test_delaunay();
    test_get_lattice_symmetry();

    println!("\n测试完成");
}