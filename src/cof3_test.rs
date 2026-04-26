//! CoF₃ 空间群测试。
//!
//! CoF₃ 是六方 R-3c (#167) 结构，含 18 F + 6 Cr = 24 原子。
//! 本模块手动硬编码 POSCAR 数据，直接构造 Cell 并判断空间群。

#[cfg(test)]
mod test_cof3 {

    use crate::cell::{Cell, TensorRank};
    use crate::mathfunc::Mat3;
    use crate::pointgroup::ptg_get_pointgroup;

    /// 手动硬编码 CoF₃ 的 POSCAR 数据，直接构造 Cell 判断空间群。
    /// 预期：R-3c (#167), 点群 -3m (D3d), 36 个对称操作。
    #[test]
    fn test_cof3_hardcoded() {
        // --- POSCAR 手动写入 ---
        // 六方晶格: a=4.979, c=13.181, γ=120°
        // POSCAR 原始行向量 → 转置为 lattice[cart][vec]（列=晶格矢量）
        let lattice: Mat3 = [
            [4.9790000916,   -2.4895000458,    0.0           ], // x 分量
            [0.0,             4.3119405647,    0.0           ], // y 分量
            [0.0,             0.0,            13.1809997559  ], // z 分量
        ];

        // 18 F 原子 + 6 Cr 原子，Direct 坐标
        let positions = [
            // F (18 atoms)
            [0.617923975, 0.0,         0.25       ],
            [0.0,         0.617923975, 0.25       ],
            [0.382076025, 0.382076025, 0.25       ],
            [0.382076025, 0.0,         0.75       ],
            [0.0,         0.382076025, 0.75       ],
            [0.617923975, 0.617923975, 0.75       ],
            [0.284590632, 0.333333343, 0.583333313],
            [0.666666687, 0.951257288, 0.583333313],
            [0.048742693, 0.715409338, 0.583333313],
            [0.048742693, 0.333333343, 0.083333336],
            [0.666666687, 0.715409338, 0.083333336],
            [0.284590632, 0.951257288, 0.083333336],
            [0.951257288, 0.666666687, 0.916666687],
            [0.333333343, 0.284590632, 0.916666687],
            [0.715409338, 0.048742693, 0.916666687],
            [0.715409338, 0.666666687, 0.416666657],
            [0.333333343, 0.048742693, 0.416666657],
            [0.951257288, 0.284590632, 0.416666657],
            // Cr (6 atoms)
            [0.0,         0.0,         0.0        ],
            [0.333333343, 0.666666687, 0.166666672],
            [0.666666687, 0.333333343, 0.333333343],
            [0.0,         0.0,         0.5        ],
            [0.333333343, 0.666666687, 0.666666687],
            [0.666666687, 0.333333343, 0.833333313],
        ];

        // F=9, Cr=24
        let mut types = Vec::with_capacity(24);
        for _ in 0..18 { types.push(9_i32); }  // F
        for _ in 0..6  { types.push(24_i32); } // Cr

        eprintln!("=== CoF₃ 手动构造 Cell ===");
        eprintln!("n_atoms = {}", positions.len());

        // 验证晶格: 六方 a=b, γ=120°
        let a_col = [lattice[0][0], lattice[1][0], lattice[2][0]];
        let b_col = [lattice[0][1], lattice[1][1], lattice[2][1]];
        let len_a = (a_col[0]*a_col[0] + a_col[1]*a_col[1] + a_col[2]*a_col[2]).sqrt();
        let len_b = (b_col[0]*b_col[0] + b_col[1]*b_col[1] + b_col[2]*b_col[2]).sqrt();
        let dot = a_col[0]*b_col[0] + a_col[1]*b_col[1] + a_col[2]*b_col[2];
        let gamma = (dot / (len_a * len_b)).acos() * 180.0 / std::f64::consts::PI;
        eprintln!("|a| = {:.6}, |b| = {:.6}, γ = {:.4}°", len_a, len_b, gamma);
        assert!((len_a - len_b).abs() < 1e-6, "a must equal b for hexagonal");
        assert!((gamma - 120.0).abs() < 0.1, "γ must be 120° for hexagonal");

        // --- 构造 Cell ---
        let mut cell = Cell::new(24, TensorRank::NoSpin);
        cell.set_cell(&lattice, &positions, &types);
        cell.aperiodic_axis = None;

        // --- 获取原胞 ---
        let primitive = crate::primitive::prm_get_primitive(&cell, 1e-5, -1.0)
            .expect("primitive failed");
        eprintln!("primitive: n_atoms={}", primitive.size);

        // --- 获取对称操作 ---
        if let Some(ref prim_cell) = primitive.cell {
            let sym = crate::symmetry::sym_get_operation(prim_cell, 1e-5, -1.0)
                .expect("sym_get_operation failed");
            eprintln!("symmetry ops on primitive cell: {}", sym.size);

            // 点群分类
            let pointsym = crate::pointgroup::ptg_get_pointsymmetry(&sym.rot);
            eprintln!("unique point symmetry rotations: {}", pointsym.size);

            // 获取点群编号（通过内部函数不可用，直接用已知映射）
            // 打印所有独特旋转用于验证
            for (i, r) in pointsym.rot.iter().enumerate() {
                let det = crate::mathfunc::mat_get_determinant_i3(r);
                let trace = crate::mathfunc::mat_get_trace_i3(r);
                eprintln!("  rot[{}]: det={}, trace={}", i, det, trace);
            }
        }

        // --- 展示原胞结构 ---
        if let Some(ref prim_cell) = primitive.cell {
            let p_lat = &prim_cell.lattice;
            eprintln!("\n=== 约化后原胞结构 ===");
            eprintln!("n_atoms = {}", prim_cell.size);
            eprintln!("lattice = [{:?},", p_lat[0]);
            eprintln!("           {:?},", p_lat[1]);
            eprintln!("           {:?}]", p_lat[2]);
            // 列矢量的模长
            for j in 0..3 {
                let len = (p_lat[0][j]*p_lat[0][j] + p_lat[1][j]*p_lat[1][j] + p_lat[2][j]*p_lat[2][j]).sqrt();
                eprintln!("  |v{}| = {:.6}", j, len);
            }
            // 列矢量之间角度
            for j in 0..3 {
                for k in (j+1)..3 {
                    let dot = p_lat[0][j]*p_lat[0][k] + p_lat[1][j]*p_lat[1][k] + p_lat[2][j]*p_lat[2][k];
                    let len_j = (p_lat[0][j]*p_lat[0][j] + p_lat[1][j]*p_lat[1][j] + p_lat[2][j]*p_lat[2][j]).sqrt();
                    let len_k = (p_lat[0][k]*p_lat[0][k] + p_lat[1][k]*p_lat[1][k] + p_lat[2][k]*p_lat[2][k]).sqrt();
                    let angle = (dot / (len_j * len_k)).acos() * 180.0 / std::f64::consts::PI;
                    eprintln!("  ∠(v{}, v{}) = {:.4}°", j, k, angle);
                }
            }
            for i in 0..prim_cell.size {
                eprintln!("  atom {:2}: pos=[{:.6}, {:.6}, {:.6}] type={}",
                    i+1, prim_cell.position[i][0], prim_cell.position[i][1],
                    prim_cell.position[i][2], prim_cell.types[i]);
            }
        }

        // --- 搜索空间群 ---
        let spg = crate::spacegroup::spa_search_spacegroup(&primitive, 0, 1e-5, -1.0)
            .expect("spacegroup search failed");

        eprintln!("\n=== 空间群结果 ===");
        eprintln!("  Number:            {}", spg.number);
        eprintln!("  International:     {}", spg.international_short.trim());
        eprintln!("  International full:{}", spg.international.trim());
        eprintln!("  Schoenflies:       {}", spg.schoenflies.trim());
        eprintln!("  Hall number:       {}", spg.hall_number);
        eprintln!("  Point group number:{}", spg.pointgroup_number);

        // 点群
        let pg = ptg_get_pointgroup(spg.pointgroup_number);
        eprintln!("\n=== 点群 ===");
        eprintln!("  Number:     {}", pg.number);
        eprintln!("  Symbol:     {}", pg.symbol.trim());
        eprintln!("  Schoenflies:{}", pg.schoenflies.trim());

        // 验证空间群
        assert_eq!(spg.number, 167, "CoF₃ should be R-3c (#167)");

        // 点群 -3m 即 D3d
        assert_eq!(spg.pointgroup_number, 20, "point group should be -3m");
    }

    /// CoF₃ 含磁性: Cr AFM 沿 [001]，交替 ±z 自旋。
    /// POSCAR 中的 Cr 磁矩: +z 和 -z 交替。
    /// 预期: 非磁 R-3c (#167) → BNS=167.103 (Type-3, UNI=1333)
    #[test]
    fn test_cof3_magnetic() {
        let symprec = 1e-5;

        // 晶格与 POSCAR 相同
        let lattice: Mat3 = [
            [4.9790000916,   -2.4895000458,    0.0           ],
            [0.0,             4.3119405647,    0.0           ],
            [0.0,             0.0,            13.1809997559  ],
        ];

        // 18 F + 6 Cr (非磁胞)
        let positions = [
            [0.617923975, 0.0,         0.25       ],
            [0.0,         0.617923975, 0.25       ],
            [0.382076025, 0.382076025, 0.25       ],
            [0.382076025, 0.0,         0.75       ],
            [0.0,         0.382076025, 0.75       ],
            [0.617923975, 0.617923975, 0.75       ],
            [0.284590632, 0.333333343, 0.583333313],
            [0.666666687, 0.951257288, 0.583333313],
            [0.048742693, 0.715409338, 0.583333313],
            [0.048742693, 0.333333343, 0.083333336],
            [0.666666687, 0.715409338, 0.083333336],
            [0.284590632, 0.951257288, 0.083333336],
            [0.951257288, 0.666666687, 0.916666687],
            [0.333333343, 0.284590632, 0.916666687],
            [0.715409338, 0.048742693, 0.916666687],
            [0.715409338, 0.666666687, 0.416666657],
            [0.333333343, 0.048742693, 0.416666657],
            [0.951257288, 0.284590632, 0.416666657],
            [0.0,         0.0,         0.0        ],
            [0.333333343, 0.666666687, 0.166666672],
            [0.666666687, 0.333333343, 0.333333343],
            [0.0,         0.0,         0.5        ],
            [0.333333343, 0.666666687, 0.666666687],
            [0.666666687, 0.333333343, 0.833333313],
        ];

        let mut types = Vec::with_capacity(24);
        for _ in 0..18 { types.push(9_i32); }
        for _ in 0..6  { types.push(24_i32); }

        // 磁矩: F=0, Cr=交替 ±z (对应 POSCAR 中的 0 0 ±1)
        // POSCAR 中的 Cr 磁矩: Cr0=(0,0,1), Cr1=(0,0,-1), Cr2=(0,0,1),
        //                     Cr3=(0,0,-1), Cr4=(0,0,1), Cr5=(0,0,-1)
        let mut moments = vec![[0.0; 3]; 24];
        moments[18] = [0.0, 0.0,  1.0];  // Cr at [0,0,0]
        moments[19] = [0.0, 0.0, -1.0];  // Cr at [1/3,2/3,1/6]
        moments[20] = [0.0, 0.0,  1.0];  // Cr at [2/3,1/3,1/3]
        moments[21] = [0.0, 0.0, -1.0];  // Cr at [0,0,1/2]
        moments[22] = [0.0, 0.0,  1.0];  // Cr at [1/3,2/3,2/3]
        moments[23] = [0.0, 0.0, -1.0];  // Cr at [2/3,1/3,5/6]

        eprintln!("=== CoF₃ AFM [001] 磁性测试 ===");

        let result = crate::spg_get_magnetic_dataset(
            &lattice, &positions, &types, Some(&moments), symprec,
        ).expect("spg_get_magnetic_dataset must succeed");

        eprintln!("{}", crate::spg_format_magnetic_symmetry(&result));

        assert_eq!(result.spacegroup_number, 167, "non-mag: R-3c");
        assert_eq!(result.uni_number, 1333, "UNI=1333");
        assert_eq!(result.bns_number.trim(), "167.103");
        assert_eq!(result.magnetic_type, crate::MagneticType::Ordinary);
        eprintln!("UNI={}, BNS={}, type={:?}", result.uni_number, result.bns_number.trim(), result.magnetic_type);
    }
}
