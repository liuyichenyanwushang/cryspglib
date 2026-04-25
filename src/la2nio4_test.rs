//! La₂NiO₄ 空间群测试。
//!
//! 四方 P4₂/ncm (#138)，D₄ₕ 点群。晶胞含 8 La + 4 Ni + 16 O = 28 原子。

#[cfg(test)]
mod test_la2nio4 {

    use crate::cell::{Cell, TensorRank};
    use crate::mathfunc::Mat3;
    use crate::pointgroup::ptg_get_pointgroup;

    /// 手动硬编码 La₂NiO₄ POSCAR 数据，构造 Cell 并判断空间群。
    ///
    /// 该结构是四方 P4₂/ncm (#138)。28 个原子全部保留在原胞中
    /// （未检测到纯平移——POSCAR 中 La 原子排列打破了 I4/mmm 的体心平移）。
    ///
    /// 对称性: 16 个操作 (D₄ₕ 点群)。
    #[test]
    fn test_la2nio4_hardcoded() {
        // --- POSCAR 手动写入 ---
        // 四方晶格: a=b=5.4995, c=12.5052, α=β=γ=90°
        // POSCAR 行向量 → lattice[cart][vec]（列=晶格矢量）
        let lattice: Mat3 = [
            [5.4994997978, 0.0,            0.0           ], // x 分量
            [0.0,          5.4994997978,   0.0           ], // y 分量
            [0.0,          0.0,           12.5052003860  ], // z 分量
        ];

        // 8 La + 4 Ni + 16 O，POSCAR Direct 坐标
        let positions = [
            // La (8 atoms), type=57
            [0.992799997, 0.992799997, 0.363900006],
            [0.492799997, 0.007200000, 0.136099994],
            [0.007200000, 0.492799997, 0.136099994],
            [0.507200003, 0.507200003, 0.363900006],
            [0.007200000, 0.007200000, 0.636099994],
            [0.507200003, 0.992799997, 0.863900006],
            [0.992799997, 0.507200003, 0.863900006],
            [0.492799997, 0.492799997, 0.636099994],
            // Ni (4 atoms), type=28
            [0.000000000, 0.000000000, 0.000000000],
            [0.500000000, 0.000000000, 0.500000000],
            [0.000000000, 0.500000000, 0.500000000],
            [0.500000000, 0.500000000, 0.000000000],
            // O (16 atoms), type=8
            [0.750000000, 0.750000000, 0.015500000],
            [0.250000000, 0.250000000, 0.484499991],
            [0.250000000, 0.250000000, 0.984499991],
            [0.750000000, 0.750000000, 0.515500009],
            [0.750000000, 0.250000000, 0.000000000],
            [0.250000000, 0.750000000, 0.500000000],
            [0.250000000, 0.750000000, 0.000000000],
            [0.750000000, 0.250000000, 0.500000000],
            [0.031399999, 0.031399999, 0.177100003],
            [0.531400025, 0.968599975, 0.322899997],
            [0.968599975, 0.531400025, 0.322899997],
            [0.468600005, 0.468600005, 0.177100003],
            [0.968599975, 0.968599975, 0.822899997],
            [0.468600005, 0.031399999, 0.677100003],
            [0.031399999, 0.468600005, 0.677100003],
            [0.531400025, 0.531400025, 0.822899997],
        ];

        let mut types = Vec::with_capacity(28);
        for _ in 0..8  { types.push(57_i32); } // La
        for _ in 0..4  { types.push(28_i32); } // Ni
        for _ in 0..16 { types.push(8_i32);  } // O

        eprintln!("=== La₂NiO₄ 手动构造 Cell ===");
        eprintln!("n_atoms = {}", positions.len());

        // 验证四方晶格: a=b≠c, 所有角 90°
        let a = [lattice[0][0], lattice[1][0], lattice[2][0]];
        let b = [lattice[0][1], lattice[1][1], lattice[2][1]];
        let c = [lattice[0][2], lattice[1][2], lattice[2][2]];
        let len_a = (a[0]*a[0] + a[1]*a[1] + a[2]*a[2]).sqrt();
        let len_b = (b[0]*b[0] + b[1]*b[1] + b[2]*b[2]).sqrt();
        let len_c = (c[0]*c[0] + c[1]*c[1] + c[2]*c[2]).sqrt();
        eprintln!("|a| = {:.6}, |b| = {:.6}, |c| = {:.6}", len_a, len_b, len_c);
        assert!((len_a - len_b).abs() < 1e-6, "a must equal b for tetragonal");

        // --- 构造 Cell ---
        let mut cell = Cell::new(28, TensorRank::NoSpin);
        cell.set_cell(&lattice, &positions, &types);
        cell.aperiodic_axis = None;

        // --- 获取原胞 ---
        let primitive = crate::primitive::prm_get_primitive(&cell, 1e-5, -1.0)
            .expect("primitive failed");
        eprintln!("\n原胞原子数: {}", primitive.size);

        // --- 原胞结构 ---
        if let Some(ref prim_cell) = primitive.cell {
            let p_lat = &prim_cell.lattice;
            eprintln!("原胞晶格:");
            eprintln!("  [{:>10.6}, {:>10.6}, {:>10.6}]", p_lat[0][0], p_lat[0][1], p_lat[0][2]);
            eprintln!("  [{:>10.6}, {:>10.6}, {:>10.6}]", p_lat[1][0], p_lat[1][1], p_lat[1][2]);
            eprintln!("  [{:>10.6}, {:>10.6}, {:>10.6}]", p_lat[2][0], p_lat[2][1], p_lat[2][2]);
            for j in 0..3 {
                let len = (p_lat[0][j]*p_lat[0][j] + p_lat[1][j]*p_lat[1][j] + p_lat[2][j]*p_lat[2][j]).sqrt();
                eprintln!("  |v{}| = {:.6}", j, len);
            }
            for j in 0..3 {
                for k in (j+1)..3 {
                    let dot = p_lat[0][j]*p_lat[0][k] + p_lat[1][j]*p_lat[1][k] + p_lat[2][j]*p_lat[2][k];
                    let lj = (p_lat[0][j]*p_lat[0][j] + p_lat[1][j]*p_lat[1][j] + p_lat[2][j]*p_lat[2][j]).sqrt();
                    let lk = (p_lat[0][k]*p_lat[0][k] + p_lat[1][k]*p_lat[1][k] + p_lat[2][k]*p_lat[2][k]).sqrt();
                    let angle = (dot / (lj * lk)).acos() * 180.0 / std::f64::consts::PI;
                    eprintln!("  ∠(v{}, v{}) = {:.4}°", j, k, angle);
                }
            }
            for i in 0..prim_cell.size {
                eprintln!("  atom {:2}: pos=[{:.6}, {:.6}, {:.6}] type={}",
                    i+1, prim_cell.position[i][0], prim_cell.position[i][1],
                    prim_cell.position[i][2], prim_cell.types[i]);
            }
        }

        // --- 获取对称操作 ---
        if let Some(ref prim_cell) = primitive.cell {
            let sym = crate::symmetry::sym_get_operation(prim_cell, 1e-5, -1.0)
                .expect("sym_get_operation failed");
            eprintln!("\n原胞对称操作数: {}", sym.size);

            let pointsym = crate::pointgroup::ptg_get_pointsymmetry(&sym.rot);
            eprintln!("独特旋转操作数: {}", pointsym.size);
            for (i, r) in pointsym.rot.iter().enumerate() {
                let det = crate::mathfunc::mat_get_determinant_i3(r);
                let trace = crate::mathfunc::mat_get_trace_i3(r);
                eprintln!("  rot[{}]: det={}, trace={}", i, det, trace);
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

        let pg = ptg_get_pointgroup(spg.pointgroup_number);
        eprintln!("\n=== 点群 ===");
        eprintln!("  Number:     {}", pg.number);
        eprintln!("  Symbol:     {}", pg.symbol.trim());
        eprintln!("  Schoenflies:{}", pg.schoenflies.trim());

        // --- 断言 ---
        // 空间群 #138 P4₂/ncm (D₄ₕ¹⁶)
        assert_eq!(spg.number, 138,
            "La₂NiO₄ expected P4₂/ncm (#138), got #{}", spg.number);
        assert_eq!(spg.hall_number, 422,
            "expected Hall 422, got {}", spg.hall_number);
        assert_eq!(spg.pointgroup_number, 15,
            "expected point group 15 (4/mmm), got {}", spg.pointgroup_number);
        assert_eq!(spg.international_short.trim(), "P4_2/ncm",
            "expected P4_2/ncm, got '{}'", spg.international_short.trim());
        assert_eq!(spg.schoenflies.trim(), "D4h^16",
            "expected D4h^16, got '{}'", spg.schoenflies.trim());

        // 点群: 4/mmm (D₄ₕ)
        assert_eq!(pg.number, 15, "point group number mismatch");
        assert_eq!(pg.symbol.trim(), "4/mmm",
            "expected 4/mmm, got '{}'", pg.symbol.trim());
        assert_eq!(pg.schoenflies.trim(), "D4h",
            "expected D4h, got '{}'", pg.schoenflies.trim());

        // 原胞: 28 原子（未约化，无纯平移）
        assert_eq!(primitive.size, 28,
            "primitive cell should have 28 atoms, got {}", primitive.size);
        if let Some(ref prim_cell) = primitive.cell {
            // 四方晶格: a=b=5.50, c=12.51
            let p = &prim_cell.lattice;
            let va2 = p[0][0]*p[0][0] + p[1][0]*p[1][0] + p[2][0]*p[2][0];
            let vb2 = p[0][1]*p[0][1] + p[1][1]*p[1][1] + p[2][1]*p[2][1];
            let vc2 = p[0][2]*p[0][2] + p[1][2]*p[1][2] + p[2][2]*p[2][2];
            let va = va2.sqrt();
            let vb = vb2.sqrt();
            let vc = vc2.sqrt();
            assert!((va - vb).abs() < 0.01, "tetragonal: |a| ≈ |b| expected");
            assert!((va - 5.50).abs() < 0.1, "|a| ≈ 5.50, got {:.4}", va);
            assert!((vc - 12.51).abs() < 0.1, "|c| ≈ 12.51, got {:.4}", vc);

            // D₄ₕ 点群: 16 个旋转
            let sym = crate::symmetry::sym_get_operation(prim_cell, 1e-5, -1.0)
                .expect("sym_get_operation failed");
            assert_eq!(sym.size, 16,
                "expected 16 symmetry ops (D4h), got {}", sym.size);

            let pointsym = crate::pointgroup::ptg_get_pointsymmetry(&sym.rot);
            assert_eq!(pointsym.size, 16,
                "expected 16 unique rotations (D4h), got {}", pointsym.size);
        }
    }
}
