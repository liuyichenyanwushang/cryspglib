//! CrPS₄ 空间群测试。
//!
//! 单斜晶系 C2 (#5)。POSCAR 含 8 Cr + 8 P + 32 S = 48 原子。
//! PPOSCAR 原胞含 2 Cr + 2 P + 8 S = 12 原子（约化因子 4）。

#[cfg(test)]
mod test_crps4 {

    use crate::cell::{Cell, TensorRank};
    use crate::mathfunc::Mat3;
    use crate::pointgroup::ptg_get_pointgroup;

    /// PPOSCAR 给出的预期原胞晶格，用于验证原胞约化结果与 phonopy 一致。
    const EXPECTED_PRIM_LATTICE: Mat3 = [
        [ 5.4254999161,  5.4254999161, -0.2169236918],
        [-3.6280000210,  3.6280000210,  0.0000000000],
        [ 0.0000000000,  0.0000000000,  6.1019455591],
    ];

    /// 手动硬编码 CrPS₄ POSCAR 数据，构造 Cell、约化原胞并与 PPOSCAR 对比。
    ///
    /// 原胞晶格与 PPOSCAR 一致（允许列交换和符号翻转）。
    /// 空间群 C2 (#5)，点群 2 (C₂)。
    #[test]
    fn test_crps4_hardcoded() {
        // --- POSCAR 手动写入 ---
        // 单斜晶格: a=10.851, b=7.256, c=12.204, β≈92.0°
        // POSCAR 行向量 → lattice[cart][vec]（列=晶格矢量）
        let lattice: Mat3 = [
            [10.8509998322,  0.0000000000, -0.4338473836], // x 分量
            [ 0.0000000000,  7.2560000420,  0.0000000000], // y 分量
            [ 0.0000000000,  0.0000000000, 12.2038911182], // z 分量
        ];

        // 8 Cr + 8 P + 32 S，POSCAR Direct 坐标
        let positions = [
            // Cr (8 atoms), type=24
            [0.000000000, 0.037500001, 0.000000000],
            [0.000000000, 0.037500001, 0.500000000],
            [0.500000000, 0.537500024, 0.000000000],
            [0.500000000, 0.537500024, 0.500000000],
            [0.000000000, 0.542699993, 0.000000000],
            [0.000000000, 0.542699993, 0.500000000],
            [0.500000000, 0.042699993, 0.000000000],
            [0.500000000, 0.042699993, 0.500000000],
            // P (8 atoms), type=15
            [0.294600010, 0.273799986, 0.082050003],
            [0.705399990, 0.273799986, 0.417950004],
            [0.294600010, 0.273799986, 0.582050025],
            [0.705399990, 0.273799986, 0.917949975],
            [0.794600010, 0.773800015, 0.082050003],
            [0.205399990, 0.773800015, 0.417950004],
            [0.794600010, 0.773800015, 0.582050025],
            [0.205399990, 0.773800015, 0.917949975],
            // S (32 atoms), type=16
            [0.143000007, 0.000000000, 0.346799999],
            [0.856999993, 0.000000000, 0.153200001],
            [0.143000007, 0.000000000, 0.846799970],
            [0.856999993, 0.000000000, 0.653200030],
            [0.643000007, 0.500000000, 0.346799999],
            [0.356999993, 0.500000000, 0.153200001],
            [0.643000007, 0.500000000, 0.846799970],
            [0.356999993, 0.500000000, 0.653200030],
            [0.127900004, 0.530399978, 0.351249993],
            [0.872099996, 0.530399978, 0.148750007],
            [0.127900004, 0.530399978, 0.851249993],
            [0.872099996, 0.530399978, 0.648750007],
            [0.627900004, 0.030399978, 0.351249993],
            [0.372099996, 0.030399978, 0.148750007],
            [0.627900004, 0.030399978, 0.851249993],
            [0.372099996, 0.030399978, 0.648750007],
            [0.104199998, 0.284500003, 0.097750001],
            [0.895799994, 0.284500003, 0.402249992],
            [0.104199998, 0.284500003, 0.597750008],
            [0.895799994, 0.284500003, 0.902249992],
            [0.604200006, 0.784500003, 0.097750001],
            [0.395799994, 0.784500003, 0.402249992],
            [0.604200006, 0.784500003, 0.597750008],
            [0.395799994, 0.784500003, 0.902249992],
            [0.127100006, 0.761200011, 0.076800004],
            [0.872900009, 0.761200011, 0.423200011],
            [0.127100006, 0.761200011, 0.576799989],
            [0.872900009, 0.761200011, 0.923200011],
            [0.627099991, 0.261200011, 0.076800004],
            [0.372900009, 0.261200011, 0.423200011],
            [0.627099991, 0.261200011, 0.576799989],
            [0.372900009, 0.261200011, 0.923200011],
        ];

        let mut types = Vec::with_capacity(48);
        for _ in 0..8  { types.push(24_i32); } // Cr
        for _ in 0..8  { types.push(15_i32); } // P
        for _ in 0..32 { types.push(16_i32); } // S

        // --- 构造 Cell ---
        let mut cell = Cell::new(48, TensorRank::NoSpin);
        cell.set_cell(&lattice, &positions, &types);
        cell.aperiodic_axis = None;

        // --- 获取原胞 ---
        let primitive = crate::primitive::prm_get_primitive(&cell, 1e-5, -1.0)
            .expect("prm_get_primitive failed");
        let prim_cell = primitive.cell.as_ref().expect("primitive cell exists");

        // --- 断言原胞 ---
        // 48→12 原子（约化因子 4：C 底心 + 沿 c 平移）
        assert_eq!(prim_cell.size, 12,
            "primitive should be 12 atoms, got {}", prim_cell.size);

        // --- 对比 PPOSCAR 原胞晶格 ---
        let p = &prim_cell.lattice;
        let e = &EXPECTED_PRIM_LATTICE;
        // 验证原胞晶格列矢量与 PPOSCAR 一致（允许符号翻转和列交换）
        for c in 0..3 {
            let mut found = false;
            for ec in 0..3 {
                let d0 = p[0][c] - e[0][ec];
                let d1 = p[1][c] - e[1][ec];
                let d2 = p[2][c] - e[2][ec];
                let neg_d0 = p[0][c] + e[0][ec];
                let neg_d1 = p[1][c] + e[1][ec];
                let neg_d2 = p[2][c] + e[2][ec];
                let norm = (d0*d0 + d1*d1 + d2*d2).sqrt();
                let neg_norm = (neg_d0*neg_d0 + neg_d1*neg_d1 + neg_d2*neg_d2).sqrt();
                if norm < 0.01 || neg_norm < 0.01 {
                    found = true;
                    break;
                }
            }
            assert!(found, "primitive lattice column {} not in PPOSCAR (up to sign)", c);
        }

        // --- 获取对称操作 ---
        let sym = crate::symmetry::sym_get_operation(prim_cell, 1e-5, -1.0)
            .expect("sym_get_operation failed");
        // C₂ 点群: identity + C₂ 旋转 = 2 个操作
        assert_eq!(sym.size, 2,
            "expected 2 symmetry ops (C2), got {}", sym.size);

        let pointsym = crate::pointgroup::ptg_get_pointsymmetry(&sym.rot);
        assert_eq!(pointsym.size, 2,
            "expected 2 unique rotations, got {}", pointsym.size);

        // --- 搜索空间群 ---
        let spg = crate::spacegroup::spa_search_spacegroup(&primitive, 0, 1e-5, -1.0)
            .expect("spacegroup search failed");

        // C2 (#5), Hall 9, 点群 2 (C₂)
        assert_eq!(spg.number, 5,
            "CrPS₄ expected C2 (#5), got #{}", spg.number);
        assert_eq!(spg.hall_number, 9,
            "expected Hall 9, got {}", spg.hall_number);
        assert_eq!(spg.pointgroup_number, 3,
            "expected point group 3 (2), got {}", spg.pointgroup_number);
        assert_eq!(spg.international_short.trim(), "C2",
            "expected C2, got '{}'", spg.international_short.trim());
        assert_eq!(spg.schoenflies.trim(), "C2^3",
            "expected C2^3, got '{}'", spg.schoenflies.trim());

        // 点群: 2 (C₂)
        let pg = ptg_get_pointgroup(spg.pointgroup_number);
        assert_eq!(pg.number, 3, "point group number mismatch");
        assert_eq!(pg.symbol.trim(), "2",
            "expected 2, got '{}'", pg.symbol.trim());
        assert_eq!(pg.schoenflies.trim(), "C2",
            "expected C2, got '{}'", pg.schoenflies.trim());
    }
}
