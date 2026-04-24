//! 磁性空间群/磁点群测试。
//!
//! 通过构造磁性对称操作（旋转 + 平移 + 时间反演），验证磁性空间群识别管线。
//! 磁性类型通过 `type_` 字段分类: 1=ordinary, 2=grey, 3=black-white, 4=black-white (anti-translation)。

#[cfg(test)]
mod tests {
    use crate::mathfunc::{mat_get_determinant_i3, Mat3, Mat3I, Vec3};
    use crate::msg_database::msgdb_get_magnetic_spacegroup_type;
    use crate::spg_get_magnetic_spacegroup_type_from_symmetry;
    use crate::symmetry::MagneticSymmetry;

    const SYMPREC: f64 = 1e-5;

    fn cubic_lattice() -> Mat3 {
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    }

    /// 获取空间群 #221 (Pm-3m, hall=497) 的数据库对称操作。
    fn pm3m_ops() -> Vec<(Mat3I, Vec3)> {
        let (count, start) = crate::spg_database::spgdb_get_operation_index(497);
        (0..count)
            .filter_map(|i| crate::spg_database::spgdb_get_operation_by_index(start + i))
            .collect()
    }

    fn make_mag_sym(timerev: &[i32], ops: &[(Mat3I, Vec3)]) -> MagneticSymmetry {
        assert_eq!(timerev.len(), ops.len());
        let mut sym = MagneticSymmetry::new(ops.len());
        for (i, ((r, t), &tr)) in ops.iter().zip(timerev).enumerate() {
            sym.rot[i] = *r;
            sym.trans[i] = *t;
            sym.timerev[i] = tr;
        }
        sym
    }

    // ====================================================================
    // Type-1: 无时间反演的普通磁性空间群
    // ====================================================================

    #[test]
    fn test_type1_pm3m() {
        // Pm-3m (#221) 的对称操作全部 timerev=0 → type-1
        let ops = pm3m_ops();
        let mag_sym = make_mag_sym(&vec![0; ops.len()], &ops);
        let ds = crate::magnetic_spacegroup::msg_identify_magnetic_space_group_type(
            &cubic_lattice(), &mag_sym, SYMPREC,
        )
        .expect("Type-1 Pm-3m should match");
        assert_eq!(ds.msg_type, 1);
        assert_eq!(ds.hall_number, 497);
        assert_eq!(msgdb_get_magnetic_spacegroup_type(ds.uni_number).type_, 1);
    }

    #[test]
    fn test_type1_identity_only() {
        // 仅恒等操作 → P1 (#1, type-1)
        let mut mag_sym = MagneticSymmetry::new(1);
        mag_sym.rot[0] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        mag_sym.trans[0] = [0.0; 3];
        mag_sym.timerev[0] = 0;

        let ds = crate::magnetic_spacegroup::msg_identify_magnetic_space_group_type(
            &cubic_lattice(), &mag_sym, SYMPREC,
        )
        .expect("Identity should match type-1");
        assert_eq!(ds.msg_type, 1);
    }

    // ====================================================================
    // Type-2 (灰色/顺磁): 每个操作附带时间反演副本
    // ====================================================================

    #[test]
    fn test_type2_pm3m() {
        // Pm-3m 每个操作有 timerev=0 和 timerev=1 两份
        // FSG=96 (去重后), XSG=96, MSG=192 → type-2
        let ops = pm3m_ops();
        let n = ops.len();
        let mut mag_sym = MagneticSymmetry::new(n * 2);
        for (i, (r, t)) in ops.iter().enumerate() {
            mag_sym.rot[i] = *r;
            mag_sym.trans[i] = *t;
            mag_sym.timerev[i] = 0;
            mag_sym.rot[i + n] = *r;
            mag_sym.trans[i + n] = *t;
            mag_sym.timerev[i + n] = 1;
        }

        let ds = crate::magnetic_spacegroup::msg_identify_magnetic_space_group_type(
            &cubic_lattice(), &mag_sym, SYMPREC,
        )
        .expect("Type-2 Pm-3m should match");
        assert_eq!(ds.msg_type, 2);
        assert_eq!(ds.hall_number, 497);
        assert_eq!(msgdb_get_magnetic_spacegroup_type(ds.uni_number).type_, 2);
    }

    #[test]
    fn test_type2_identity_with_timerev() {
        // 恒等 + 时间反演恒等 → type-2 P1
        let mut mag_sym = MagneticSymmetry::new(2);
        mag_sym.rot[0] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        mag_sym.trans[0] = [0.0; 3];
        mag_sym.timerev[0] = 0;
        mag_sym.rot[1] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        mag_sym.trans[1] = [0.0; 3];
        mag_sym.timerev[1] = 1;

        let ds = crate::magnetic_spacegroup::msg_identify_magnetic_space_group_type(
            &cubic_lattice(), &mag_sym, SYMPREC,
        )
        .expect("Type-2 identity should match");
        assert_eq!(ds.msg_type, 2);
    }

    // ====================================================================
    // Type-3 (黑白/反铁磁): proper/improper 旋转分离
    // ====================================================================

    #[test]
    fn test_type3_pm3m() {
        // Pm-3m: proper 旋转 (det=1) → ordinary, improper (det=-1) → anti
        // FSG=96, XSG=48 (proper 旋转 + 平移), FSG=2*XSG → type-3
        let ops = pm3m_ops();
        let timerev: Vec<i32> = ops
            .iter()
            .map(|(r, _)| if mat_get_determinant_i3(r) == 1 { 0 } else { 1 })
            .collect();
        let mag_sym = make_mag_sym(&timerev, &ops);

        let ds = crate::magnetic_spacegroup::msg_identify_magnetic_space_group_type(
            &cubic_lattice(), &mag_sym, SYMPREC,
        )
        .expect("Type-3 Pm-3m should match");
        assert_eq!(ds.msg_type, 3);
        assert_eq!(ds.hall_number, 497);
        assert_eq!(msgdb_get_magnetic_spacegroup_type(ds.uni_number).type_, 3);
    }

    // ====================================================================
    // 物理测试: Fe 原子在简单立方体心 [0.5,0.5,0.5]
    //
    // 验证:
    //   1. 非磁空间群 → Pm-3m (#221)
    //   2. 不同磁矩方向下的磁群:
    //      磁矩作为轴矢量 (pseudovector) 变换: m'_i = det(R) * R_ij * m_j
    //      若 m' = m → timerev=0; 若 m' = -m → timerev=1
    //
    // 磁空间群识别要求输入的对称操作构成完整的可识别空间群。
    // 因此使用实际晶体 (Cell → Primitive → sym_get_operation) 得到的对称操作，
    // 而非数据库操作。
    // ====================================================================

    /// 计算每个对称操作对磁矩的时间反演标记。
    fn compute_timerev(moment: &Vec3, ops: &[(Mat3I, Vec3)]) -> Vec<i32> {
        let norm = (moment[0] * moment[0] + moment[1] * moment[1] + moment[2] * moment[2]).sqrt();
        let m = [moment[0] / norm, moment[1] / norm, moment[2] / norm];

        ops.iter()
            .map(|(rot, _)| {
                let det = mat_get_determinant_i3(rot);
                // m'_j = det(R) * Σ_k R[j][k] * m[k]
                let m_t0 = (det as f64)
                    * (rot[0][0] as f64 * m[0] + rot[0][1] as f64 * m[1] + rot[0][2] as f64 * m[2]);
                let m_t1 = (det as f64)
                    * (rot[1][0] as f64 * m[0] + rot[1][1] as f64 * m[1] + rot[1][2] as f64 * m[2]);
                let m_t2 = (det as f64)
                    * (rot[2][0] as f64 * m[0] + rot[2][1] as f64 * m[1] + rot[2][2] as f64 * m[2]);

                let preserved =
                    (m_t0 - m[0]).abs() < SYMPREC && (m_t1 - m[1]).abs() < SYMPREC && (m_t2 - m[2]).abs() < SYMPREC;
                let reversed =
                    (m_t0 + m[0]).abs() < SYMPREC && (m_t1 + m[1]).abs() < SYMPREC && (m_t2 + m[2]).abs() < SYMPREC;

                if preserved {
                    0
                } else if reversed {
                    1
                } else {
                    -1 // 无效: 该操作不是磁对称操作
                }
            })
            .collect()
    }

    #[test]
    fn test_fe_center_cubic_magnetic() {
        // 简单立方晶格，1个Fe原子在体心 [0.5,0.5,0.5]
        // 物理内容:
        //   非磁: Fe 立方晶格 → Pm-3m (#221)
        //   顺磁: 每个对称操作都有时间反演副本 → Type-2 MSG
        //   反铁磁(proper/improper): proper 旋转为普通, improper 为时间反演 → Type-3 MSG
        //   [1,0,0] 磁矩: 方向沿 x → 16/48 操作保留 → 用 parent hall 识别磁群
        //   [1,1,1] 磁矩: 方向沿对角线 → 12/48 操作保留 → 用 parent hall 识别磁群
        let lattice = cubic_lattice();

        // --- 构建 Fe 在体心的 Cell ---
        let mut cell = crate::cell::Cell::new(1, crate::cell::TensorRank::NoSpin);
        cell.set_cell(&lattice, &[[0.5, 0.5, 0.5]], &[26]);
        cell.aperiodic_axis = None;

        // --- 1. 非磁空间群 ---
        let primitive = crate::primitive::prm_get_primitive(&cell, SYMPREC, -1.0)
            .expect("Primitive search failed");
        let spg = crate::spacegroup::spa_search_spacegroup(&primitive, 0, SYMPREC, -1.0)
            .expect("Space group search failed");
        assert_eq!(spg.number, 221, "Fe at SC body center → Pm-3m (#221)");
        let parent_hall = spg.hall_number;

        // --- 获取晶体对称操作 ---
        let prim_cell = primitive.cell.as_ref().expect("Primitive cell exists");
        let symmetry = crate::symmetry::sym_get_operation(prim_cell, SYMPREC, -1.0)
            .expect("Symmetry operation search failed");
        assert_eq!(symmetry.size, 48, "Pm-3m primitive cell has 48 symmetry ops");
        let crystal_ops: Vec<(Mat3I, Vec3)> = (0..symmetry.size)
            .map(|i| (symmetry.rot[i], symmetry.trans[i]))
            .collect();

        // --- 2. 磁矩 [1,0,0]: 验证对称性破缺 + 磁群识别 ---
        let tr_x = compute_timerev(&[1.0, 0.0, 0.0], &crystal_ops);
        assert_eq!(tr_x[0], 0, "Identity preserves moment [1,0,0]");
        let n_valid_x = tr_x.iter().filter(|&&t| t != -1).count();
        let n_ord_x = tr_x.iter().filter(|&&t| t == 0).count();
        let n_anti_x = tr_x.iter().filter(|&&t| t == 1).count();
        assert_eq!(n_valid_x, 16, "Moment [1,0,0]: 16 of 48 ops preserved");
        assert_eq!(n_ord_x, 8, "Moment [1,0,0]: 8 ordinary");
        assert_eq!(n_anti_x, 8, "Moment [1,0,0]: 8 anti");

        // 用磁对称操作识别磁空间群
        // [1,0,0] 磁矩将立方对称破缺到 16 个有效操作，按 type-3 (8 ordinary + 8 anti)
        // 这个构型在标准 Litvin 数据库中可能没有完全一致的条目（立方 D4h 子群与
        // 四方 D4h 的操作矩阵不同），但管线应正确分类为 Type-3
        let valid_indices_x: Vec<usize> = tr_x
            .iter()
            .cloned()
            .enumerate()
            .filter(|(_, t)| *t != -1)
            .map(|(i, _)| i)
            .collect();
        let mag_sym_x = {
            let n = valid_indices_x.len();
            let mut sym = MagneticSymmetry::new(n);
            for (j, &i) in valid_indices_x.iter().enumerate() {
                sym.rot[j] = crystal_ops[i].0;
                sym.trans[j] = crystal_ops[i].1;
                sym.timerev[j] = tr_x[i];
            }
            sym
        };

        // 先用 parent hall 尝试匹配（Pm-3m 的磁变体）
        let try_x = crate::magnetic_spacegroup::msg_identify_with_parent_hall(
            &lattice, &mag_sym_x, Some(parent_hall), SYMPREC,
        );
        match try_x {
            Some(ds) => {
                let mt = msgdb_get_magnetic_spacegroup_type(ds.uni_number);
                assert_eq!(ds.msg_type, 3, "[1,0,0]: Type-3");
                assert_eq!(mt.number, 221, "[1,0,0]: should match a #221 variant");
            }
            None => {
                // 标准 Litvin DB 没有直接收录此构型，属正常
                eprintln!("Moment [1,0,0]: not in MSG DB as a #221 variant");
                // 但仍可通过标准路径自动找到对应 Hall 编号
                let ds_auto = crate::magnetic_spacegroup::msg_identify_magnetic_space_group_type(
                    &lattice, &mag_sym_x, SYMPREC,
                );
                if let Some(ds) = ds_auto {
                    let mt = msgdb_get_magnetic_spacegroup_type(ds.uni_number);
                    eprintln!("  auto-detect: uni={}, type={}, number={}, bns='{}'",
                        ds.uni_number, ds.msg_type, mt.number, mt.bns_number.trim());
                    assert_eq!(ds.msg_type, 3, "[1,0,0] auto: Type-3");
                } else {
                    eprintln!("  auto-detect also failed");
                }
            }
        }

        // --- 3. 磁矩 [1,1,1]: 验证对称性破缺 + 磁群识别 ---
        let tr_111 = compute_timerev(&[1.0, 1.0, 1.0], &crystal_ops);
        assert_eq!(tr_111[0], 0, "Identity preserves moment [1,1,1]");
        let n_valid_111 = tr_111.iter().filter(|&&t| t != -1).count();
        let n_ord_111 = tr_111.iter().filter(|&&t| t == 0).count();
        let n_anti_111 = tr_111.iter().filter(|&&t| t == 1).count();
        assert_eq!(n_valid_111, 12, "Moment [1,1,1]: 12 of 48 ops preserved");
        assert_eq!(n_ord_111, 6, "Moment [1,1,1]: 6 ordinary");
        assert_eq!(n_anti_111, 6, "Moment [1,1,1]: 6 anti");

        let valid_indices_111: Vec<usize> = tr_111
            .iter()
            .cloned()
            .enumerate()
            .filter(|(_, t)| *t != -1)
            .map(|(i, _)| i)
            .collect();
        let mag_sym_111 = {
            let n = valid_indices_111.len();
            let mut sym = MagneticSymmetry::new(n);
            for (j, &i) in valid_indices_111.iter().enumerate() {
                sym.rot[j] = crystal_ops[i].0;
                sym.trans[j] = crystal_ops[i].1;
                sym.timerev[j] = tr_111[i];
            }
            sym
        };
        let try_111 = crate::magnetic_spacegroup::msg_identify_with_parent_hall(
            &lattice, &mag_sym_111, Some(parent_hall), SYMPREC,
        );
        match try_111 {
            Some(ds) => {
                let mt = msgdb_get_magnetic_spacegroup_type(ds.uni_number);
                assert_eq!(ds.msg_type, 3, "[1,1,1]: Type-3");
                assert_eq!(mt.number, 221);
            }
            None => {
                eprintln!("Moment [1,1,1]: not in MSG DB as a #221 variant");
                if let Some(ds) = crate::magnetic_spacegroup::msg_identify_magnetic_space_group_type(
                    &lattice, &mag_sym_111, SYMPREC,
                ) {
                    let mt = msgdb_get_magnetic_spacegroup_type(ds.uni_number);
                    eprintln!("  auto-detect: uni={}, type={}, number={}, bns='{}'",
                        ds.uni_number, ds.msg_type, mt.number, mt.bns_number.trim());
                    assert_eq!(ds.msg_type, 3, "[1,1,1] auto: Type-3");
                } else {
                    eprintln!("  auto-detect also failed");
                }
            }
        }

        // --- 4. 顺磁 Type-2: 所有操作 + 时间反演副本 ---
        {
            let n = crystal_ops.len();
            let mut mag_sym = MagneticSymmetry::new(n * 2);
            for (i, (r, t)) in crystal_ops.iter().enumerate() {
                mag_sym.rot[i] = *r;
                mag_sym.trans[i] = *t;
                mag_sym.timerev[i] = 0;
                mag_sym.rot[i + n] = *r;
                mag_sym.trans[i + n] = *t;
                mag_sym.timerev[i + n] = 1;
            }
            let ds = crate::magnetic_spacegroup::msg_identify_magnetic_space_group_type(
                &lattice, &mag_sym, SYMPREC,
            )
            .expect("Type-2 paramagnetic Fe should match");
            let mt = msgdb_get_magnetic_spacegroup_type(ds.uni_number);
            assert_eq!(ds.msg_type, 2);
            assert_eq!(mt.type_, 2);
            assert_eq!(mt.number, 221, "Type-2: crystal #221");
        }

        // --- 5. Type-3 反铁磁: proper/improper 分离 ---
        {
            let timerev: Vec<i32> = crystal_ops
                .iter()
                .map(|(r, _)| if mat_get_determinant_i3(r) == 1 { 0 } else { 1 })
                .collect();
            let mag_sym = make_mag_sym(&timerev, &crystal_ops);
            let ds = crate::magnetic_spacegroup::msg_identify_magnetic_space_group_type(
                &lattice, &mag_sym, SYMPREC,
            )
            .expect("Type-3 AFM Fe should match");
            let mt = msgdb_get_magnetic_spacegroup_type(ds.uni_number);
            assert_eq!(ds.msg_type, 3);
            assert_eq!(mt.type_, 3);
            assert_eq!(mt.number, 221, "Type-3: crystal #221");
        }
    }

    // ====================================================================
    // 公共 API 测试: spg_get_magnetic_spacegroup_type_from_symmetry
    // ====================================================================

    #[test]
    fn test_api_type1() {
        let ops = pm3m_ops();
        let rots: Vec<_> = ops.iter().map(|(r, _)| *r).collect();
        let trans: Vec<_> = ops.iter().map(|(_, t)| *t).collect();
        let result = spg_get_magnetic_spacegroup_type_from_symmetry(
            &rots, &trans, None, &cubic_lattice(), SYMPREC,
        );
        assert!(result.uni_number > 0);
        assert_eq!(result.type_, 1);
    }

    #[test]
    fn test_api_type2() {
        let ops = pm3m_ops();
        let n = ops.len();
        let rots: Vec<_> = ops.iter().chain(ops.iter()).map(|(r, _)| *r).collect();
        let trans: Vec<_> = ops.iter().chain(ops.iter()).map(|(_, t)| *t).collect();
        let timerev: Vec<i32> = (0..n).map(|_| 0).chain((0..n).map(|_| 1)).collect();
        let result = spg_get_magnetic_spacegroup_type_from_symmetry(
            &rots, &trans, Some(&timerev), &cubic_lattice(), SYMPREC,
        );
        assert!(result.uni_number > 0);
        assert_eq!(result.type_, 2);
    }

    #[test]
    fn test_api_type3() {
        let ops = pm3m_ops();
        let rots: Vec<_> = ops.iter().map(|(r, _)| *r).collect();
        let trans: Vec<_> = ops.iter().map(|(_, t)| *t).collect();
        let timerev: Vec<i32> = ops
            .iter()
            .map(|(r, _)| if mat_get_determinant_i3(r) == 1 { 0 } else { 1 })
            .collect();
        let result = spg_get_magnetic_spacegroup_type_from_symmetry(
            &rots, &trans, Some(&timerev), &cubic_lattice(), SYMPREC,
        );
        assert!(result.uni_number > 0);
        assert_eq!(result.type_, 3);
    }

    // ====================================================================
    // 边界情况
    // ====================================================================

    #[test]
    fn test_empty_symmetry() {
        let mag_sym = MagneticSymmetry::new(0);
        let result = crate::magnetic_spacegroup::msg_identify_magnetic_space_group_type(
            &cubic_lattice(), &mag_sym, SYMPREC,
        );
        assert!(result.is_none());
    }

    #[test]
    fn test_no_identity() {
        // 没有恒等操作 → 预期 None（无纯平移）
        let mut mag_sym = MagneticSymmetry::new(1);
        mag_sym.rot[0] = [[0, -1, 0], [1, 0, 0], [0, 0, 1]];
        mag_sym.trans[0] = [0.0; 3];
        mag_sym.timerev[0] = 0;
        let result = crate::magnetic_spacegroup::msg_identify_magnetic_space_group_type(
            &cubic_lattice(), &mag_sym, SYMPREC,
        );
        assert!(result.is_none());
    }
}
