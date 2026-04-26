//! 磁性空间群测试。
//!
//! 所有测试走公共 API `spg_get_magnetic_dataset`，不绕过真实管线。
//! 数据库层测试保留（验证 `msg_identify_magnetic_space_group_type` 算法本身）。

#[cfg(test)]
mod tests {
    use crate::MagneticType;
    use crate::mathfunc::{Mat3, Mat3I, Vec3};
    use crate::msg_database::msgdb_get_magnetic_spacegroup_type;
    use crate::spg_get_magnetic_spacegroup_type_from_symmetry;
    use crate::symmetry::MagneticSymmetry;

    const SYMPREC: f64 = 1e-5;

    fn cubic_lattice() -> Mat3 {
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    }

    /// 运行公共 API 并打印结果。
    fn run_dataset(
        label: &str,
        lattice: &Mat3,
        positions: &[Vec3],
        types: &[i32],
        moments: Option<&[[f64; 3]]>,
    ) -> crate::SpglibMagneticSymmetry {
        let result = crate::spg_get_magnetic_dataset(lattice, positions, types, moments, SYMPREC)
            .unwrap_or_else(|| panic!("{}: spg_get_magnetic_dataset returned None", label));
        eprintln!("=== {} ===", label);
        eprintln!("{}", crate::spg_format_magnetic_symmetry(&result));
        result
    }

    // ====================================================================
    // 数据库层: 用 Pm-3m 数据库操作直接测试识别算法
    // ====================================================================

    fn pm3m_ops() -> Vec<(Mat3I, Vec3)> {
        let (count, start) = crate::spg_database::spgdb_get_operation_index(517);
        (0..count)
            .filter_map(|i| crate::spg_database::spgdb_get_operation_by_index(start + i))
            .collect()
    }

    fn make_mag_sym(timerev: &[bool], ops: &[(Mat3I, Vec3)]) -> MagneticSymmetry {
        assert_eq!(timerev.len(), ops.len());
        let mut sym = MagneticSymmetry::new(ops.len());
        for (i, ((r, t), &tr)) in ops.iter().zip(timerev).enumerate() {
            sym.rot[i] = *r;
            sym.trans[i] = *t;
            sym.timerev[i] = tr;
        }
        sym
    }

    #[test]
    fn test_db_type1() {
        let ops = pm3m_ops();
        let mag_sym = make_mag_sym(&vec![false; ops.len()], &ops);
        let ds = crate::magnetic_spacegroup::msg_identify_magnetic_space_group_type(
            &cubic_lattice(), &mag_sym, SYMPREC,
        )
        .expect("must match");
        assert_eq!(ds.msg_type, MagneticType::Ordinary);
        assert_eq!(ds.hall_number, 517);
        assert_eq!(msgdb_get_magnetic_spacegroup_type(ds.uni_number).type_, MagneticType::Ordinary);
    }

    #[test]
    fn test_db_type2() {
        let ops = pm3m_ops();
        let n = ops.len();
        let mut mag_sym = MagneticSymmetry::new(n * 2);
        for (i, (r, t)) in ops.iter().enumerate() {
            mag_sym.rot[i] = *r;
            mag_sym.trans[i] = *t;
            mag_sym.timerev[i] = false;
            mag_sym.rot[i + n] = *r;
            mag_sym.trans[i + n] = *t;
            mag_sym.timerev[i + n] = true;
        }
        let ds = crate::magnetic_spacegroup::msg_identify_magnetic_space_group_type(
            &cubic_lattice(), &mag_sym, SYMPREC,
        )
        .expect("must match");
        assert_eq!(ds.msg_type, MagneticType::Grey);
        assert_eq!(ds.hall_number, 517);
        assert_eq!(msgdb_get_magnetic_spacegroup_type(ds.uni_number).type_, MagneticType::Grey);
    }

    #[test]
    fn test_db_type3() {
        let ops = pm3m_ops();
        let timerev: Vec<bool> = ops.iter().map(|(r, _)| !crate::mathfunc::is_proper(r)).collect();
        let mag_sym = make_mag_sym(&timerev, &ops);
        let ds = crate::magnetic_spacegroup::msg_identify_magnetic_space_group_type(
            &cubic_lattice(), &mag_sym, SYMPREC,
        )
        .expect("must match");
        assert_eq!(ds.msg_type, MagneticType::BlackWhite);
        assert_eq!(ds.hall_number, 517);
        assert_eq!(msgdb_get_magnetic_spacegroup_type(ds.uni_number).type_, MagneticType::BlackWhite);
    }

    // ====================================================================
    // 公共 API: spg_get_magnetic_spacegroup_type_from_symmetry
    // ====================================================================

    #[test]
    fn test_api_type1() {
        let ops = pm3m_ops();
        let rots: Vec<_> = ops.iter().map(|(r, _)| *r).collect();
        let trans: Vec<_> = ops.iter().map(|(_, t)| *t).collect();
        let result = spg_get_magnetic_spacegroup_type_from_symmetry(
            &rots, &trans, None, &cubic_lattice(), SYMPREC,
        );
        assert_eq!(result.type_, MagneticType::Ordinary);
        assert!(result.uni_number > 0);
    }

    #[test]
    fn test_api_type2() {
        let ops = pm3m_ops();
        let n = ops.len();
        let rots: Vec<_> = ops.iter().chain(ops.iter()).map(|(r, _)| *r).collect();
        let trans: Vec<_> = ops.iter().chain(ops.iter()).map(|(_, t)| *t).collect();
        let timerev: Vec<bool> = (0..n).map(|_| false).chain((0..n).map(|_| true)).collect();
        let result = spg_get_magnetic_spacegroup_type_from_symmetry(
            &rots, &trans, Some(&timerev), &cubic_lattice(), SYMPREC,
        );
        assert_eq!(result.type_, MagneticType::Grey);
        assert!(result.uni_number > 0);
    }

    #[test]
    fn test_api_type3() {
        let ops = pm3m_ops();
        let rots: Vec<_> = ops.iter().map(|(r, _)| *r).collect();
        let trans: Vec<_> = ops.iter().map(|(_, t)| *t).collect();
        let timerev: Vec<bool> = ops.iter().map(|(r, _)| !crate::mathfunc::is_proper(r)).collect();
        let result = spg_get_magnetic_spacegroup_type_from_symmetry(
            &rots, &trans, Some(&timerev), &cubic_lattice(), SYMPREC,
        );
        assert_eq!(result.type_, MagneticType::BlackWhite);
        assert!(result.uni_number > 0);
    }

    // ====================================================================
    // 边界情况
    // ====================================================================

    #[test]
    fn test_empty_symmetry() {
        let mag_sym = MagneticSymmetry::new(0);
        assert!(crate::magnetic_spacegroup::msg_identify_magnetic_space_group_type(
            &cubic_lattice(), &mag_sym, SYMPREC,
        ).is_none());
    }

    #[test]
    fn test_no_identity() {
        let mut mag_sym = MagneticSymmetry::new(1);
        mag_sym.rot[0] = [[0, -1, 0], [1, 0, 0], [0, 0, 1]];
        mag_sym.trans[0] = [0.0; 3];
        mag_sym.timerev[0] = false;
        assert!(crate::magnetic_spacegroup::msg_identify_magnetic_space_group_type(
            &cubic_lattice(), &mag_sym, SYMPREC,
        ).is_none());
    }

    // ====================================================================
    // 物理测试: 通过公共 API，POSCAR 式输入
    // ====================================================================

    /// Fe SC 体心, 非磁 (moments=None)
    #[test]
    fn test_fe_sc_nonmagnetic() {
        let lattice = cubic_lattice();
        let positions = [[0.5, 0.5, 0.5]];
        let types = [26];
        let r = run_dataset("Fe SC, non-magnetic", &lattice, &positions, &types, None);
        assert_eq!(r.spacegroup_number, 221);
        assert_eq!(r.magnetic_type, MagneticType::NonMagnetic);
        assert_eq!(r.uni_number, 0);
        assert!(r.num_operations > 0);
    }

    /// Fe SC 体心, 磁矩 [001] → P4/mmm (#123) BNS=123.345, UNI=1005
    #[test]
    fn test_fe_sc_001() {
        let lattice = cubic_lattice();
        let positions = [[0.5, 0.5, 0.5]];
        let types = [26];
        let moments = [[0.0, 0.0, 1.0]];
        let r = run_dataset("Fe SC [001]", &lattice, &positions, &types, Some(&moments));
        assert_eq!(r.spacegroup_number, 221, "non-mag: Pm-3m");
        assert_eq!(r.magnetic_type, MagneticType::BlackWhite);
        assert!(r.uni_number > 0, "must match DB entry, not fallback");
        assert!(!r.bns_number.is_empty(), "must have BNS");
        assert_eq!(r.bns_number.trim(), "123.345");
        assert_eq!(r.uni_number, 1005);
    }

    /// Fe SC 体心, 磁矩 [100] → 应与 [001] 一样匹配 UNI=1005
    #[test]
    fn test_fe_sc_100() {
        let lattice = cubic_lattice();
        let positions = [[0.5, 0.5, 0.5]];
        let types = [26];
        let moments = [[1.0, 0.0, 0.0]];
        let r = run_dataset("Fe SC [100]", &lattice, &positions, &types, Some(&moments));
        assert_eq!(r.spacegroup_number, 221, "non-mag: Pm-3m");
        assert_eq!(r.uni_number, 1005, "[100] and [001] must match the same UNI");
        assert_eq!(r.bns_number.trim(), "123.345");
        assert_eq!(r.magnetic_type, MagneticType::BlackWhite);
    }

    /// Fe BCC AFM [111]: 2 个 Fe 在 [0,0,0] 和 [0.5,0.5,0.5], 磁矩相反沿 [111]
    /// 预期: R-3m (#166) type-3, BNS=166.101, UNI=1331
    #[test]
    fn test_fe_bcc_afm_111() {
        let lattice = cubic_lattice();
        let n = (3.0_f64).sqrt();
        let positions = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]];
        let types = [26, 26];
        let moments = [
            [1.0 / n, 1.0 / n, 1.0 / n],
            [-1.0 / n, -1.0 / n, -1.0 / n],
        ];
        let r = run_dataset("Fe BCC AFM [111]", &lattice, &positions, &types, Some(&moments));
        assert_eq!(r.spacegroup_number, 229, "non-mag: Im-3m");
        assert!(r.uni_number > 0, "AFM [111] must match a DB entry");
        assert_eq!(r.magnetic_type, MagneticType::BlackWhite);
        assert_eq!(r.bns_number.trim(), "166.101");
        assert_eq!(r.uni_number, 1331);
    }

    /// FCC FM [001]: 4 个原子, 全部磁矩沿 [001]
    #[test]
    fn test_fcc_fm_001() {
        let lattice = cubic_lattice();
        let positions = [
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.0],
            [0.5, 0.0, 0.5],
            [0.0, 0.5, 0.5],
        ];
        let types = [26, 26, 26, 26];
        let moments = [
            [0.0, 0.0, 1.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, 1.0],
        ];
        let r = run_dataset("FCC FM [001]", &lattice, &positions, &types, Some(&moments));
        assert_eq!(r.spacegroup_number, 225, "non-mag: Fm-3m");
        assert_eq!(r.magnetic_type, MagneticType::BlackWhite);
        // FM [001] → tetragonal P4/mmm Type-3, UNI=1005
        assert_eq!(r.uni_number, 1005);
        assert_eq!(r.bns_number.trim(), "123.345");
    }

    /// FCC FM [111]: 4 个原子, 全部磁矩沿 [111]
    /// 预期: R-3m type-3, BNS=166.101, UNI=1331
    #[test]
    fn test_fcc_fm_111() {
        let lattice = cubic_lattice();
        let n = (3.0_f64).sqrt();
        let positions = [
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.0],
            [0.5, 0.0, 0.5],
            [0.0, 0.5, 0.5],
        ];
        let types = [26, 26, 26, 26];
        let m = [1.0 / n, 1.0 / n, 1.0 / n];
        let moments = [m, m, m, m];
        let r = run_dataset("FCC FM [111]", &lattice, &positions, &types, Some(&moments));
        assert_eq!(r.spacegroup_number, 225, "non-mag: Fm-3m");
        assert_eq!(r.magnetic_type, MagneticType::BlackWhite);
        assert!(r.uni_number > 0, "FCC FM [111] must match DB entry");
        assert_eq!(r.bns_number.trim(), "166.101");
        assert_eq!(r.uni_number, 1331);
    }
}
