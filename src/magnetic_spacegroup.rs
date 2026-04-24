//! 磁性空间群识别。
//!
//! 使用磁性对称操作数据库识别磁性空间群类型。
//! 参考: Litvin, "Magnetic Group Tables" 2013

use crate::debug;
use crate::mathfunc::{
    mat_cast_matrix_3i_to_3d, mat_check_identity_matrix_i3, mat_copy_matrix_d3,
    mat_copy_vector_d3, mat_inverse_matrix_d3, mat_multiply_matrix_d3,
    mat_multiply_matrix_vector_d3, mat_multiply_matrix_vector_id3, mat_nint, Mat3, Mat3I, Vec3,
};
use crate::msg_database::{
    msgdb_get_magnetic_spacegroup_type, msgdb_get_spacegroup_operations,
    msgdb_get_std_transformations, msgdb_get_uni_candidates,
};
use crate::primitive::prm_get_primitive_symmetry;
use crate::refinement::ref_find_similar_bravais_lattice;
use crate::spacegroup::{spa_copy_spacegroup, spa_search_spacegroup_with_symmetry, Spacegroup};
use crate::spg_database::spgdb_get_spacegroup_type;
use crate::symmetry::{MagneticSymmetry, Symmetry};

const MAX_DENOMINATOR: f64 = 100.0;

/// 磁性空间群识别结果的中间数据结构。
pub struct MagneticDataset {
    pub uni_number: i32,
    pub msg_type: i32,
    pub hall_number: i32,
    pub transformation_matrix: Mat3,
    pub origin_shift: Vec3,
    pub std_rotation_matrix: Mat3,
}

/// 识别磁性空间群类型。
///
/// 给定晶格和磁性对称操作，返回识别出的磁性数据集。
/// 如果标准路径失败（磁对称性为非完整空间群），在 `parent_hall_number` 提供时
/// 使用非磁母空间群的 Hall 编号作为 fallback。
pub fn msg_identify_magnetic_space_group_type(
    lattice: &Mat3,
    magnetic_symmetry: &MagneticSymmetry,
    symprec: f64,
) -> Option<MagneticDataset> {
    msg_identify_with_parent_hall(lattice, magnetic_symmetry, None, symprec)
}

/// 与 [`msg_identify_magnetic_space_group_type`] 相同，但可指定非磁母空间群的
/// Hall 编号作为 fallback，适用于磁对称性显著低于母空间群的场景。
pub fn msg_identify_with_parent_hall(
    lattice: &Mat3,
    magnetic_symmetry: &MagneticSymmetry,
    parent_hall_number: Option<i32>,
    symprec: f64,
) -> Option<MagneticDataset> {
    // 标准路径: 从磁对称性中提取 FSG/XSG 并搜索空间群
    let (ref_sg, changed_symmetry, mut tmat, mut shift, msgtype_num) =
        match get_reference_space_group(magnetic_symmetry, symprec) {
            Some(result) => result,
            None => {
                // 标准路径失败 → 尝试 fallback（用母空间群的 Hall 编号）
                let hall_number = parent_hall_number?;
                build_fallback_reference(lattice, magnetic_symmetry, hall_number, symprec)?
            }
        };

    let hall_number = ref_sg.hall_number;
    debug::debug_print(format_args!("Search UNI number over between...\n"));

    let range = msgdb_get_uni_candidates(hall_number)?;
    let min_uni = range[0];
    let max_uni = range[1];

    let mut best_uni = 0;
    let mut best_msg_type = 0;

    for uni_number in min_uni..=max_uni {
        // Check type and order
        let msgtype_db = msgdb_get_magnetic_spacegroup_type(uni_number);
        if msgtype_db.type_ != msgtype_num {
            continue;
        }
        let msg_uni = msgdb_get_spacegroup_operations(uni_number, hall_number)?;
        if msg_uni.size != changed_symmetry.size {
            continue;
        }

        let transformations = msgdb_get_std_transformations(uni_number, hall_number)?;

        let mut same = false;
        for trans_idx in 0..transformations.size {
            let tmat_cor = transformations.rot[trans_idx];
            let shift_cor = transformations.trans[trans_idx];

            let tmat_cor_d = mat_cast_matrix_3i_to_3d(&tmat_cor);

            let symmetry_cor = get_distinct_changed_magnetic_symmetry(
                &tmat_cor_d, &shift_cor, &changed_symmetry,
            )?;

            if is_equal(&symmetry_cor, &msg_uni, symprec) {
                same = true;
                // Update transformation: tmat = tmat_cor @ tmat
                tmat = mat_multiply_matrix_d3(&tmat_cor_d, &tmat);
                let shift_tmp = mat_multiply_matrix_vector_d3(&tmat_cor_d, &shift);
                for s in 0..3 {
                    shift[s] = shift_tmp[s] + shift_cor[s];
                }
                break;
            }
        }

        if same {
            best_uni = uni_number;
            best_msg_type = msgtype_db.type_;
            break;
        }
    }

    if best_uni == 0 {
        return None;
    }

    let _msgtype = msgdb_get_magnetic_spacegroup_type(best_uni);
    let mut ret = MagneticDataset {
        uni_number: best_uni,
        msg_type: best_msg_type,
        hall_number,
        transformation_matrix: [[0.0; 3]; 3],
        origin_shift: [0.0; 3],
        std_rotation_matrix: [[0.0; 3]; 3],
    };

    get_rigid_rotation(&mut ret.std_rotation_matrix, lattice, &tmat, &ref_sg);
    mat_copy_matrix_d3(&mut ret.transformation_matrix, &tmat);
    mat_copy_vector_d3(&mut ret.origin_shift, &shift);

    Some(ret)
}

/// 获取参考空间群和变换后的磁性对称操作。
fn get_reference_space_group(
    magnetic_symmetry: &MagneticSymmetry,
    symprec: f64,
) -> Option<(Spacegroup, MagneticSymmetry, Mat3, Vec3, i32)> {
    // 1. Get family space group (FSG) and its symmetry
    let (mut fsg, sym_fsg) =
        get_family_space_group_with_magnetic_symmetry(magnetic_symmetry, symprec)?;

    // 2. Get maximal subspace group (XSG) and its symmetry
    let (mut xsg, sym_xsg) =
        get_maximal_subspace_group_with_magnetic_symmetry(magnetic_symmetry, symprec)?;

    // 3. Determine type of MSG
    let msgtype_num =
        get_magnetic_space_group_type(magnetic_symmetry, sym_fsg.size, sym_xsg.size)?;

    // 4. Choose reference setting
    // For type-IV, use setting from Hall symbol of XSG.
    // For other types, use setting from Hall symbol of FSG.
    let ref_sg = if msgtype_num == 4 { &mut xsg } else { &mut fsg };

    // 5. Compute tmat and shift from reference spacegroup
    // tmat = bravais_lattice^{-1}, shift = origin_shift
    let tmat = mat_inverse_matrix_d3(&ref_sg.bravais_lattice, 0.0)?;
    let shift = ref_sg.origin_shift;

    // 6. Compute changed symmetry by transforming magnetic_symmetry to reference setting
    let changed_symmetry =
        get_distinct_changed_magnetic_symmetry(&tmat, &shift, magnetic_symmetry)?;

    // 7. Copy ref_sg for return
    let mut ref_sg_copy = Spacegroup::new();
    spa_copy_spacegroup(&mut ref_sg_copy, ref_sg);

    Some((ref_sg_copy, changed_symmetry, tmat, shift, msgtype_num))
}

/// Fallback when `get_reference_space_group` fails.
/// 使用提供的非磁 Hall 编号构建参考空间群，绕过空间群搜索。
fn build_fallback_reference(
    lattice: &Mat3,
    magnetic_symmetry: &MagneticSymmetry,
    parent_hall_number: i32,
    symprec: f64,
) -> Option<(Spacegroup, MagneticSymmetry, Mat3, Vec3, i32)> {
    // 1. 提取 FSG/XSG symmetry (只取操作，不搜索空间群)
    let sym_fsg = extract_symmetry(magnetic_symmetry, true, symprec)?;
    let sym_xsg = extract_symmetry(magnetic_symmetry, false, symprec)?;

    // 2. 确定磁性类型
    let msgtype_num =
        get_magnetic_space_group_type(magnetic_symmetry, sym_fsg.size, sym_xsg.size)?;

    // 3. 用非磁 Hall 编号构建参考 Spacegroup
    let spg_type = spgdb_get_spacegroup_type(parent_hall_number);
    let mut ref_sg = Spacegroup::new();
    ref_sg.hall_number = parent_hall_number;
    ref_sg.number = spg_type.number;
    ref_sg.pointgroup_number = spg_type.pointgroup_number;
    ref_sg.schoenflies = spg_type.schoenflies;
    ref_sg.hall_symbol = spg_type.hall_symbol;
    ref_sg.international = spg_type.international;
    ref_sg.international_long = spg_type.international_full;
    ref_sg.international_short = spg_type.international_short;
    ref_sg.choice = spg_type.choice;
    mat_copy_matrix_d3(&mut ref_sg.bravais_lattice, lattice);
    // origin_shift 保持为 [0;3]（默认值）

    // 4. 计算 changed_symmetry: 将磁对称性变换到参考设置
    let tmat = mat_inverse_matrix_d3(&ref_sg.bravais_lattice, 0.0)?;
    let shift = ref_sg.origin_shift;
    let changed_symmetry =
        get_distinct_changed_magnetic_symmetry(&tmat, &shift, magnetic_symmetry)?;

    // 5. 复制 ref_sg 用于返回
    let mut ref_sg_copy = Spacegroup::new();
    spa_copy_spacegroup(&mut ref_sg_copy, &ref_sg);

    Some((ref_sg_copy, changed_symmetry, tmat, shift, msgtype_num))
}

/// 从磁性对称操作中提取普通对称操作（不搜索空间群）。
/// `ignore_time_reversal=true` → FSG (所有操作)
/// `ignore_time_reversal=false` → XSG (仅 timerev=0)
fn extract_symmetry(
    magnetic_symmetry: &MagneticSymmetry,
    ignore_time_reversal: bool,
    symprec: f64,
) -> Option<Symmetry> {
    let identity: Mat3I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    let num_sym_msg = magnetic_symmetry.size;

    // Check if MSG is type-II
    let is_type2 = magnetic_symmetry
        .rot
        .iter()
        .zip(magnetic_symmetry.trans.iter())
        .zip(magnetic_symmetry.timerev.iter())
        .any(|((rot, trans), &timerev)| {
            mat_check_identity_matrix_i3(&identity, rot)
                && trans[0].abs() < symprec
                && trans[1].abs() < symprec
                && trans[2].abs() < symprec
                && timerev == 1
        });

    let mut sym = Symmetry::new(num_sym_msg);
    let mut num_sym = 0;
    for i in 0..num_sym_msg {
        if !ignore_time_reversal && magnetic_symmetry.timerev[i] == 1 {
            continue;
        }
        if is_type2 && magnetic_symmetry.timerev[i] == 1 {
            continue;
        }
        sym.rot[num_sym] = magnetic_symmetry.rot[i];
        sym.trans[num_sym] = magnetic_symmetry.trans[i];
        num_sym += 1;
    }
    sym.size = num_sym;

    if num_sym == 0 {
        None
    } else {
        Some(sym)
    }
}

/// Get family space group (FSG) and its symmetry.
fn get_family_space_group_with_magnetic_symmetry(
    magnetic_symmetry: &MagneticSymmetry,
    symprec: f64,
) -> Option<(Spacegroup, Symmetry)> {
    get_space_group_with_magnetic_symmetry(magnetic_symmetry, true, symprec)
}

/// Get maximal subspace group (XSG) and its symmetry.
fn get_maximal_subspace_group_with_magnetic_symmetry(
    magnetic_symmetry: &MagneticSymmetry,
    symprec: f64,
) -> Option<(Spacegroup, Symmetry)> {
    get_space_group_with_magnetic_symmetry(magnetic_symmetry, false, symprec)
}

/// Get space group from magnetic symmetry.
///
/// `ignore_time_reversal=true` → FSG (family space group, all ops regardless of timerev)
/// `ignore_time_reversal=false` → XSG (maximal subspace group, only ordinary ops)
///
/// Returns (Spacegroup, Symmetry) pair.
fn get_space_group_with_magnetic_symmetry(
    magnetic_symmetry: &MagneticSymmetry,
    ignore_time_reversal: bool,
    symprec: f64,
) -> Option<(Spacegroup, Symmetry)> {
    let identity: Mat3I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    let unit_lat: Mat3 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

    let num_sym_msg = magnetic_symmetry.size;

    // Check if MSG is type-II (has pure time-reversal operation (I, 0)1')
    let is_type2 = magnetic_symmetry
        .rot
        .iter()
        .zip(magnetic_symmetry.trans.iter())
        .zip(magnetic_symmetry.timerev.iter())
        .any(|((rot, trans), &timerev)| {
            mat_check_identity_matrix_i3(&identity, rot)
                && trans[0].abs() < symprec
                && trans[1].abs() < symprec
                && trans[2].abs() < symprec
                && timerev == 1
        });

    // Extract operations
    let mut sym = Symmetry::new(num_sym_msg);
    let mut num_sym = 0;
    for i in 0..num_sym_msg {
        if !ignore_time_reversal && magnetic_symmetry.timerev[i] == 1 {
            continue;
        }
        // For type-II MSG, deduplicate: skip timerev=1 copies
        if is_type2 && magnetic_symmetry.timerev[i] == 1 {
            continue;
        }
        sym.rot[num_sym] = magnetic_symmetry.rot[i];
        sym.trans[num_sym] = magnetic_symmetry.trans[i];
        num_sym += 1;
    }
    sym.size = num_sym;

    if num_sym == 0 {
        return None;
    }

    // Get primitive symmetry: (a, b, c) = (a_prim, b_prim, c_prim) @ tmat
    let mut tmat = [[0.0; 3]; 3];
    let prim_sym = match prm_get_primitive_symmetry(&mut tmat, &sym, symprec) {
        Some(ps) => ps,
        None => {
            debug::debug_print(format_args!(
                "get_space_group_with_magnetic_symmetry: prm_get_primitive_symmetry failed (n_sym={})\n",
                num_sym
            ));
            return None;
        }
    };

    // Search space group with primitive symmetry and unit lattice
    let mut spacegroup = match spa_search_spacegroup_with_symmetry(&prim_sym, &unit_lat, symprec) {
        Some(sg) => sg,
        None => {
            // Fallback: try direct Hall matching by point group
            debug::debug_print(format_args!(
                "get_space_group_with_magnetic_symmetry: trying fallback (n_prim={}, n_sym={})\n",
                prim_sym.size, num_sym
            ));
            find_spacegroup_by_symmetry(&sym, &unit_lat, symprec)?
        }
    };

    // Refine bravais lattice and origin_shift
    ref_find_similar_bravais_lattice(&mut spacegroup, symprec);

    // Change basis from primitive to original:
    // x = (tmat, 0)^-1 x_prim
    // => x_std = (P^-1, p) (tmat, 0) x = ( P^-1 @ tmat, p) x
    //    (a_std, b_std, c_std) = (a, b, c) @ tmat^-1 @ P
    let inv_tmat = mat_inverse_matrix_d3(&tmat, 0.0)?;
    spacegroup.bravais_lattice = mat_multiply_matrix_d3(&inv_tmat, &spacegroup.bravais_lattice);

    Some((spacegroup, sym))
}

/// Determine MSG type. Returns None if failed.
fn get_magnetic_space_group_type(
    magnetic_symmetry: &MagneticSymmetry,
    num_sym_fsg: usize,
    num_sym_xsg: usize,
) -> Option<i32> {
    let identity: Mat3I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

    if num_sym_fsg == num_sym_xsg {
        let num_sym_msg = magnetic_symmetry.size;
        if num_sym_msg == num_sym_fsg {
            // Type-I: all operations are ordinary
            Some(1)
        } else if num_sym_msg == 2 * num_sym_fsg {
            // Type-II: has pure time-reversal operation
            Some(2)
        } else {
            None
        }
    } else if num_sym_fsg == 2 * num_sym_xsg {
        // Check if anti-operation is translation (type-IV) or rotation (type-III)
        let has_anti_translation = magnetic_symmetry
            .rot
            .iter()
            .zip(magnetic_symmetry.timerev.iter())
            .any(|(rot, &timerev)| mat_check_identity_matrix_i3(&identity, rot) && timerev == 1);

        if has_anti_translation {
            Some(4) // Type-IV: anti-translation
        } else {
            Some(3) // Type-III: anti-rotation
        }
    } else {
        None
    }
}

/// Get coset representative of XSG in MSG.
/// For type-III and type-IV MSGs. Returns identity for type-I/II.
fn get_representative(magnetic_symmetry: &MagneticSymmetry) -> Option<MagneticSymmetry> {
    let identity: Mat3I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

    let mut representative = MagneticSymmetry::new(2);

    // Set the first representative as identity
    representative.rot[0] = identity;
    representative.trans[0] = [0.0; 3];
    representative.timerev[0] = 0;

    // Type-IV: anti-translation with identity rotation
    for i in 0..magnetic_symmetry.size {
        if mat_check_identity_matrix_i3(&identity, &magnetic_symmetry.rot[i])
            && magnetic_symmetry.timerev[i] == 1
        {
            representative.rot[1] = magnetic_symmetry.rot[i];
            representative.trans[1] = magnetic_symmetry.trans[i];
            representative.timerev[1] = 1;
            representative.size = 2;
            return Some(representative);
        }
    }

    // Type-III: anti-rotation
    for i in 0..magnetic_symmetry.size {
        if magnetic_symmetry.timerev[i] != 1 {
            continue;
        }
        representative.rot[1] = magnetic_symmetry.rot[i];
        representative.trans[1] = magnetic_symmetry.trans[i];
        representative.timerev[1] = 1;
        representative.size = 2;
        return Some(representative);
    }

    // Type-I or II: no anti-operations
    representative.size = 1;
    Some(representative)
}

/// Apply transformation (tmat, shift) to magnetic symmetry, deduplicating.
///
/// R' = T^-1 * R * T
/// t' = T^-1 * (t + (R - I) * shift)
fn get_distinct_changed_magnetic_symmetry(
    tmat: &Mat3,
    shift: &Vec3,
    sym_msg: &MagneticSymmetry,
) -> Option<MagneticSymmetry> {
    let inv_tmat = mat_inverse_matrix_d3(tmat, 0.0)?;

    let mut changed = MagneticSymmetry::new(sym_msg.size);
    let mut size = 0usize;

    for i in 0..sym_msg.size {
        // R' = T^-1 * R * T
        let rot_f64 = mat_cast_matrix_3i_to_3d(&sym_msg.rot[i]);
        let tmp = mat_multiply_matrix_d3(&inv_tmat, &rot_f64);
        let r_new = mat_multiply_matrix_d3(&tmp, tmat);

        // Round to integer rotation matrix
        let rot_i = [
            [mat_nint(r_new[0][0]), mat_nint(r_new[0][1]), mat_nint(r_new[0][2])],
            [mat_nint(r_new[1][0]), mat_nint(r_new[1][1]), mat_nint(r_new[1][2])],
            [mat_nint(r_new[2][0]), mat_nint(r_new[2][1]), mat_nint(r_new[2][2])],
        ];

        // t' = T^-1 * (t + (R - I) * shift)
        let mut tmp_rot = sym_msg.rot[i];
        tmp_rot[0][0] -= 1;
        tmp_rot[1][1] -= 1;
        tmp_rot[2][2] -= 1;
        let shift_term = mat_multiply_matrix_vector_id3(&tmp_rot, shift);

        let mut t_new = [0.0; 3];
        for j in 0..3 {
            t_new[j] = sym_msg.trans[i][j] + shift_term[j];
        }
        t_new = mat_multiply_matrix_vector_d3(&inv_tmat, &t_new);

        // Check for uniqueness (same rotation, same translation, same timerev)
        let mut is_dup = false;
        for j in 0..size {
            if mat_check_identity_matrix_i3(&changed.rot[j], &rot_i) {
                let mut diff = [0.0; 3];
                for k in 0..3 {
                    diff[k] = changed.trans[j][k] - t_new[k];
                    diff[k] -= mat_nint(diff[k]) as f64;
                }
                if diff[0].abs() < 1e-5 && diff[1].abs() < 1e-5 && diff[2].abs() < 1e-5 {
                    if changed.timerev[j] == sym_msg.timerev[i] {
                        is_dup = true;
                        break;
                    }
                }
            }
        }

        if !is_dup {
            changed.rot[size] = rot_i;
            changed.trans[size] = t_new;
            changed.timerev[size] = sym_msg.timerev[i];
            size += 1;
        }
    }

    changed.size = size;
    Some(changed)
}

/// 检查两个磁性对称操作集合是否等价。
fn is_equal(
    sym1: &MagneticSymmetry,
    sym2: &MagneticSymmetry,
    symprec: f64,
) -> bool {
    if sym1.size != sym2.size {
        return false;
    }

    let mut found = vec![false; sym2.size];
    for i in 0..sym1.size {
        let mut matched = false;
        for j in 0..sym2.size {
            if found[j] {
                continue;
            }
            if !mat_check_identity_matrix_i3(&sym1.rot[i], &sym2.rot[j]) {
                continue;
            }
            if sym1.timerev[i] != sym2.timerev[j] {
                continue;
            }
            let mut diff = [0.0; 3];
            for k in 0..3 {
                diff[k] = sym1.trans[i][k] - sym2.trans[j][k];
                diff[k] -= mat_nint(diff[k]) as f64;
            }
            if diff[0].abs() < symprec && diff[1].abs() < symprec && diff[2].abs() < symprec {
                found[j] = true;
                matched = true;
                break;
            }
        }
        if !matched {
            return false;
        }
    }
    true
}

/// 计算刚性旋转矩阵。
fn get_rigid_rotation(
    rigid_rot: &mut Mat3,
    lattice: &Mat3,
    tmat: &Mat3,
    ref_sg: &Spacegroup,
) {
    let inv_tmat = mat_inverse_matrix_d3(tmat, 0.0);
    if let Some(inv) = inv_tmat {
        let tmp = mat_multiply_matrix_d3(&ref_sg.bravais_lattice, &inv);
        let inv_lat = mat_inverse_matrix_d3(lattice, 0.0);
        if let Some(inv_l) = inv_lat {
            let result = mat_multiply_matrix_d3(&inv_l, &tmp);
            *rigid_rot = result;
        }
    }
}

