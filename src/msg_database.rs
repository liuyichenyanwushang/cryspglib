//! 磁性空间群数据库。
//!
//! 包含磁性空间群类型、对称操作、Hall 映射等数据以及查询函数。
//! 数据部分由 `msg_database_gen.rs` 提供（自动转换自 C）。
//! 逻辑函数手动移植自 `msg_database.txt`。

include!("msg_database_gen.rs");

use crate::debug;
use crate::mathfunc::Mat3I;
use crate::symmetry::{MagneticSymmetry, Symmetry};

/// 获取磁性空间群类型。
pub fn msgdb_get_magnetic_spacegroup_type(uni_number: i32) -> &'static MagneticSpacegroupType {
    if uni_number > 0 && (uni_number as usize) < MAGNETIC_SPACEGROUP_TYPES.len() {
        &MAGNETIC_SPACEGROUP_TYPES[uni_number as usize]
    } else {
        &MAGNETIC_SPACEGROUP_TYPES[0]
    }
}

/// 返回给定 Hall 编号的最小和最大 UNI 编号。
/// `uni_number_range[0]` = 最小, `[1]` = 最大。
pub fn msgdb_get_uni_candidates(hall_number: i32) -> Option<[i32; 2]> {
    let hall = hall_number as usize;
    if hall >= MAGNETIC_SPACEGROUP_HALL_MAPPING.len() {
        return None;
    }
    Some(MAGNETIC_SPACEGROUP_HALL_MAPPING[hall])
}

/// 获取磁性对称操作。
pub fn msgdb_get_spacegroup_operations(
    uni_number: i32,
    hall_number: i32,
) -> Option<MagneticSymmetry> {
    let hall_number_offset = get_hall_number_offset(uni_number, hall_number)?;
    let order = MAGNETIC_SPACEGROUP_OPERATION_INDEX[uni_number as usize][hall_number_offset as usize][0];
    let start = MAGNETIC_SPACEGROUP_OPERATION_INDEX[uni_number as usize][hall_number_offset as usize][1];

    debug::debug_print(format_args!(
        "Load operations with UNI={}, hall={}: hall_number_offset={} order={} starts={}\n",
        uni_number, hall_number, hall_number_offset, order, start
    ));

    let mut sym = MagneticSymmetry::new(order as usize);
    for i in 0..order {
        let idx = (start + i) as usize;
        let (rot, trans, timerev) = msgdb_get_magnetic_operation(idx);
        sym.rot[i as usize] = rot;
        sym.trans[i as usize] = trans;
        sym.timerev[i as usize] = timerev;
    }

    Some(sym)
}

/// 获取标准变换。
pub fn msgdb_get_std_transformations(
    uni_number: i32,
    hall_number: i32,
) -> Option<Symmetry> {
    let hall_number_offset = get_hall_number_offset(uni_number, hall_number)?;
    let mut sym = Symmetry::new(7);
    // Identity transformation as first element
    sym.rot[0] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    sym.trans[0] = [0.0; 3];

    for i in 0..7 {
        let enc = ALTERNATIVE_TRANSFORMATIONS[uni_number as usize][hall_number_offset as usize][i];
        if enc == 0 {
            sym.size = i + 1;
            break;
        }
        let (tmat, origin_shift) = crate::spg_database::spgdb_decode_symmetry(enc);
        sym.rot[i + 1] = tmat;
        sym.trans[i + 1] = origin_shift;
    }

    Some(sym)
}

/// 解码磁性对称操作。
fn msgdb_get_magnetic_operation(op_number: usize) -> (Mat3I, [f64; 3], i32) {
    let enc = MAGNETIC_SYMMETRY_OPERATIONS[op_number];
    // timerev=1 for anti operation, timerev=0 for ordinary operation
    // 34012224 = 3^9 * 12^3
    let timerev = enc / 34012224;
    let (rot, trans) = crate::spg_database::spgdb_decode_symmetry(enc % 34012224);
    (rot, trans, timerev)
}

/// 获取 Hall 编号偏移量。
fn get_hall_number_offset(uni_number: i32, hall_number: i32) -> Option<i32> {
    if uni_number <= 0 || uni_number > 1651 {
        return None;
    }

    let num_halls = MAGNETIC_SPACEGROUP_UNI_MAPPING[uni_number as usize][0];
    let first_hall_number = MAGNETIC_SPACEGROUP_UNI_MAPPING[uni_number as usize][1];

    let hall_number_offset = if hall_number > 0 && hall_number <= 530 {
        hall_number - first_hall_number
    } else if hall_number == 0 {
        // Use the smallest Hall number as default
        0
    } else {
        return None;
    };

    // check if hall_number_offset is in range
    if hall_number_offset < 0 || hall_number_offset >= num_halls {
        return None;
    }

    Some(hall_number_offset)
}
