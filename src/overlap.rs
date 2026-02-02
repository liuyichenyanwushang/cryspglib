// overlap.rs

use crate::cell::Cell;
use crate::mathfunc::{
    Mat3, Vec3, mat_copy_matrix_d3, mat_multiply_matrix_vector_d3, mat_multiply_matrix_vector_id3,
    mat_nint, mat_norm_squared_d3,
};
use std::cmp::Ordering;

/// 用于排序的辅助结构体
#[derive(Clone, Copy)]
struct ValueWithIndex {
    value: f64,
    type_: i32,
    index: usize,
}

impl Default for ValueWithIndex {
    fn default() -> Self {
        Self {
            value: 0.0,
            type_: 0,
            index: 0,
        }
    }
}

/// 重叠检查器结构体
/// 对应 C 代码中的 OverlapChecker
pub struct OverlapChecker {
    pub size: usize,
    pub lattice: Mat3,
    pub types_sorted: Vec<i32>,
    pub pos_sorted: Vec<Vec3>,     // 已排序的原子位置
    pub periodic_axes: [usize; 2], // 用于层状结构
    
    // 缓存区，避免重复分配内存 (对应 C 中的 blob 和 argsort_work)
    perm_temp: Vec<usize>,
    distance_temp: Vec<f64>,
    pos_temp_1: Vec<Vec3>,
    pos_temp_2: Vec<Vec3>,
    argsort_work: Vec<ValueWithIndex>, // 优化：复用排序工作区
}

impl OverlapChecker {
    /// 初始化 OverlapChecker
    /// 对应 C: ovl_overlap_checker_init
    pub fn new(cell: &Cell) -> Option<Self> {
        let size = cell.size;
        // 预分配所有 Vec，避免后续 realloc
        let mut checker = OverlapChecker {
            size,
            lattice: [[0.0; 3]; 3],
            types_sorted: vec![0; size],
            pos_sorted: vec![[0.0; 3]; size],
            periodic_axes: [0, 0],
            perm_temp: vec![0; size],
            distance_temp: vec![0.0; size],
            pos_temp_1: vec![[0.0; 3]; size],
            pos_temp_2: vec![[0.0; 3]; size],
            argsort_work: vec![ValueWithIndex::default(); size],
        };

        mat_copy_matrix_d3(&mut checker.lattice, &cell.lattice);

        // 获取排序排列
        // 注意：这里直接使用 checker 内部的 buffer
        if !argsort_by_lattice_point_distance(
            &mut checker.perm_temp,
            &cell.lattice,
            &cell.position,
            Some(&cell.types),
            &mut checker.distance_temp,
            &mut checker.argsort_work,
        ) {
            return None;
        }

        // 应用排序
        permute_vec3(&mut checker.pos_sorted, &cell.position, &checker.perm_temp);
        permute_i32(&mut checker.types_sorted, &cell.types, &checker.perm_temp);

        // 设置周期性轴 (用于层状结构)
        let mut lattice_rank = 0;
        for i in 0..3 {
            if i as i32 != cell.aperiodic_axis {
                if lattice_rank < 2 {
                    checker.periodic_axes[lattice_rank] = i;
                }
                lattice_rank += 1;
            }
        }

        Some(checker)
    }

    /// 检查完全重叠
    /// 返回: -1 (Error), 0 (False), 1 (True)
    /// 对应 C: ovl_check_total_overlap
    pub fn check_total_overlap(
        &mut self,
        test_trans: &Vec3,
        rot: &[[i32; 3]; 3],
        symprec: f64,
        is_identity: bool,
    ) -> i32 {
        // 快速检查
        if !self.check_possible_overlap(test_trans, rot, symprec) {
            return 0;
        }

        // 计算旋转和平移后的位置
        for i in 0..self.size {
            if is_identity {
                self.pos_temp_1[i] = self.pos_sorted[i];
            } else {
                self.pos_temp_1[i] = mat_multiply_matrix_vector_id3(rot, &self.pos_sorted[i]);
            }

            for k in 0..3 {
                self.pos_temp_1[i][k] += test_trans[k];
            }
        }

        // 对变换后的位置进行排序
        if !argsort_by_lattice_point_distance(
            &mut self.perm_temp,
            &self.lattice,
            &self.pos_temp_1,
            Some(&self.types_sorted),
            &mut self.distance_temp,
            &mut self.argsort_work,
        ) {
            return -1;
        }

        permute_vec3(&mut self.pos_temp_2, &self.pos_temp_1, &self.perm_temp);

        // 检查排序后的重叠
        check_total_overlap_for_sorted(
            &self.lattice,
            &self.pos_sorted,
            &self.pos_temp_2,
            &self.types_sorted,
            &self.types_sorted,
            self.size,
            symprec,
        )
    }

    /// 检查层状结构的完全重叠
    /// 对应 C: ovl_check_layer_total_overlap
    pub fn check_layer_total_overlap(
        &mut self,
        test_trans: &Vec3,
        rot: &[[i32; 3]; 3],
        symprec: f64,
        is_identity: bool,
    ) -> i32 {
        if !self.check_possible_overlap(test_trans, rot, symprec) {
            return 0;
        }

        for i in 0..self.size {
            if is_identity {
                self.pos_temp_1[i] = self.pos_sorted[i];
            } else {
                self.pos_temp_1[i] = mat_multiply_matrix_vector_id3(rot, &self.pos_sorted[i]);
            }
            for k in 0..3 {
                self.pos_temp_1[i][k] += test_trans[k];
            }
        }

        if !argsort_by_lattice_point_distance(
            &mut self.perm_temp,
            &self.lattice,
            &self.pos_temp_1,
            Some(&self.types_sorted),
            &mut self.distance_temp,
            &mut self.argsort_work,
        ) {
            return -1;
        }

        permute_vec3(&mut self.pos_temp_2, &self.pos_temp_1, &self.perm_temp);

        check_layer_total_overlap_for_sorted(
            &self.lattice,
            &self.pos_sorted,
            &self.pos_temp_2,
            &self.types_sorted,
            &self.types_sorted,
            self.size,
            &self.periodic_axes,
            symprec,
        )
    }

    /// 快速预检查
    /// 对应 C: check_possible_overlap
    fn check_possible_overlap(&self, test_trans: &Vec3, rot: &[[i32; 3]; 3], symprec: f64) -> bool {
        let max_search_num = 3;
        let search_num = if self.size <= max_search_num {
            self.size
        } else {
            max_search_num
        };

        for i_test in 0..search_num {
            let type_rot = self.types_sorted[i_test];
            let mut pos_rot = mat_multiply_matrix_vector_id3(rot, &self.pos_sorted[i_test]);
            for k in 0..3 {
                pos_rot[k] += test_trans[k];
            }

            let mut is_found = false;
            // 暴力搜索，因为只检查前几个原子，开销可控
            for i in 0..self.size {
                if has_overlap_with_same_type(
                    &pos_rot,
                    &self.pos_sorted[i],
                    type_rot,
                    self.types_sorted[i],
                    &self.lattice,
                    symprec,
                ) {
                    is_found = true;
                    break;
                }
            }

            if !is_found {
                return false;
            }
        }
        true
    }
}

// --- Helper Functions ---

#[inline]
fn cartesian_norm(lat: &Mat3, v: &Vec3) -> f64 {
    let mut temp = [0.0; 3];
    for i in 0..3 {
        temp[i] = lat[i][0] * v[0] + lat[i][1] * v[1] + lat[i][2] * v[2];
    }
    (temp[0] * temp[0] + temp[1] * temp[1] + temp[2] * temp[2]).sqrt()
}

#[inline]
fn has_overlap(a: &Vec3, b: &Vec3, lattice: &Mat3, symprec: f64) -> bool {
    let mut v_diff = [0.0; 3];
    for i in 0..3 {
        v_diff[i] = a[i] - b[i];
        v_diff[i] -= mat_nint(v_diff[i]) as f64;
    }
    cartesian_norm(lattice, &v_diff) <= symprec
}

#[inline]
fn has_overlap_with_same_type(
    a: &Vec3,
    b: &Vec3,
    type_a: i32,
    type_b: i32,
    lattice: &Mat3,
    symprec: f64,
) -> bool {
    if type_a == type_b {
        has_overlap(a, b, lattice, symprec)
    } else {
        false
    }
}

#[inline]
fn layer_has_overlap(
    a: &Vec3,
    b: &Vec3,
    lattice: &Mat3,
    periodic_axes: &[usize; 2],
    symprec: f64,
) -> bool {
    let mut v_diff = [0.0; 3];
    for i in 0..3 {
        v_diff[i] = a[i] - b[i];
    }
    // 仅在周期性方向上应用最小镜像约定
    v_diff[periodic_axes[0]] -= mat_nint(v_diff[periodic_axes[0]]) as f64;
    v_diff[periodic_axes[1]] -= mat_nint(v_diff[periodic_axes[1]]) as f64;

    cartesian_norm(lattice, &v_diff) <= symprec
}

#[inline]
fn layer_has_overlap_with_same_type(
    a: &Vec3,
    b: &Vec3,
    type_a: i32,
    type_b: i32,
    lattice: &Mat3,
    periodic_axes: &[usize; 2],
    symprec: f64,
) -> bool {
    if type_a == type_b {
        layer_has_overlap(a, b, lattice, periodic_axes, symprec)
    } else {
        false
    }
}

// --- Sorting Logic ---

/// 根据到最近格点的距离对原子进行排序
/// 对应 C: argsort_by_lattice_point_distance
fn argsort_by_lattice_point_distance(
    perm: &mut [usize],
    lattice: &Mat3,
    positions: &[Vec3],
    types: Option<&[i32]>,
    distance_temp: &mut [f64],
    work: &mut [ValueWithIndex], // 传入预分配的工作区
) -> bool {
    let size = positions.len();
    let mut diff = [0.0; 3];

    for i in 0..size {
        for k in 0..3 {
            let x = positions[i][k];
            diff[k] = x - mat_nint(x) as f64;
        }
        let diff_cart = mat_multiply_matrix_vector_d3(lattice, &diff);
        distance_temp[i] = mat_norm_squared_d3(&diff_cart);
    }

    // 填充工作区
    for i in 0..size {
        work[i].value = distance_temp[i];
        work[i].type_ = if let Some(t) = types { t[i] } else { 0 };
        work[i].index = i;
    }

    // 排序：先按类型降序，再按距离升序
    // 对应 C: ValueWithIndex_comparator
    work.sort_by(|a, b| {
        let type_cmp = b.type_.cmp(&a.type_); // Descending
        if type_cmp != Ordering::Equal {
            type_cmp
        } else {
            a.value.partial_cmp(&b.value).unwrap_or(Ordering::Equal) // Ascending
        }
    });

    for i in 0..size {
        perm[i] = work[i].index;
    }

    true
}

fn permute_vec3(data_out: &mut [Vec3], data_in: &[Vec3], perm: &[usize]) {
    for (i, &idx) in perm.iter().enumerate() {
        data_out[i] = data_in[idx];
    }
}

fn permute_i32(data_out: &mut [i32], data_in: &[i32], perm: &[usize]) {
    for (i, &idx) in perm.iter().enumerate() {
        data_out[i] = data_in[idx];
    }
}

/// 检查两个已排序的原子列表是否重叠
/// 对应 C: check_total_overlap_for_sorted
fn check_total_overlap_for_sorted(
    lattice: &Mat3,
    pos_original: &[Vec3],
    pos_rotated: &[Vec3],
    types_original: &[i32],
    types_rotated: &[i32],
    num_pos: usize,
    symprec: f64,
) -> i32 {
    let mut found = vec![false; num_pos];
    let mut search_start = 0;

    for i_orig in 0..num_pos {
        // 跳过开头已经匹配过的原子
        while search_start < num_pos && found[search_start] {
            search_start += 1;
        }

        let mut matched = false;
        for i_rot in search_start..num_pos {
            if found[i_rot] {
                continue;
            }

            if has_overlap_with_same_type(
                &pos_original[i_orig],
                &pos_rotated[i_rot],
                types_original[i_orig],
                types_rotated[i_rot],
                lattice,
                symprec,
            ) {
                found[i_rot] = true;
                matched = true;
                break;
            }
        }

        if !matched {
            return 0;
        }
    }
    1
}

/// 检查两个已排序的原子列表是否重叠 (层状结构)
/// 对应 C: check_layer_total_overlap_for_sorted
fn check_layer_total_overlap_for_sorted(
    lattice: &Mat3,
    pos_original: &[Vec3],
    pos_rotated: &[Vec3],
    types_original: &[i32],
    types_rotated: &[i32],
    num_pos: usize,
    periodic_axes: &[usize; 2],
    symprec: f64,
) -> i32 {
    let mut found = vec![false; num_pos];
    let mut search_start = 0;

    for i_orig in 0..num_pos {
        while search_start < num_pos && found[search_start] {
            search_start += 1;
        }

        let mut matched = false;
        for i_rot in search_start..num_pos {
            if found[i_rot] {
                continue;
            }

            if layer_has_overlap_with_same_type(
                &pos_original[i_orig],
                &pos_rotated[i_rot],
                types_original[i_orig],
                types_rotated[i_rot],
                lattice,
                periodic_axes,
                symprec,
            ) {
                found[i_rot] = true;
                matched = true;
                break;
            }
        }

        if !matched {
            return 0;
        }
    }
    1
}
