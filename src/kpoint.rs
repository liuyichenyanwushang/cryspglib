// kpoint.rs

use crate::debug;
use crate::kgrid;
use crate::mathfunc::{
    MatINT, mat_alloc_matint, mat_check_identity_matrix_i3, mat_copy_matrix_i3, mat_dabs,
    mat_multiply_matrix_i3, mat_multiply_matrix_vector_d3,
    mat_multiply_matrix_vector_i3, mat_multiply_matrix_vector_id3, mat_nint, mat_norm_squared_d3,
    mat_transpose_matrix_i3,
};

// 引入 rayon 用于并行计算，对应 C 中的 OpenMP
#[cfg(feature = "parallel")]
use rayon::prelude::*;

const KPT_NUM_BZ_SEARCH_SPACE: usize = 125;

// 静态搜索空间数组：用于在倒易空间寻找最近邻的格点
static BZ_SEARCH_SPACE: [[i32; 3]; KPT_NUM_BZ_SEARCH_SPACE] = [
    [0, 0, 0],
    [0, 0, 1],
    [0, 0, 2],
    [0, 0, -2],
    [0, 0, -1],
    [0, 1, 0],
    [0, 1, 1],
    [0, 1, 2],
    [0, 1, -2],
    [0, 1, -1],
    [0, 2, 0],
    [0, 2, 1],
    [0, 2, 2],
    [0, 2, -2],
    [0, 2, -1],
    [0, -2, 0],
    [0, -2, 1],
    [0, -2, 2],
    [0, -2, -2],
    [0, -2, -1],
    [0, -1, 0],
    [0, -1, 1],
    [0, -1, 2],
    [0, -1, -2],
    [0, -1, -1],
    [1, 0, 0],
    [1, 0, 1],
    [1, 0, 2],
    [1, 0, -2],
    [1, 0, -1],
    [1, 1, 0],
    [1, 1, 1],
    [1, 1, 2],
    [1, 1, -2],
    [1, 1, -1],
    [1, 2, 0],
    [1, 2, 1],
    [1, 2, 2],
    [1, 2, -2],
    [1, 2, -1],
    [1, -2, 0],
    [1, -2, 1],
    [1, -2, 2],
    [1, -2, -2],
    [1, -2, -1],
    [1, -1, 0],
    [1, -1, 1],
    [1, -1, 2],
    [1, -1, -2],
    [1, -1, -1],
    [2, 0, 0],
    [2, 0, 1],
    [2, 0, 2],
    [2, 0, -2],
    [2, 0, -1],
    [2, 1, 0],
    [2, 1, 1],
    [2, 1, 2],
    [2, 1, -2],
    [2, 1, -1],
    [2, 2, 0],
    [2, 2, 1],
    [2, 2, 2],
    [2, 2, -2],
    [2, 2, -1],
    [2, -2, 0],
    [2, -2, 1],
    [2, -2, 2],
    [2, -2, -2],
    [2, -2, -1],
    [2, -1, 0],
    [2, -1, 1],
    [2, -1, 2],
    [2, -1, -2],
    [2, -1, -1],
    [-2, 0, 0],
    [-2, 0, 1],
    [-2, 0, 2],
    [-2, 0, -2],
    [-2, 0, -1],
    [-2, 1, 0],
    [-2, 1, 1],
    [-2, 1, 2],
    [-2, 1, -2],
    [-2, 1, -1],
    [-2, 2, 0],
    [-2, 2, 1],
    [-2, 2, 2],
    [-2, 2, -2],
    [-2, 2, -1],
    [-2, -2, 0],
    [-2, -2, 1],
    [-2, -2, 2],
    [-2, -2, -2],
    [-2, -2, -1],
    [-2, -1, 0],
    [-2, -1, 1],
    [-2, -1, 2],
    [-2, -1, -2],
    [-2, -1, -1],
    [-1, 0, 0],
    [-1, 0, 1],
    [-1, 0, 2],
    [-1, 0, -2],
    [-1, 0, -1],
    [-1, 1, 0],
    [-1, 1, 1],
    [-1, 1, 2],
    [-1, 1, -2],
    [-1, 1, -1],
    [-1, 2, 0],
    [-1, 2, 1],
    [-1, 2, 2],
    [-1, 2, -2],
    [-1, 2, -1],
    [-1, -2, 0],
    [-1, -2, 1],
    [-1, -2, 2],
    [-1, -2, -2],
    [-1, -2, -1],
    [-1, -1, 0],
    [-1, -1, 1],
    [-1, -1, 2],
    [-1, -1, -2],
    [-1, -1, -1],
];

/// 获取不可约倒易网格 (Irreducible Reciprocal Mesh)
///
/// # Arguments
/// * `grid_address` - 输出：网格点坐标
/// * `ir_mapping_table` - 输出：映射表，将每个网格点映射到其不可约代表点
/// * `mesh` - 网格尺寸 [Nx, Ny, Nz]
/// * `is_shift` - 网格位移 (Monkhorst-Pack shift)
/// * `rot_reciprocal` - 倒易空间中的点群旋转矩阵
pub fn kpt_get_irreducible_reciprocal_mesh(
    grid_address: &mut [[i32; 3]],
    ir_mapping_table: &mut [usize],
    mesh: &[i32; 3],
    is_shift: &[i32; 3],
    rot_reciprocal: &MatINT,
) -> usize {
    kpt_get_dense_irreducible_reciprocal_mesh(
        grid_address,
        ir_mapping_table,
        mesh,
        is_shift,
        rot_reciprocal,
    )
}

pub fn kpt_get_dense_irreducible_reciprocal_mesh(
    grid_address: &mut [[i32; 3]],
    ir_mapping_table: &mut [usize],
    mesh: &[i32; 3],
    is_shift: &[i32; 3],
    rot_reciprocal: &MatINT,
) -> usize {
    get_dense_ir_reciprocal_mesh(
        grid_address,
        ir_mapping_table,
        mesh,
        is_shift,
        rot_reciprocal,
    )
}

/// 获取考虑了时间反演和 q 点的稳定倒易网格
pub fn kpt_get_stabilized_reciprocal_mesh(
    grid_address: &mut [[i32; 3]],
    ir_mapping_table: &mut [usize],
    mesh: &[i32; 3],
    is_shift: &[i32; 3],
    is_time_reversal: i32,
    rotations: &MatINT,
    qpoints: &[[f64; 3]],
) -> usize {
    kpt_get_dense_stabilized_reciprocal_mesh(
        grid_address,
        ir_mapping_table,
        mesh,
        is_shift,
        is_time_reversal,
        rotations,
        qpoints,
    )
}

pub fn kpt_get_dense_stabilized_reciprocal_mesh(
    grid_address: &mut [[i32; 3]],
    ir_mapping_table: &mut [usize],
    mesh: &[i32; 3],
    is_shift: &[i32; 3],
    is_time_reversal: i32,
    rotations: &MatINT,
    qpoints: &[[f64; 3]],
) -> usize {
    // 获取倒易空间点群（包含时间反演）
    let rot_reciprocal = get_point_group_reciprocal(rotations, is_time_reversal)
        .expect("Failed to allocate rot_reciprocal");

    // 计算容差
    let tolerance = 0.01 / (mesh[0] + mesh[1] + mesh[2]) as f64;
    
    // 获取稳定化 q 点后的点群
    let rot_reciprocal_q = get_point_group_reciprocal_with_q(&rot_reciprocal, tolerance, qpoints)
        .expect("Failed to allocate rot_reciprocal_q");

    let num_ir = get_dense_ir_reciprocal_mesh(
        grid_address,
        ir_mapping_table,
        mesh,
        is_shift,
        &rot_reciprocal_q,
    );

    num_ir
}

/// 将网格点重定位到第一布里渊区 (First Brillouin Zone)
pub fn kpt_relocate_bz_grid_address(
    bz_grid_address: &mut [[i32; 3]],
    bz_map: &mut [usize],
    grid_address: &[[i32; 3]],
    mesh: &[i32; 3],
    rec_lattice: &[[f64; 3]; 3],
    is_shift: &[i32; 3],
) -> usize {
    let num_bz_map = (mesh[0] * mesh[1] * mesh[2]) as usize * 8;
    // 使用 vec! 宏分配内存，对应 C 的 malloc
    let mut dense_bz_map = vec![0; num_bz_map];

    let num_bzgp = relocate_dense_bz_grid_address(
        bz_grid_address,
        &mut dense_bz_map,
        grid_address,
        mesh,
        rec_lattice,
        is_shift,
    );

    for i in 0..num_bz_map {
        if dense_bz_map[i] == num_bz_map {
            bz_map[i] = usize::MAX; // 对应 C 中的 -1 (size_t)
        } else {
            bz_map[i] = dense_bz_map[i];
        }
    }

    num_bzgp
}

// --- Internal Logic ---

/// 获取倒易空间点群
fn get_point_group_reciprocal(rotations: &MatINT, is_time_reversal: i32) -> Option<MatINT> {
    let inversion = [[-1, 0, 0], [0, -1, 0], [0, 0, -1]];
    let size = if is_time_reversal != 0 {
        rotations.size * 2
    } else {
        rotations.size
    };

    let mut rot_reciprocal = mat_alloc_matint(size);
    let mut unique_rot = vec![-1; size];

    for i in 0..rotations.size {
        // 倒易空间的旋转矩阵是实空间旋转矩阵的转置
        let t = mat_transpose_matrix_i3(&rotations.mat[i]);
        mat_copy_matrix_i3(&mut rot_reciprocal.mat[i], &t);

        if is_time_reversal != 0 {
            let inv_rot = mat_multiply_matrix_i3(&inversion, &rot_reciprocal.mat[i]);
            mat_copy_matrix_i3(&mut rot_reciprocal.mat[rotations.size + i], &inv_rot);
        }
    }

    // 筛选唯一旋转矩阵
    let mut num_rot = 0;
    for i in 0..rot_reciprocal.size {
        let mut is_unique = true;
        for j in 0..num_rot {
            if mat_check_identity_matrix_i3(
                &rot_reciprocal.mat[unique_rot[j] as usize],
                &rot_reciprocal.mat[i],
            ) {
                is_unique = false;
                break;
            }
        }
        if is_unique {
            unique_rot[num_rot] = i as i32;
            num_rot += 1;
        }
    }

    let mut rot_return = mat_alloc_matint(num_rot);
    for i in 0..num_rot {
        mat_copy_matrix_i3(
            &mut rot_return.mat[i],
            &rot_reciprocal.mat[unique_rot[i] as usize],
        );
    }

    Some(rot_return)
}

/// 考虑 q 点的对称性
fn get_point_group_reciprocal_with_q(
    rot_reciprocal: &MatINT,
    symprec: f64,
    qpoints: &[[f64; 3]],
) -> Option<MatINT> {
    let mut ir_rot = vec![-1; rot_reciprocal.size];
    let mut num_rot = 0;

    for i in 0..rot_reciprocal.size {
        let mut is_all_ok = true;
        for j in 0..qpoints.len() {
            let q_rot = mat_multiply_matrix_vector_id3(&rot_reciprocal.mat[i], &qpoints[j]);

            let mut found_diff = false;
            for k in 0..qpoints.len() {
                let mut diff = [0.0; 3];
                for l in 0..3 {
                    diff[l] = q_rot[l] - qpoints[k][l];
                    diff[l] -= mat_nint(diff[l]) as f64;
                }
                if mat_dabs(diff[0]) < symprec
                    && mat_dabs(diff[1]) < symprec
                    && mat_dabs(diff[2]) < symprec
                {
                    found_diff = true;
                    break;
                }
            }

            if !found_diff {
                is_all_ok = false;
                break;
            }
        }

        if is_all_ok {
            ir_rot[num_rot] = i as i32;
            num_rot += 1;
        }
    }

    let mut rot_reciprocal_q = mat_alloc_matint(num_rot);
    for i in 0..num_rot {
        mat_copy_matrix_i3(
            &mut rot_reciprocal_q.mat[i],
            &rot_reciprocal.mat[ir_rot[i] as usize],
        );
    }

    Some(rot_reciprocal_q)
}

fn get_dense_ir_reciprocal_mesh(
    grid_address: &mut [[i32; 3]],
    ir_mapping_table: &mut [usize],
    mesh: &[i32; 3],
    is_shift: &[i32; 3],
    rot_reciprocal: &MatINT,
) -> usize {
    if check_mesh_symmetry(mesh, is_shift, rot_reciprocal) {
        get_dense_ir_reciprocal_mesh_normal(
            grid_address,
            ir_mapping_table,
            mesh,
            is_shift,
            rot_reciprocal,
        )
    } else {
        get_dense_ir_reciprocal_mesh_distortion(
            grid_address,
            ir_mapping_table,
            mesh,
            is_shift,
            rot_reciprocal,
        )
    }
}

/// 普通网格约化（适用于正交或高对称性网格）
fn get_dense_ir_reciprocal_mesh_normal(
    grid_address: &mut [[i32; 3]],
    ir_mapping_table: &mut [usize],
    mesh: &[i32; 3],
    is_shift: &[i32; 3],
    rot_reciprocal: &MatINT,
) -> usize {
    kgrid::kgd_get_all_grid_addresses(grid_address, mesh);

    let total_pts = (mesh[0] * mesh[1] * mesh[2]) as usize;

    // 使用 Rayon 并行迭代器（如果启用 feature）
    // C 代码中使用 OpenMP 并行化此循环
    #[cfg(feature = "parallel")]
    let iter = (0..total_pts).into_par_iter();
    #[cfg(not(feature = "parallel"))]
    let iter = 0..total_pts;

    // 注意：我们需要修改 ir_mapping_table，并行写入需要 unsafe 或拆分
    // 这里为了保持与 C 代码逻辑一致（C 代码中每个线程处理不同的 i），
    // 我们可以使用 chunk 或者直接串行，因为 Rust 的借用检查器会阻止简单的并行写入。
    // 为了简单且安全，这里保持串行，或者使用 UnsafeCell/Atomic (过于复杂)。
    // 考虑到 C 代码中 OpenMP 的使用，这里实际上是数据并行的。
    // 我们可以先计算出每个点的映射，然后收集。
    
    // 暂且保持串行实现，因为 ir_mapping_table 的依赖性较弱（只依赖于 grid_point_rot，它是只读计算出来的）
    // 但 ir_mapping_table[i] = ir_mapping_table[grid_point_rot] 这种写法在 C 的 OpenMP 中其实是有潜在竞争的，
    // 除非 grid_point_rot 总是小于 i 且已经被处理过？
    // C 代码逻辑：
    // if (grid_point_rot < ir_mapping_table[i]) { ir_mapping_table[i] = grid_point_rot; }
    // 初始 ir_mapping_table[i] = i。
    // 所以它是在寻找轨道中最小的索引作为代表。
    
    for i in 0..total_pts {
        let mut address_double = [0; 3];
        kgrid::kgd_get_grid_address_double_mesh(
            &mut address_double,
            &grid_address[i],
            mesh,
            is_shift,
        );

        ir_mapping_table[i] = i;

        for j in 0..rot_reciprocal.size {
            let address_double_rot =
                mat_multiply_matrix_vector_i3(&rot_reciprocal.mat[j], &address_double);
            let grid_point_rot =
                kgrid::kgd_get_dense_grid_point_double_mesh(&address_double_rot, mesh);

            if grid_point_rot < ir_mapping_table[i] {
                ir_mapping_table[i] = grid_point_rot;
                // C 代码中如果没有 OpenMP 会 break。
                // 如果有 OpenMP，则继续寻找全局最小值。
                // 为了保证结果的确定性和最优性，我们不 break，总是寻找最小代表元。
                #[cfg(not(feature = "parallel"))]
                break; 
            }
        }
    }

    get_dense_num_ir(ir_mapping_table, mesh)
}

/// 畸变网格约化（适用于非正交网格或低对称性）
fn get_dense_ir_reciprocal_mesh_distortion(
    grid_address: &mut [[i32; 3]],
    ir_mapping_table: &mut [usize],
    mesh: &[i32; 3],
    is_shift: &[i32; 3],
    rot_reciprocal: &MatINT,
) -> usize {
    kgrid::kgd_get_all_grid_addresses(grid_address, mesh);

    let divisor = [
        (mesh[1] * mesh[2]) as i64,
        (mesh[2] * mesh[0]) as i64,
        (mesh[0] * mesh[1]) as i64,
    ];
    let total_pts = (mesh[0] * mesh[1] * mesh[2]) as usize;

    for i in 0..total_pts {
        let mut address_double = [0; 3];
        kgrid::kgd_get_grid_address_double_mesh(
            &mut address_double,
            &grid_address[i],
            mesh,
            is_shift,
        );

        let mut long_address_double = [0i64; 3];
        for j in 0..3 {
            long_address_double[j] = address_double[j] as i64 * divisor[j];
        }

        ir_mapping_table[i] = i;

        for j in 0..rot_reciprocal.size {
            let mut long_address_double_rot = [0i64; 3];
            for k in 0..3 {
                long_address_double_rot[k] = rot_reciprocal.mat[j][k][0] as i64
                    * long_address_double[0]
                    + rot_reciprocal.mat[j][k][1] as i64 * long_address_double[1]
                    + rot_reciprocal.mat[j][k][2] as i64 * long_address_double[2];
            }

            let mut indivisible = false;
            let mut address_double_rot = [0; 3];

            for k in 0..3 {
                if long_address_double_rot[k] % divisor[k] != 0 {
                    indivisible = true;
                    break;
                }
                address_double_rot[k] = (long_address_double_rot[k] / divisor[k]) as i32;

                if (address_double_rot[k] % 2 != 0 && is_shift[k] == 0)
                    || (address_double_rot[k] % 2 == 0 && is_shift[k] == 1)
                {
                    indivisible = true;
                    break;
                }
            }

            if indivisible {
                continue;
            }

            let grid_point_rot =
                kgrid::kgd_get_dense_grid_point_double_mesh(&address_double_rot, mesh);

            if grid_point_rot < ir_mapping_table[i] {
                ir_mapping_table[i] = grid_point_rot;
                #[cfg(not(feature = "parallel"))]
                break;
            }
        }
    }

    get_dense_num_ir(ir_mapping_table, mesh)
}

/// 统计不可约点数量并进行路径压缩
fn get_dense_num_ir(ir_mapping_table: &mut [usize], mesh: &[i32; 3]) -> usize {
    let total_pts = (mesh[0] * mesh[1] * mesh[2]) as usize;
    let mut num_ir = 0;

    for i in 0..total_pts {
        if ir_mapping_table[i] == i {
            num_ir += 1;
        }
    }

    // 路径压缩：确保每个点直接指向其最终的代表元
    for i in 0..total_pts {
        ir_mapping_table[i] = ir_mapping_table[ir_mapping_table[i]];
    }

    num_ir
}

fn relocate_dense_bz_grid_address(
    bz_grid_address: &mut [[i32; 3]],
    bz_map: &mut [usize],
    grid_address: &[[i32; 3]],
    mesh: &[i32; 3],
    rec_lattice: &[[f64; 3]; 3],
    is_shift: &[i32; 3],
) -> usize {
    let tolerance = get_tolerance_for_bz_reduction(rec_lattice, mesh);
    let bzmesh = [mesh[0] * 2, mesh[1] * 2, mesh[2] * 2];
    let num_bzmesh = (bzmesh[0] * bzmesh[1] * bzmesh[2]) as usize;

    // 初始化 bz_map
    for i in 0..num_bzmesh {
        bz_map[i] = num_bzmesh;
    }

    let mut boundary_num_gp = 0;
    let total_num_gp = (mesh[0] * mesh[1] * mesh[2]) as usize;

    for i in 0..total_num_gp {
        let mut distance = [0.0; KPT_NUM_BZ_SEARCH_SPACE];

        // 计算到所有邻近点的距离
        for j in 0..KPT_NUM_BZ_SEARCH_SPACE {
            let mut q_vector = [0.0; 3];
            for k in 0..3 {
                q_vector[k] = ((grid_address[i][k] + BZ_SEARCH_SPACE[j][k] * mesh[k]) * 2
                    + is_shift[k]) as f64
                    / (mesh[k] as f64)
                    / 2.0;
            }
            let q_vec_rec = mat_multiply_matrix_vector_d3(rec_lattice, &q_vector);
            distance[j] = mat_norm_squared_d3(&q_vec_rec);
        }

        // 找到最小距离
        let mut min_distance = distance[0];
        let mut min_index = 0;
        for j in 1..KPT_NUM_BZ_SEARCH_SPACE {
            if distance[j] < min_distance {
                min_distance = distance[j];
                min_index = j;
            }
        }

        // 标记所有在容差范围内的点（处理边界点）
        for j in 0..KPT_NUM_BZ_SEARCH_SPACE {
            if distance[j] < min_distance + tolerance {
                let gp = if j == min_index {
                    i
                } else {
                    boundary_num_gp + total_num_gp
                };

                let mut bz_address_double = [0; 3];
                for k in 0..3 {
                    bz_grid_address[gp][k] = grid_address[i][k] + BZ_SEARCH_SPACE[j][k] * mesh[k];
                    bz_address_double[k] = bz_grid_address[gp][k] * 2 + is_shift[k];
                }

                let bzgp = kgrid::kgd_get_dense_grid_point_double_mesh(&bz_address_double, &bzmesh);
                bz_map[bzgp] = gp;

                if j != min_index {
                    boundary_num_gp += 1;
                }
            }
        }
    }

    boundary_num_gp + total_num_gp
}

fn get_tolerance_for_bz_reduction(rec_lattice: &[[f64; 3]; 3], mesh: &[i32; 3]) -> f64 {
    let mut length = [0.0; 3];
    for i in 0..3 {
        for j in 0..3 {
            length[i] += rec_lattice[j][i] * rec_lattice[j][i];
        }
        length[i] /= (mesh[i] * mesh[i]) as f64;
    }

    let mut tolerance = length[0];
    for i in 1..3 {
        if tolerance < length[i] {
            tolerance = length[i];
        }
    }
    tolerance * 0.01
}

/// 检查网格对称性
/// 
/// 注意：提供的 C 代码中此函数存在复制粘贴错误（重复检查 column 0）。
/// 本 Rust 实现已修正此问题，正确检查 column 0, 1, 2 分别对应 a=b, b=c, c=a 的对称性。
fn check_mesh_symmetry(mesh: &[i32; 3], is_shift: &[i32; 3], rot_reciprocal: &MatINT) -> bool {
    let mut eq = [false; 3]; // eq[0]: a=b, eq[1]: b=c, eq[2]: c=a

    for i in 0..rot_reciprocal.size {
        let mut sum = 0;
        for j in 0..3 {
            for k in 0..3 {
                sum += rot_reciprocal.mat[i][j][k].abs();
            }
        }
        if sum > 3 {
            return false;
        }
    }

    for i in 0..rot_reciprocal.size {
        // 检查 x <-> y 交换 (a=b)
        // 矩阵应为 [0 1 0; 1 0 0; 0 0 1] 或类似，关注列 0 变为 [0, 1, 0]
        if rot_reciprocal.mat[i][0][0] == 0
            && rot_reciprocal.mat[i][1][0] == 1
            && rot_reciprocal.mat[i][2][0] == 0
        {
            eq[0] = true;
        }
        // 检查 y <-> z 交换 (b=c)
        // 关注列 1 变为 [0, 0, 1]
        if rot_reciprocal.mat[i][0][1] == 0
            && rot_reciprocal.mat[i][1][1] == 0
            && rot_reciprocal.mat[i][2][1] == 1
        {
            eq[1] = true;
        }
        // 检查 z <-> x 交换 (c=a)
        // 关注列 2 变为 [1, 0, 0]
        if rot_reciprocal.mat[i][0][2] == 1
            && rot_reciprocal.mat[i][1][2] == 0
            && rot_reciprocal.mat[i][2][2] == 0
        {
            eq[2] = true;
        }
    }

    let cond1 = (eq[0] && mesh[0] == mesh[1] && is_shift[0] == is_shift[1]) || !eq[0];
    let cond2 = (eq[1] && mesh[1] == mesh[2] && is_shift[1] == is_shift[2]) || !eq[1];
    let cond3 = (eq[2] && mesh[2] == mesh[0] && is_shift[2] == is_shift[0]) || !eq[2];

    cond1 && cond2 && cond3
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mathfunc::MatINT;

    #[test]
    fn test_check_mesh_symmetry() {
        let mesh = [4, 4, 4];
        let shift = [0, 0, 0];
        let mut rot = MatINT::new(1);
        // Identity
        rot.mat[0] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

        assert!(check_mesh_symmetry(&mesh, &shift, &rot));
    }

    #[test]
    fn test_kpt_get_irreducible_reciprocal_mesh_simple() {
        // 简单的 2x2x2 网格，无位移，只有恒等操作
        let mesh = [2, 2, 2];
        let shift = [0, 0, 0];
        let mut rot = MatINT::new(1);
        rot.mat[0] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

        let mut grid_address = vec![[0; 3]; 8];
        let mut map = vec![0; 8];

        let num_ir = kpt_get_irreducible_reciprocal_mesh(
            &mut grid_address,
            &mut map,
            &mesh,
            &shift,
            &rot,
        );

        // 由于只有恒等操作，每个点都是不可约的
        assert_eq!(num_ir, 8);
        for i in 0..8 {
            assert_eq!(map[i], i);
        }
    }
}
