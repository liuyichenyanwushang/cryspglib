// mathfunc.rs

use crate::debug;

// 定义常量，与 C 代码保持一致
const ZERO_PREC: f64 = 1e-10;

// 类型别名，方便阅读
pub type Mat3 = [[f64; 3]; 3];
pub type Mat3I = [[i32; 3]; 3];
pub type Vec3 = [f64; 3];
pub type Vec3I = [i32; 3];

/// 计算 3x3 double 矩阵的行列式
#[inline]
pub fn mat_get_determinant_d3(a: &Mat3) -> f64 {
    a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        + a[0][1] * (a[1][2] * a[2][0] - a[1][0] * a[2][2])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
}

/// 计算 3x3 int 矩阵的行列式
#[inline]
pub fn mat_get_determinant_i3(a: &Mat3I) -> i32 {
    a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        + a[0][1] * (a[1][2] * a[2][0] - a[1][0] * a[2][2])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
}

/// 计算 3x3 int 矩阵的迹 (Trace)
#[inline]
pub fn mat_get_trace_i3(a: &Mat3I) -> i32 {
    a[0][0] + a[1][1] + a[2][2]
}

/// 复制 3x3 double 矩阵
#[inline]
pub fn mat_copy_matrix_d3(a: &mut Mat3, b: &Mat3) {
    *a = *b;
}

/// 复制 3x3 int 矩阵
#[inline]
pub fn mat_copy_matrix_i3(a: &mut Mat3I, b: &Mat3I) {
    *a = *b;
}

/// 复制 double 向量
#[inline]
pub fn mat_copy_vector_d3(a: &mut Vec3, b: &Vec3) {
    *a = *b;
}

/// 复制 int 向量
#[inline]
pub fn mat_copy_vector_i3(a: &mut Vec3I, b: &Vec3I) {
    *a = *b;
}

/// 检查两个 int 矩阵是否相同
#[inline]
pub fn mat_check_identity_matrix_i3(a: &Mat3I, b: &Mat3I) -> bool {
    a == b
}

/// 检查两个 double 矩阵是否在误差范围内相同
pub fn mat_check_identity_matrix_d3(a: &Mat3, b: &Mat3, symprec: f64) -> bool {
    for i in 0..3 {
        for j in 0..3 {
            if (a[i][j] - b[i][j]).abs() > symprec {
                return false;
            }
        }
    }
    true
}

/// 检查 int 矩阵和 double 矩阵是否在误差范围内相同
pub fn mat_check_identity_matrix_id3(a: &Mat3I, b: &Mat3, symprec: f64) -> bool {
    for i in 0..3 {
        for j in 0..3 {
            if (a[i][j] as f64 - b[i][j]).abs() > symprec {
                return false;
            }
        }
    }
    true
}

/// 矩阵乘法 m = a x b (double)
#[inline]
pub fn mat_multiply_matrix_d3(a: &Mat3, b: &Mat3) -> Mat3 {
    let mut c = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j];
        }
    }
    c
}

/// 矩阵乘法 m = a x b (int)
#[inline]
pub fn mat_multiply_matrix_i3(a: &Mat3I, b: &Mat3I) -> Mat3I {
    let mut c = [[0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j];
        }
    }
    c
}

/// 矩阵乘法 m = a x b (double x int -> double)
#[inline]
pub fn mat_multiply_matrix_di3(a: &Mat3, b: &Mat3I) -> Mat3 {
    let mut c = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = a[i][0] * (b[0][j] as f64)
                + a[i][1] * (b[1][j] as f64)
                + a[i][2] * (b[2][j] as f64);
        }
    }
    c
}

/// 矩阵乘法 m = a x b (int x double -> double)
#[inline]
pub fn mat_multiply_matrix_id3(a: &Mat3I, b: &Mat3) -> Mat3 {
    let mut c = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = (a[i][0] as f64) * b[0][j]
                + (a[i][1] as f64) * b[1][j]
                + (a[i][2] as f64) * b[2][j];
        }
    }
    c
}

/// 矩阵向量乘法 v = a x b (int)
#[inline]
pub fn mat_multiply_matrix_vector_i3(a: &Mat3I, b: &Vec3I) -> Vec3I {
    let mut c = [0; 3];
    for i in 0..3 {
        c[i] = a[i][0] * b[0] + a[i][1] * b[1] + a[i][2] * b[2];
    }
    c
}

/// 矩阵向量乘法 v = a x b (double)
#[inline]
pub fn mat_multiply_matrix_vector_d3(a: &Mat3, b: &Vec3) -> Vec3 {
    let mut c = [0.0; 3];
    for i in 0..3 {
        c[i] = a[i][0] * b[0] + a[i][1] * b[1] + a[i][2] * b[2];
    }
    c
}

/// 矩阵向量乘法 v = a x b (int x double -> double)
#[inline]
pub fn mat_multiply_matrix_vector_id3(a: &Mat3I, b: &Vec3) -> Vec3 {
    let mut c = [0.0; 3];
    for i in 0..3 {
        c[i] = (a[i][0] as f64) * b[0] + (a[i][1] as f64) * b[1] + (a[i][2] as f64) * b[2];
    }
    c
}

/// 矩阵向量乘法 v = a x b (double x int -> double)
#[inline]
pub fn mat_multiply_matrix_vector_di3(a: &Mat3, b: &Vec3I) -> Vec3 {
    let mut c = [0.0; 3];
    for i in 0..3 {
        c[i] = a[i][0] * (b[0] as f64) + a[i][1] * (b[1] as f64) + a[i][2] * (b[2] as f64);
    }
    c
}

/// 矩阵加法 m = a + b (int)
#[inline]
pub fn mat_add_matrix_i3(a: &Mat3I, b: &Mat3I) -> Mat3I {
    let mut m = [[0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            m[i][j] = a[i][j] + b[i][j];
        }
    }
    m
}

/// 类型转换 int -> double
#[inline]
pub fn mat_cast_matrix_3i_to_3d(a: &Mat3I) -> Mat3 {
    let mut m = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            m[i][j] = a[i][j] as f64;
        }
    }
    m
}

/// 类型转换 double -> int (四舍五入)
#[inline]
pub fn mat_cast_matrix_3d_to_3i(a: &Mat3) -> Mat3I {
    let mut m = [[0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            m[i][j] = mat_nint(a[i][j]);
        }
    }
    m
}

/// 计算逆矩阵 m = a^-1
/// 如果行列式接近 0，返回 None
pub fn mat_inverse_matrix_d3(a: &Mat3, precision: f64) -> Option<Mat3> {
    let det = mat_get_determinant_d3(a);
    if det.abs() < precision {
        debug::debug_print(format_args!("spglib: No inverse matrix (det={})\n", det));
        return None;
    }

    let mut c = [[0.0; 3]; 3];
    c[0][0] = (a[1][1] * a[2][2] - a[1][2] * a[2][1]) / det;
    c[1][0] = (a[1][2] * a[2][0] - a[1][0] * a[2][2]) / det;
    c[2][0] = (a[1][0] * a[2][1] - a[1][1] * a[2][0]) / det;
    c[0][1] = (a[2][1] * a[0][2] - a[2][2] * a[0][1]) / det;
    c[1][1] = (a[2][2] * a[0][0] - a[2][0] * a[0][2]) / det;
    c[2][1] = (a[2][0] * a[0][1] - a[2][1] * a[0][0]) / det;
    c[0][2] = (a[0][1] * a[1][2] - a[0][2] * a[1][1]) / det;
    c[1][2] = (a[0][2] * a[1][0] - a[0][0] * a[1][2]) / det;
    c[2][2] = (a[0][0] * a[1][1] - a[0][1] * a[1][0]) / det;

    Some(c)
}

/// 计算相似矩阵 m = b^-1 a b
pub fn mat_get_similar_matrix_d3(a: &Mat3, b: &Mat3, precision: f64) -> Option<Mat3> {
    let inv_b = mat_inverse_matrix_d3(b, precision)?;
    // m = a * b
    let temp = mat_multiply_matrix_d3(a, b);
    // m = inv_b * temp
    Some(mat_multiply_matrix_d3(&inv_b, &temp))
}

/// 矩阵转置 (double)
#[inline]
pub fn mat_transpose_matrix_d3(b: &Mat3) -> Mat3 {
    let mut c = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = b[j][i];
        }
    }
    c
}

/// 矩阵转置 (int)
#[inline]
pub fn mat_transpose_matrix_i3(b: &Mat3I) -> Mat3I {
    let mut c = [[0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = b[j][i];
        }
    }
    c
}

/// 计算度量张量 metric = lattice^T * lattice
pub fn mat_get_metric(lattice: &Mat3) -> Mat3 {
    let lattice_t = mat_transpose_matrix_d3(lattice);
    mat_multiply_matrix_d3(&lattice_t, lattice)
}

/// 向量范数平方 (double)
#[inline]
pub fn mat_norm_squared_d3(a: &Vec3) -> f64 {
    a[0] * a[0] + a[1] * a[1] + a[2] * a[2]
}

/// 向量范数平方 (int)
#[inline]
pub fn mat_norm_squared_i3(a: &Vec3I) -> i32 {
    a[0] * a[0] + a[1] * a[1] + a[2] * a[2]
}

/// 向量叉积 v = a x b
#[inline]
pub fn mat_cross_product_d3(a: &Vec3, b: &Vec3) -> Vec3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// 绝对值
#[inline]
pub fn mat_dabs(a: f64) -> f64 {
    a.abs()
}

/// 四舍五入到最近的整数
/// C 代码逻辑: a < 0 ? (int)(a - 0.5) : (int)(a + 0.5)
/// Rust 的 as i32 是向零取整，所以逻辑一致
#[inline]
pub fn mat_nint(a: f64) -> i32 {
    if a < 0.0 {
        (a - 0.5) as i32
    } else {
        (a + 0.5) as i32
    }
}

/// 对 1 取模，返回 [0, 1) 范围内的值
/// 注意：C 代码中如果结果 < -ZERO_PREC 会 +1.0
#[inline]
pub fn mat_dmod1(a: f64) -> f64 {
    let b = a - mat_nint(a) as f64;
    if b < -ZERO_PREC { b + 1.0 } else { b }
}

/// 返回 a 除以 1 的余数，范围 [-0.5, 0.5]
#[inline]
pub fn mat_rem1(a: f64) -> f64 {
    a - mat_nint(a) as f64
}

/// 检查矩阵是否为整数矩阵（在误差范围内）
pub fn mat_is_int_matrix(mat: &Mat3, symprec: f64) -> bool {
    for i in 0..3 {
        for j in 0..3 {
            if (mat_nint(mat[i][j]) as f64 - mat[i][j]).abs() > symprec {
                return false;
            }
        }
    }
    true
}

// --- 动态数组结构体封装 (对应 C 的 MatINT 和 VecDBL) ---

/// 对应 C 的 MatINT 结构
/// 在 Rust 中通常直接使用 Vec<Mat3I>，这里为了兼容性保留结构体
#[derive(Debug, Clone)]
pub struct MatINT {
    pub size: usize,
    pub mat: Vec<Mat3I>,
}

impl MatINT {
    pub fn new(size: usize) -> Self {
        // 预分配内存，初始化为 0
        MatINT {
            size,
            mat: vec![[[0; 3]; 3]; size],
        }
    }
}

/// 对应 C 的 VecDBL 结构
#[derive(Debug, Clone)]
pub struct VecDBL {
    pub size: usize,
    pub vec: Vec<Vec3>,
}

impl VecDBL {
    pub fn new(size: usize) -> Self {
        VecDBL {
            size,
            vec: vec![[0.0; 3]; size],
        }
    }
}

// 辅助函数：模拟 C 的 mat_alloc_MatINT
// 建议直接使用 MatINT::new
pub fn mat_alloc_matint(size: usize) -> MatINT {
    MatINT::new(size)
}

// 辅助函数：模拟 C 的 mat_alloc_VecDBL
// 建议直接使用 VecDBL::new
pub fn mat_alloc_vecdbl(size: usize) -> VecDBL {
    VecDBL::new(size)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_determinant() {
        let a: Mat3 = [[1.0, 2.0, 3.0], [0.0, 1.0, 4.0], [5.0, 6.0, 0.0]];
        // det = 1*(0-24) - 2*(0-20) + 3*(0-5) = -24 + 40 - 15 = 1
        assert!((mat_get_determinant_d3(&a) - 1.0).abs() < 1e-10);

        let b: Mat3I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        assert_eq!(mat_get_determinant_i3(&b), 1);
    }

    #[test]
    fn test_trace() {
        let a: Mat3I = [[1, 2, 3], [4, 5, 6], [7, 8, 9]];
        assert_eq!(mat_get_trace_i3(&a), 1 + 5 + 9);
    }

    #[test]
    fn test_copy() {
        let src: Mat3 = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        let mut dst: Mat3 = [[0.0; 3]; 3];
        mat_copy_matrix_d3(&mut dst, &src);
        assert_eq!(src, dst);
    }

    #[test]
    fn test_check_identity() {
        let a: Mat3I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        let b: Mat3I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        let c: Mat3I = [[1, 0, 0], [0, 0, 0], [0, 0, 1]];
        assert!(mat_check_identity_matrix_i3(&a, &b));
        assert!(!mat_check_identity_matrix_i3(&a, &c));

        let d: Mat3 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        assert!(mat_check_identity_matrix_d3(&d, &d, 1e-5));
        assert!(mat_check_identity_matrix_id3(&a, &d, 1e-5));
    }

    #[test]
    fn test_multiply_matrix() {
        let a: Mat3I = [[1, 2, 3], [4, 5, 6], [7, 8, 9]];
        let b: Mat3I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        let res = mat_multiply_matrix_i3(&a, &b);
        assert_eq!(res, a);

        let ad: Mat3 = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        let bd: Mat3 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let res_d = mat_multiply_matrix_d3(&ad, &bd);
        assert_eq!(res_d, ad);
    }

    #[test]
    fn test_multiply_vector() {
        let a: Mat3I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        let v: Vec3I = [1, 2, 3];
        let res = mat_multiply_matrix_vector_i3(&a, &v);
        assert_eq!(res, v);
    }

    #[test]
    fn test_inverse() {
        let a: Mat3 = [[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 4.0]];
        let inv = mat_inverse_matrix_d3(&a, 1e-10).unwrap();
        assert!((inv[0][0] - 1.0).abs() < 1e-10);
        assert!((inv[1][1] - 0.5).abs() < 1e-10);
        assert!((inv[2][2] - 0.25).abs() < 1e-10);

        let singular: Mat3 = [[0.0; 3]; 3];
        assert!(mat_inverse_matrix_d3(&singular, 1e-10).is_none());
    }

    #[test]
    fn test_similar_matrix() {
        // m = b^-1 a b
        // Let b = I, then m = a
        let a: Mat3 = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        let b: Mat3 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let m = mat_get_similar_matrix_d3(&a, &b, 1e-10).unwrap();
        assert_eq!(m, a);
    }

    #[test]
    fn test_transpose() {
        let a: Mat3I = [[1, 2, 3], [4, 5, 6], [7, 8, 9]];
        let t = mat_transpose_matrix_i3(&a);
        assert_eq!(t[0][1], 4);
        assert_eq!(t[1][0], 2);
    }

    #[test]
    fn test_norm_cross() {
        let v: Vec3 = [1.0, 2.0, 2.0];
        assert_eq!(mat_norm_squared_d3(&v), 9.0);

        let v1: Vec3 = [1.0, 0.0, 0.0];
        let v2: Vec3 = [0.0, 1.0, 0.0];
        let cross = mat_cross_product_d3(&v1, &v2);
        assert_eq!(cross, [0.0, 0.0, 1.0]);
    }

    #[test]
    fn test_nint() {
        assert_eq!(mat_nint(1.1), 1);
        assert_eq!(mat_nint(1.5), 2);
        assert_eq!(mat_nint(-1.1), -1);
        assert_eq!(mat_nint(-1.5), -2); // C logic: -1.5 - 0.5 = -2.0 -> -2
    }

    #[test]
    fn test_dmod1() {
        // 1.2 -> 0.2
        assert!((mat_dmod1(1.2) - 0.2).abs() < 1e-10);
        // -0.2 -> -0.2 - (-0) = -0.2 -> < 0 -> -0.2 + 1.0 = 0.8
        assert!((mat_dmod1(-0.2) - 0.8).abs() < 1e-10);
    }

    #[test]
    fn test_alloc_structs() {
        let m = mat_alloc_matint(5);
        assert_eq!(m.size, 5);
        assert_eq!(m.mat.len(), 5);

        let v = mat_alloc_vecdbl(3);
        assert_eq!(v.size, 3);
        assert_eq!(v.vec.len(), 3);
        
        // No explicit free needed
    }
}
