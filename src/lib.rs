//! rspglib — Rust 晶体空间群识别库。
//!
//! 基于 [spglib](https://github.com/spglib/spglib) 的纯 Rust 移植，
//! 提供晶体对称性分析、空间群识别、标准晶胞构造和 k 点网格生成。
//!
//! # 主要入口
//!
//! - [`spg_get_dataset`] — 获取完整空间群信息数据集
//! - [`spg_get_symmetry`] — 仅获取对称操作
//! - [`spg_standardize_cell`] — 获取理想化标准晶胞
//! - [`spg_find_primitive`] — 获取原胞
//!
//! # 晶格矩阵约定
//!
//! 所有 3x3 矩阵采用 `lattice[cart][vec]` 布局（行=笛卡尔分量，列=晶格矢量）。
//! 详见 [`mathfunc`] 模块文档。

pub mod arithmetic;
pub mod cell;
pub mod debug;
pub mod delaunay;
pub mod determination;
pub mod hall_symbol;
pub mod kgrid;
pub mod kpoint;
pub mod magnetic_spacegroup;
pub mod mathfunc;
pub mod msg_database;
#[cfg(test)]
pub mod magnetic_spacegroup_test;
pub mod niggli;
pub mod overlap;
pub mod pointgroup;
pub mod primitive;
pub mod refinement;
pub mod site_symmetry;
pub mod sitesym_database;
pub mod spacegroup;
pub mod spg_database;
pub mod spin;
pub mod symmetry;

use crate::cell::{cel_any_overlap_with_same_type, cel_layer_any_overlap_with_same_type, AperiodicAxis, Cell, TensorRank};
use crate::delaunay::del_delaunay_reduce;
use crate::determination::det_determine_all;
use crate::mathfunc::{
    mat_copy_matrix_d3, mat_copy_vector_d3,
    mat_inverse_matrix_d3, mat_multiply_matrix_d3, Mat3, Mat3I, Vec3,
};
use crate::niggli::niggli_reduce;
use crate::pointgroup::{ptg_get_pointgroup, ptg_get_transformation_matrix};
use crate::primitive::{Primitive, prm_get_primitive_symmetry};
use crate::spacegroup::{
    Spacegroup, spa_search_spacegroup_with_symmetry, spa_transform_from_primitive,
    spa_transform_to_primitive,
};
use crate::spg_database::{Centering, spgdb_get_spacegroup_operations, spgdb_get_spacegroup_type};
use crate::symmetry::Symmetry;

// ---------------------------------------------------------------------------
// Version constants
// ---------------------------------------------------------------------------
/// 主版本号
pub const SPGLIB_MAJOR_VERSION: i32 = 2;
/// 次版本号
pub const SPGLIB_MINOR_VERSION: i32 = 5;
/// 补丁版本号
pub const SPGLIB_MICRO_VERSION: i32 = 4;
/// 版本字符串
pub const SPGLIB_VERSION: &str = "2.5.4";
/// 完整版本字符串（含构建信息）
pub const SPGLIB_VERSION_FULL: &str = "2.5.4";
/// Git 提交哈希
pub const SPGLIB_COMMIT: &str = "unknown";

// ---------------------------------------------------------------------------
// Error codes
// ---------------------------------------------------------------------------
/// spglib 错误码。
#[derive(thiserror::Error, Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum SpglibError {
    /// 无错误
    #[error("no error")]
    Success = 0,
    /// 空间群搜索失败
    #[error("spacegroup search failed")]
    SpacegroupSearchFailed = 1,
    /// 晶胞标准化失败
    #[error("cell standardization failed")]
    CellStandardizationFailed = 2,
    /// 对称操作搜索失败
    #[error("symmetry operation search failed")]
    SymmetryOperationSearchFailed = 3,
    /// 原子间距过近
    #[error("too close distance between atoms")]
    AtomsTooClose = 4,
    /// 点群未找到
    #[error("pointgroup not found")]
    PointgroupNotFound = 5,
    /// Niggli 约化失败
    #[error("Niggli reduction failed")]
    NiggliFailed = 6,
    /// Delaunay 约化失败
    #[error("Delaunay reduction failed")]
    DelaunayFailed = 7,
    /// 数组大小不足
    #[error("array size shortage")]
    ArraySizeShortage = 8,
}

// ---------------------------------------------------------------------------
// Public data structures
// ---------------------------------------------------------------------------

/// 空间群数据集的完整结构。
///
/// 包含标准晶胞、对称操作、Wyckoff 位置标记和映射信息。
/// 所有动态数据由 Rust 的 [`Vec`] 所有权管理，无需手动释放。
#[derive(Debug, Clone)]
pub struct SpglibDataset {
    /// 空间群编号 (1–230)
    pub spacegroup_number: i32,
    /// Hall 编号 (1–530)
    pub hall_number: i32,
    /// 国际符号 (最多 11 字符)
    pub international_symbol: String,
    /// Hall 符号 (最多 17 字符)
    pub hall_symbol: String,
    /// 选择 (最多 6 字符)
    pub choice: String,
    /// 变换矩阵 (Bravais → 原始晶胞)
    pub transformation_matrix: Mat3,
    /// 原点平移
    pub origin_shift: Vec3,
    /// 对称操作数量
    pub n_operations: usize,
    /// 旋转矩阵 [n_operations][3][3]
    pub rotations: Vec<Mat3I>,
    /// 平移矢量 [n_operations][3]
    pub translations: Vec<Vec3>,
    /// 原子数
    pub n_atoms: usize,
    /// Wyckoff 字母编码 (0=a, 1=b, ..., 26=z)
    pub wyckoffs: Vec<i32>,
    /// 位点对称性符号
    pub site_symmetry_symbols: Vec<String>,
    /// 对等原子映射
    pub equivalent_atoms: Vec<i32>,
    /// 晶体学轨道
    pub crystallographic_orbits: Vec<i32>,
    /// 原子 → 原胞映射
    pub mapping_to_primitive: Vec<i32>,
    /// 标准晶胞原子数
    pub n_std_atoms: usize,
    /// 标准晶胞晶格
    pub std_lattice: Mat3,
    /// 标准晶胞原子位置
    pub std_positions: Vec<Vec3>,
    /// 标准晶胞原子类型
    pub std_types: Vec<i32>,
    /// 标准晶胞旋转矩阵
    pub std_rotation_matrix: Mat3,
    /// 标准晶胞 → 原胞映射
    pub std_mapping_to_primitive: Vec<i32>,
    /// 原胞晶格
    pub primitive_lattice: Mat3,
    /// 点群符号 (最多 6 字符)
    pub pointgroup_symbol: String,
}

/// 空间群类型信息（从数据库查询）。
#[derive(Debug, Clone)]
pub struct SpglibSpacegroupType {
    /// 空间群编号 (1–230)
    pub number: i32,
    /// Hall 编号
    pub hall_number: i32,
    /// Schoenflies 符号
    pub schoenflies: String,
    /// Hall 符号
    pub hall_symbol: String,
    /// 选择
    pub choice: String,
    /// 国际符号（完整）
    pub international: String,
    /// 国际符号（完整，多行格式）
    pub international_full: String,
    /// 国际符号（短格式）
    pub international_short: String,
    /// 点群国际符号
    pub pointgroup_international: String,
    /// 点群 Schoenflies 符号
    pub pointgroup_schoenflies: String,
    /// 算术晶体类编号
    pub arithmetic_crystal_class_number: i32,
    /// 算术晶体类符号
    pub arithmetic_crystal_class_symbol: String,
}

/// 磁性空间群数据集。
#[derive(Debug, Clone)]
pub struct SpglibMagneticDataset {
    /// UNI 编号 (1–1651)
    pub uni_number: i32,
    /// 磁性空间群类型
    pub msg_type: i32,
    /// Hall 编号
    pub hall_number: i32,
    /// 张量秩
    pub tensor_rank: i32,
    /// 对称操作数
    pub n_operations: usize,
    /// 旋转矩阵
    pub rotations: Vec<Mat3I>,
    /// 平移矢量
    pub translations: Vec<Vec3>,
    /// 时间反演 (±1)
    pub time_reversals: Vec<i32>,
    /// 原子数
    pub n_atoms: usize,
    /// 对等原子
    pub equivalent_atoms: Vec<i32>,
    /// 标准晶胞原子数
    pub n_std_atoms: usize,
    /// 标准晶胞类型
    pub std_types: Vec<i32>,
    /// 标准晶胞位置
    pub std_positions: Vec<Vec3>,
    /// 标准晶胞张量
    pub std_tensors: Vec<f64>,
    /// 原点平移
    pub origin_shift: Vec3,
    /// 变换矩阵
    pub transformation_matrix: Mat3,
    /// 标准晶胞晶格
    pub std_lattice: Mat3,
    /// 原胞晶格
    pub primitive_lattice: Mat3,
    /// 标准旋转矩阵
    pub std_rotation_matrix: Mat3,
}

/// 磁性空间群类型（从数据库查询）。
#[derive(Debug, Clone)]
pub struct SpglibMagneticSpacegroupType {
    /// UNI 编号
    pub uni_number: i32,
    /// Litvin 编号
    pub litvin_number: i32,
    /// BNS 符号 (最多 8 字符)
    pub bns_number: String,
    /// OG 符号 (最多 12 字符)
    pub og_number: String,
    /// 晶体学编号 (1–230)
    pub number: i32,
    /// 磁性类型 (1-4)
    pub type_: i32,
}

// ========================================================================
// Public API
// ========================================================================

// ---------------------------------------------------------------------------
// Version
// ---------------------------------------------------------------------------

/// 获取 spglib 版本字符串。
pub fn spg_get_version() -> &'static str {
    SPGLIB_VERSION
}

/// 获取完整版本字符串。
pub fn spg_get_version_full() -> &'static str {
    SPGLIB_VERSION_FULL
}

/// 获取 Git 提交哈希。
pub fn spg_get_commit() -> &'static str {
    SPGLIB_COMMIT
}

/// 获取主版本号。
pub fn spg_get_major_version() -> i32 {
    SPGLIB_MAJOR_VERSION
}

/// 获取次版本号。
pub fn spg_get_minor_version() -> i32 {
    SPGLIB_MINOR_VERSION
}

/// 获取补丁版本号。
pub fn spg_get_micro_version() -> i32 {
    SPGLIB_MICRO_VERSION
}

// ---------------------------------------------------------------------------
// Error
// ---------------------------------------------------------------------------

/// 获取错误码对应的消息。
pub fn spg_get_error_message(error: SpglibError) -> &'static str {
    match error {
        SpglibError::Success => "no error",
        SpglibError::SpacegroupSearchFailed => "spacegroup search failed",
        SpglibError::CellStandardizationFailed => "cell standardization failed",
        SpglibError::SymmetryOperationSearchFailed => "symmetry operation search failed",
        SpglibError::AtomsTooClose => "too close distance between atoms",
        SpglibError::PointgroupNotFound => "pointgroup not found",
        SpglibError::NiggliFailed => "Niggli reduction failed",
        SpglibError::DelaunayFailed => "Delaunay reduction failed",
        SpglibError::ArraySizeShortage => "array size shortage",
    }
}

// ---------------------------------------------------------------------------
// Dataset (核心 API)
// ---------------------------------------------------------------------------

/// 获取空间群数据集。
///
/// 这是最常用的入口，返回完整的空间群信息包括标准晶胞、
/// 对称操作、Wyckoff 位置等。
pub fn spg_get_dataset(
    lattice: &Mat3,
    position: &[Vec3],
    types: &[i32],
    symprec: f64,
) -> Result<SpglibDataset, SpglibError> {
    get_dataset(lattice, position, types, None, 0, symprec, -1.0)
}

/// 获取空间群数据集（带角度容差）。
pub fn spgat_get_dataset(
    lattice: &Mat3,
    position: &[Vec3],
    types: &[i32],
    symprec: f64,
    angle_tolerance: f64,
) -> Result<SpglibDataset, SpglibError> {
    get_dataset(lattice, position, types, None, 0, symprec, angle_tolerance)
}

/// 获取空间群数据集（指定 Hall 编号）。
pub fn spg_get_dataset_with_hall_number(
    lattice: &Mat3,
    position: &[Vec3],
    types: &[i32],
    hall_number: i32,
    symprec: f64,
) -> Result<SpglibDataset, SpglibError> {
    get_dataset(lattice, position, types, None, hall_number, symprec, -1.0)
}

/// 获取空间群数据集（指定 Hall 编号，带角度容差）。
pub fn spgat_get_dataset_with_hall_number(
    lattice: &Mat3,
    position: &[Vec3],
    types: &[i32],
    hall_number: i32,
    symprec: f64,
    angle_tolerance: f64,
) -> Result<SpglibDataset, SpglibError> {
    get_dataset(lattice, position, types, None, hall_number, symprec, angle_tolerance)
}

/// 获取层状空间群数据集。
///
/// `aperiodic_axis` 指定无周期性方向的轴 (0, 1, 2)。
pub fn spg_get_layer_dataset(
    lattice: &Mat3,
    position: &[Vec3],
    types: &[i32],
    aperiodic_axis: i32,
    symprec: f64,
) -> Result<SpglibDataset, SpglibError> {
    use crate::cell::AperiodicAxis;
    let ap = match aperiodic_axis {
        0 => Some(AperiodicAxis::X),
        1 => Some(AperiodicAxis::Y),
        2 => Some(AperiodicAxis::Z),
        _ => return Err(SpglibError::SpacegroupSearchFailed),
    };
    get_dataset(lattice, position, types, ap, 0, symprec, -1.0)
}

// ---------------------------------------------------------------------------
// Symmetry operations
// ---------------------------------------------------------------------------

/// 获取对称操作（旋转和平移）。
///
/// 返回 `Some(Symmetry)` 包含旋转矩阵和对应的平移矢量。
/// 返回 `None` 表示搜索失败。
pub fn spg_get_symmetry(
    lattice: &Mat3,
    position: &[Vec3],
    types: &[i32],
    symprec: f64,
) -> Result<Symmetry, SpglibError> {
    get_symmetry_from_dataset(lattice, position, types, symprec, -1.0)
}

/// 获取对称操作（带角度容差）。
pub fn spgat_get_symmetry(
    lattice: &Mat3,
    position: &[Vec3],
    types: &[i32],
    symprec: f64,
    angle_tolerance: f64,
) -> Result<Symmetry, SpglibError> {
    get_symmetry_from_dataset(lattice, position, types, symprec, angle_tolerance)
}

/// 从空间群数据库获取对称操作。
///
/// 根据 Hall 编号直接返回所有空间群操作。
pub fn spg_get_symmetry_from_database(hall_number: i32) -> Option<Symmetry> {
    spgdb_get_spacegroup_operations(hall_number)
}

/// 从对称操作确定 Hall 编号。
///
/// 给定一组旋转和平移操作，搜索匹配的空间群 Hall 编号。
/// 返回 `None` 表示未找到匹配。
pub fn spg_get_hall_number_from_symmetry(
    rotations: &[Mat3I],
    translations: &[Vec3],
    symprec: f64,
) -> Result<i32, SpglibError> {
    let lattice: Mat3 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let hall_number = get_hall_number_from_symmetry(
        rotations, translations, &lattice, false, symprec,
    )?;
    if hall_number > 0 {
        Ok(hall_number)
    } else {
        Err(SpglibError::SpacegroupSearchFailed)
    }
}

/// 从对称操作确定空间群类型。
pub fn spg_get_spacegroup_type_from_symmetry(
    rotations: &[Mat3I],
    translations: &[Vec3],
    lattice: &Mat3,
    symprec: f64,
) -> Result<SpglibSpacegroupType, SpglibError> {
    let hall_number = get_hall_number_from_symmetry(
        rotations, translations, lattice, true, symprec,
    )?;
    if hall_number > 0 {
        get_spacegroup_type(hall_number)
    } else {
        Err(SpglibError::SpacegroupSearchFailed)
    }
}

// ---------------------------------------------------------------------------
// Standardization / refinement
// ---------------------------------------------------------------------------

/// 标准化晶胞。
///
/// 返回理想化的标准晶胞。参数：
/// - `to_primitive`: 若为 true，返回原胞而非常规晶胞
/// - `no_idealize`: 若为 true，跳过原子位置理想化
/// 返回 `None` 表示失败。
pub fn spg_standardize_cell(
    lattice: &Mat3,
    position: &[Vec3],
    types: &[i32],
    to_primitive: bool,
    no_idealize: bool,
    symprec: f64,
) -> Result<Cell, SpglibError> {
    spgat_standardize_cell(lattice, position, types, to_primitive, no_idealize, symprec, -1.0)
}

/// 标准化晶胞（带角度容差）。
pub fn spgat_standardize_cell(
    lattice: &Mat3,
    position: &[Vec3],
    types: &[i32],
    to_primitive: bool,
    no_idealize: bool,
    symprec: f64,
    angle_tolerance: f64,
) -> Result<Cell, SpglibError> {
    if to_primitive {
        if no_idealize {
            get_standardized_cell(lattice, position, types, true, symprec, angle_tolerance)
        } else {
            standardize_primitive(lattice, position, types, symprec, angle_tolerance)
        }
    } else {
        if no_idealize {
            get_standardized_cell(lattice, position, types, false, symprec, angle_tolerance)
        } else {
            standardize_cell(lattice, position, types, symprec, angle_tolerance)
        }
    }
}

/// 寻找原胞。
///
/// 将任意晶胞约化为其原胞，返回原胞结构和新的原子数。
pub fn spg_find_primitive(
    lattice: &Mat3,
    position: &[Vec3],
    types: &[i32],
    symprec: f64,
) -> Result<Cell, SpglibError> {
    standardize_primitive(lattice, position, types, symprec, -1.0)
}

/// 寻找原胞（带角度容差）。
pub fn spgat_find_primitive(
    lattice: &Mat3,
    position: &[Vec3],
    types: &[i32],
    symprec: f64,
    angle_tolerance: f64,
) -> Result<Cell, SpglibError> {
    standardize_primitive(lattice, position, types, symprec, angle_tolerance)
}

/// 精细化晶胞。
///
/// 对输入晶胞进行理想化处理，返回标准化的常规晶胞。
pub fn spg_refine_cell(
    lattice: &Mat3,
    position: &[Vec3],
    types: &[i32],
    symprec: f64,
) -> Result<Cell, SpglibError> {
    standardize_cell(lattice, position, types, symprec, -1.0)
}

/// 精细化晶胞（带角度容差）。
pub fn spgat_refine_cell(
    lattice: &Mat3,
    position: &[Vec3],
    types: &[i32],
    symprec: f64,
    angle_tolerance: f64,
) -> Result<Cell, SpglibError> {
    standardize_cell(lattice, position, types, symprec, angle_tolerance)
}

// ---------------------------------------------------------------------------
// Information retrieval
// ---------------------------------------------------------------------------

/// 获取空间群的国际符号。
///
/// 返回 `Some((spacegroup_number, international_symbol))`。
pub fn spg_get_international(
    lattice: &Mat3,
    position: &[Vec3],
    types: &[i32],
    symprec: f64,
) -> Result<(i32, String), SpglibError> {
    get_international(lattice, position, types, symprec, -1.0)
}

/// 获取空间群的国际符号（带角度容差）。
pub fn spgat_get_international(
    lattice: &Mat3,
    position: &[Vec3],
    types: &[i32],
    symprec: f64,
    angle_tolerance: f64,
) -> Result<(i32, String), SpglibError> {
    get_international(lattice, position, types, symprec, angle_tolerance)
}

/// 获取空间群的 Schoenflies 符号。
///
/// 返回 `Some((spacegroup_number, schoenflies_symbol))`。
pub fn spg_get_schoenflies(
    lattice: &Mat3,
    position: &[Vec3],
    types: &[i32],
    symprec: f64,
) -> Result<(i32, String), SpglibError> {
    get_schoenflies(lattice, position, types, symprec, -1.0)
}

/// 获取空间群的 Schoenflies 符号（带角度容差）。
pub fn spgat_get_schoenflies(
    lattice: &Mat3,
    position: &[Vec3],
    types: &[i32],
    symprec: f64,
    angle_tolerance: f64,
) -> Result<(i32, String), SpglibError> {
    get_schoenflies(lattice, position, types, symprec, angle_tolerance)
}

/// 获取对称操作的多重数（即对称操作的个数）。
///
/// 返回 `None` 表示搜索失败。
pub fn spg_get_multiplicity(
    lattice: &Mat3,
    position: &[Vec3],
    types: &[i32],
    symprec: f64,
) -> Result<usize, SpglibError> {
    get_multiplicity(lattice, position, types, symprec, -1.0)
}

/// 获取对称操作的多重数（带角度容差）。
pub fn spgat_get_multiplicity(
    lattice: &Mat3,
    position: &[Vec3],
    types: &[i32],
    symprec: f64,
    angle_tolerance: f64,
) -> Result<usize, SpglibError> {
    get_multiplicity(lattice, position, types, symprec, angle_tolerance)
}

/// 根据 Hall 编号获取空间群类型信息。
pub fn spg_get_spacegroup_type(hall_number: i32) -> Result<SpglibSpacegroupType, SpglibError> {
    if hall_number > 0 && hall_number < 531 {
        get_spacegroup_type(hall_number)
    } else {
        Err(SpglibError::SpacegroupSearchFailed)
    }
}

/// 获取点群信息。
///
/// 给定一组旋转操作，确定对应的晶体学点群。
/// 返回 `Some((symbol, transform_matrix, pointgroup_number))`。
/// symbol 最多 6 字符。
pub fn spg_get_pointgroup(
    rotations: &[Mat3I],
) -> Result<(String, Mat3I, i32), SpglibError> {
    let mut transform_mat = [[0; 3]; 3];
    let pointgroup = ptg_get_transformation_matrix(
        &mut transform_mat, rotations, None, /* aperiodic_axis */
    );

    if pointgroup.number == 0 {
        return Err(SpglibError::PointgroupNotFound);
    }

    Ok((pointgroup.symbol.to_string(), transform_mat, pointgroup.number))
}

// ---------------------------------------------------------------------------
// Magnetic space groups
// ---------------------------------------------------------------------------

/// 获取磁性空间群类型。
///
/// 根据 UNI 编号查询磁性空间群类型信息。返回默认值（全零）表示未找到。
pub fn spg_get_magnetic_spacegroup_type(
    uni_number: i32,
) -> SpglibMagneticSpacegroupType {
    let msgtype = crate::msg_database::msgdb_get_magnetic_spacegroup_type(uni_number);
    SpglibMagneticSpacegroupType {
        uni_number: msgtype.uni_number,
        litvin_number: msgtype.litvin_number,
        bns_number: msgtype.bns_number.to_string(),
        og_number: msgtype.og_number.to_string(),
        number: msgtype.number,
        type_: msgtype.type_,
    }
}

/// 获取磁性空间群类型（从对称操作）。
///
/// `time_reversals` 为 `None` 时全部视为 0（无时间反演）。
pub fn spg_get_magnetic_spacegroup_type_from_symmetry(
    rotations: &[Mat3I],
    translations: &[Vec3],
    time_reversals: Option<&[i32]>,
    lattice: &Mat3,
    symprec: f64,
) -> SpglibMagneticSpacegroupType {
    let n_ops = rotations.len();
    let mut mag_sym = crate::symmetry::MagneticSymmetry::new(n_ops);
    for i in 0..n_ops {
        mag_sym.rot[i] = rotations[i];
        mag_sym.trans[i] = translations[i];
        mag_sym.timerev[i] = time_reversals.map_or(0, |tr| tr[i]);
    }

    match crate::magnetic_spacegroup::msg_identify_magnetic_space_group_type(
        lattice, &mag_sym, symprec,
    ) {
        Some(dataset) => spg_get_magnetic_spacegroup_type(dataset.uni_number),
        None => SpglibMagneticSpacegroupType {
            uni_number: 0,
            litvin_number: 0,
            bns_number: String::new(),
            og_number: String::new(),
            number: 0,
            type_: 0,
        },
    }
}

// ---------------------------------------------------------------------------
// Lattice reduction
// ---------------------------------------------------------------------------

/// Delaunay 晶格约化。
///
pub fn spg_delaunay_reduce(lattice: &mut Mat3, symprec: f64) -> Result<(), SpglibError> {
    let mut red_lattice = [[0.0; 3]; 3];
    if del_delaunay_reduce(&mut red_lattice, lattice, symprec) {
        mat_copy_matrix_d3(lattice, &red_lattice);
        Ok(())
    } else {
        Err(SpglibError::DelaunayFailed)
    }
}

pub fn spg_niggli_reduce(lattice: &mut Mat3, symprec: f64) -> Result<(), SpglibError> {
    if niggli_reduce(lattice, symprec, None) {
        Ok(())
    } else {
        Err(SpglibError::NiggliFailed)
    }
}

// ---------------------------------------------------------------------------
// K-point grid
// ---------------------------------------------------------------------------

/// 从网格地址获取网格点索引。
pub fn spg_get_grid_point_from_address(grid_address: &[i32; 3], mesh: &[i32; 3]) -> usize {
    let mut address_double = [0i32; 3];
    let is_shift = [0i32; 3];
    crate::kgrid::kgd_get_grid_address_double_mesh(&mut address_double, grid_address, mesh, &is_shift);
    crate::kgrid::kgd_get_dense_grid_point_double_mesh(&address_double, mesh)
}

/// 获取不可约倒易网格。
///
/// 返回不可约网格点的数量。`grid_address` 和 `ir_mapping_table`
/// 需预分配足够空间（通常为 `mesh[0]*mesh[1]*mesh[2]`）。
pub fn spg_get_ir_reciprocal_mesh(
    grid_address: &mut [[i32; 3]],
    ir_mapping_table: &mut [usize],
    mesh: &[i32; 3],
    is_shift: &[i32; 3],
    is_time_reversal: bool,
    lattice: &Mat3,
    position: &[Vec3],
    types: &[i32],
    symprec: f64,
) -> Result<usize, SpglibError> {
    get_ir_reciprocal_mesh(
        grid_address, ir_mapping_table, mesh, is_shift,
        is_time_reversal, lattice, position, types, symprec, -1.0,
    )
}

/// 获取不可约倒易网格（密集版本，使用 usize 映射表）。
pub fn spg_get_dense_ir_reciprocal_mesh(
    grid_address: &mut [[i32; 3]],
    ir_mapping_table: &mut [usize],
    mesh: &[i32; 3],
    is_shift: &[i32; 3],
    is_time_reversal: bool,
    lattice: &Mat3,
    position: &[Vec3],
    types: &[i32],
    symprec: f64,
) -> Result<usize, SpglibError> {
    get_dense_ir_reciprocal_mesh(
        grid_address, ir_mapping_table, mesh, is_shift,
        is_time_reversal, lattice, position, types, symprec, -1.0,
    )
}

/// 获取稳定化倒易网格（给定对称操作和 q 点）。
pub fn spg_get_stabilized_reciprocal_mesh(
    grid_address: &mut [[i32; 3]],
    ir_mapping_table: &mut [usize],
    mesh: &[i32; 3],
    is_shift: &[i32; 3],
    is_time_reversal: bool,
    rotations: &[Mat3I],
    qpoints: &[[f64; 3]],
) -> usize {
    use crate::mathfunc::MatINT;
    let mut rot_real = MatINT::new(rotations.len());
    for (i, r) in rotations.iter().enumerate() {
        rot_real.mat[i] = *r;
    }
    crate::kpoint::kpt_get_stabilized_reciprocal_mesh(
        grid_address, ir_mapping_table, mesh, is_shift,
        if is_time_reversal { 1 } else { 0 },
        &rot_real, qpoints,
    )
}

/// 通过旋转矩阵获取密集网格点。
pub fn spg_get_dense_grid_points_by_rotations(
    rot_grid_points: &mut [usize],
    address_orig: &[i32; 3],
    rot_reciprocal: &[Mat3I],
    mesh: &[i32; 3],
    is_shift: &[i32; 3],
) {
    use crate::mathfunc::MatINT;
    let mut rot = MatINT::new(rot_reciprocal.len());
    for (i, r) in rot_reciprocal.iter().enumerate() {
        rot.mat[i] = *r;
    }
    crate::kpoint::kpt_get_dense_grid_points_by_rotations(
        rot_grid_points, address_orig, &rot, mesh, is_shift,
    )
}

/// 通过旋转矩阵获取 BZ 网格点。
pub fn spg_get_dense_BZ_grid_points_by_rotations(
    rot_grid_points: &mut [usize],
    address_orig: &[i32; 3],
    rot_reciprocal: &[Mat3I],
    mesh: &[i32; 3],
    is_shift: &[i32; 3],
    bz_map: &[usize],
) {
    use crate::mathfunc::MatINT;
    let mut rot = MatINT::new(rot_reciprocal.len());
    for (i, r) in rot_reciprocal.iter().enumerate() {
        rot.mat[i] = *r;
    }
    crate::kpoint::kpt_get_dense_BZ_grid_points_by_rotations(
        rot_grid_points, address_orig, &rot, mesh, is_shift, bz_map,
    )
}

/// 将网格点重新定位到第一布里渊区。
///
/// 返回 BZ 网格点的数量。`bz_map` 中未映射的条目设为 `usize::MAX`（对应 C 中的 -1）。
pub fn spg_relocate_BZ_grid_address(
    bz_grid_address: &mut [[i32; 3]],
    bz_map: &mut [usize],
    grid_address: &[[i32; 3]],
    mesh: &[i32; 3],
    rec_lattice: &Mat3,
    is_shift: &[i32; 3],
) -> usize {
    crate::kpoint::kpt_relocate_bz_grid_address(
        bz_grid_address, bz_map, grid_address, mesh, rec_lattice, is_shift,
    )
}

// ========================================================================
// Internal functions
// ========================================================================

/// 内部：核心数据集获取逻辑。
fn get_dataset(
    lattice: &Mat3,
    position: &[Vec3],
    types: &[i32],
    aperiodic_axis: Option<AperiodicAxis>,
    hall_number: i32,
    symprec: f64,
    angle_tolerance: f64,
) -> Result<SpglibDataset, SpglibError> {
    if hall_number > 530 {
        return Err(SpglibError::SpacegroupSearchFailed);
    }

    let num_atom = position.len();
    let mut cell = Cell::new(num_atom, TensorRank::NoSpin);
    if aperiodic_axis.is_none() {
        cel_set_cell(&mut cell, lattice, position, types);
        if cel_any_overlap_with_same_type(&cell, symprec) {
            return Err(SpglibError::AtomsTooClose);
        }
    } else {
        cel_set_layer_cell(&mut cell, lattice, position, types, aperiodic_axis);
        if cel_layer_any_overlap_with_same_type(&cell, aperiodic_axis.unwrap(), symprec) {
            return Err(SpglibError::AtomsTooClose);
        }
    }

    let container = det_determine_all(&cell, hall_number, symprec, angle_tolerance)
        .ok_or(SpglibError::SpacegroupSearchFailed)?;

    let spacegroup = container.spacegroup.as_ref()
        .ok_or(SpglibError::SpacegroupSearchFailed)?;
    let primitive = container.primitive.as_ref()
        .ok_or(SpglibError::SpacegroupSearchFailed)?;
    let exstr = container.exact_structure.as_ref()
        .ok_or(SpglibError::SpacegroupSearchFailed)?;

    let dataset = set_dataset(&cell, primitive, spacegroup, exstr)
        .ok_or(SpglibError::SpacegroupSearchFailed)?;
    Ok(dataset)
}

/// 将 Cell 数据设置到输入晶胞。
fn cel_set_cell(cell: &mut Cell, lattice: &Mat3, position: &[Vec3], types: &[i32]) {
    mat_copy_matrix_d3(&mut cell.lattice, lattice);
    for i in 0..cell.size {
        cell.types[i] = types[i];
        mat_copy_vector_d3(&mut cell.position[i], &position[i]);
    }
}

/// 将层状 Cell 数据设置到输入晶胞。
fn cel_set_layer_cell(
    cell: &mut Cell,
    lattice: &Mat3,
    position: &[Vec3],
    types: &[i32],
    aperiodic_axis: Option<AperiodicAxis>,
) {
    mat_copy_matrix_d3(&mut cell.lattice, lattice);
    cell.aperiodic_axis = aperiodic_axis;
    for i in 0..cell.size {
        cell.types[i] = types[i];
        mat_copy_vector_d3(&mut cell.position[i], &position[i]);
    }
}

/// 从内部结构填充 SpglibDataset。
fn set_dataset(
    cell: &Cell,
    primitive: &Primitive,
    spacegroup: &Spacegroup,
    exstr: &crate::refinement::ExactStructure,
) -> Option<SpglibDataset> {
    let n_atoms = cell.size;
    let n_operations = exstr.symmetry.size;

    let mut dataset = SpglibDataset {
        spacegroup_number: spacegroup.number,
        hall_number: spacegroup.hall_number,
        international_symbol: spacegroup.international_short.clone(),
        hall_symbol: spacegroup.hall_symbol.clone(),
        choice: spacegroup.choice.clone(),
        transformation_matrix: [[0.0; 3]; 3],
        origin_shift: spacegroup.origin_shift,
        n_operations,
        rotations: vec![[[0; 3]; 3]; n_operations],
        translations: vec![[0.0; 3]; n_operations],
        n_atoms,
        wyckoffs: vec![0i32; n_atoms],
        site_symmetry_symbols: vec![String::new(); n_atoms],
        equivalent_atoms: vec![0i32; n_atoms],
        crystallographic_orbits: vec![0i32; n_atoms],
        mapping_to_primitive: vec![0i32; n_atoms],
        n_std_atoms: exstr.bravais.size,
        std_lattice: exstr.bravais.lattice,
        std_positions: exstr.bravais.position.clone(),
        std_types: exstr.bravais.types.clone(),
        std_rotation_matrix: [[0.0; 3]; 3],
        std_mapping_to_primitive: vec![0i32; exstr.bravais.size],
        primitive_lattice: [[0.0; 3]; 3],
        pointgroup_symbol: String::new(),
    };

    // Transformation matrix: inv(brv_lat) * cell_lat
    let inv_lat = mat_inverse_matrix_d3(&spacegroup.bravais_lattice, 0.0)?;
    dataset.transformation_matrix = mat_multiply_matrix_d3(&inv_lat, &cell.lattice);

    // Copy symmetry operations
    for i in 0..n_operations {
        dataset.rotations[i] = exstr.symmetry.rot[i];
        dataset.translations[i] = exstr.symmetry.trans[i];
    }

    // Copy Wyckoff, site symmetry, equivalent atoms, crystallographic orbits
    for i in 0..n_atoms {
        dataset.wyckoffs[i] = exstr.wyckoffs[i];
        dataset.site_symmetry_symbols[i] = exstr.site_symmetry_symbols[i].clone();
        dataset.equivalent_atoms[i] = exstr.equivalent_atoms[i];
        dataset.crystallographic_orbits[i] = exstr.crystallographic_orbits[i];
    }

    // Mapping to primitive
    if let Some(prim_cell) = &primitive.cell {
        dataset.primitive_lattice = prim_cell.lattice;
    }
    for i in 0..n_atoms {
        dataset.mapping_to_primitive[i] = primitive.mapping_table[i];
    }

    // Standardized cell data
    for i in 0..dataset.n_std_atoms {
        dataset.std_mapping_to_primitive[i] = exstr.std_mapping_to_primitive[i];
    }

    // Standard rotation matrix
    dataset.std_rotation_matrix = exstr.rotation;

    // Point group symbol
    let pointgroup = ptg_get_pointgroup(spacegroup.pointgroup_number);
    dataset.pointgroup_symbol = pointgroup.symbol.to_string();

    Some(dataset)
}

/// 从数据集获取对称操作。
fn get_symmetry_from_dataset(
    lattice: &Mat3,
    position: &[Vec3],
    types: &[i32],
    symprec: f64,
    angle_tolerance: f64,
) -> Result<Symmetry, SpglibError> {
    let dataset = get_dataset(lattice, position, types, None, 0, symprec, angle_tolerance)?;
    let n_ops = dataset.n_operations;
    let mut sym = Symmetry::new(n_ops);
    for i in 0..n_ops {
        sym.rot[i] = dataset.rotations[i];
        sym.trans[i] = dataset.translations[i];
    }
    Ok(sym)
}

/// 获取多重数。
fn get_multiplicity(
    lattice: &Mat3,
    position: &[Vec3],
    types: &[i32],
    symprec: f64,
    angle_tolerance: f64,
) -> Result<usize, SpglibError> {
    let dataset = get_dataset(lattice, position, types, None, 0, symprec, angle_tolerance)?;
    Ok(dataset.n_operations)
}

/// 寻找原胞。
fn standardize_primitive(
    lattice: &Mat3,
    position: &[Vec3],
    types: &[i32],
    symprec: f64,
    angle_tolerance: f64,
) -> Result<Cell, SpglibError> {
    let dataset = get_dataset(lattice, position, types, None, 0, symprec, angle_tolerance)?;

    let centering = get_centering(dataset.hall_number)
        .ok_or(SpglibError::CellStandardizationFailed)?;

    let mut bravais = Cell::new(dataset.n_std_atoms, TensorRank::NoSpin);
    mat_copy_matrix_d3(&mut bravais.lattice, &dataset.std_lattice);
    for i in 0..dataset.n_std_atoms {
        bravais.types[i] = dataset.std_types[i];
        mat_copy_vector_d3(&mut bravais.position[i], &dataset.std_positions[i]);
    }

    let mut mapping_table: Vec<usize> = vec![0; bravais.size];
    let identity: Mat3 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

    let primitive = spa_transform_to_primitive(&mut mapping_table, &bravais, &identity, centering, symprec)
        .ok_or(SpglibError::CellStandardizationFailed)?;

    // Validation: mapping_table should be identity for standardize_primitive
    for i in 0..primitive.size {
        if mapping_table[i] != i {
            debug::warning_print(format_args!(
                "spglib: spa_transform_to_primitive failed ({} != {})\n",
                mapping_table[i], i
            ));
            return Err(SpglibError::CellStandardizationFailed);
        }
    }

    Ok(primitive)
}

/// 标准化晶胞。
fn standardize_cell(
    lattice: &Mat3,
    position: &[Vec3],
    types: &[i32],
    symprec: f64,
    angle_tolerance: f64,
) -> Result<Cell, SpglibError> {
    let dataset = get_dataset(lattice, position, types, None, 0, symprec, angle_tolerance)?;

    let n_std = dataset.n_std_atoms;
    let mut cell = Cell::new(n_std, TensorRank::NoSpin);
    mat_copy_matrix_d3(&mut cell.lattice, &dataset.std_lattice);
    for i in 0..n_std {
        cell.types[i] = dataset.std_types[i];
        mat_copy_vector_d3(&mut cell.position[i], &dataset.std_positions[i]);
    }
    Ok(cell)
}

/// 获取标准化晶胞（无理想化）。
fn get_standardized_cell(
    lattice: &Mat3,
    position: &[Vec3],
    types: &[i32],
    to_primitive: bool,
    symprec: f64,
    angle_tolerance: f64,
) -> Result<Cell, SpglibError> {
    let dataset = get_dataset(lattice, position, types, None, 0, symprec, angle_tolerance)?;
    let centering = get_centering(dataset.hall_number)
        .ok_or(SpglibError::CellStandardizationFailed)?;

    let num_atom = position.len();
    let mut cell = Cell::new(num_atom, TensorRank::NoSpin);
    mat_copy_matrix_d3(&mut cell.lattice, lattice);
    for i in 0..num_atom {
        cell.types[i] = types[i];
        mat_copy_vector_d3(&mut cell.position[i], &position[i]);
    }

    let mut mapping_table: Vec<usize> = vec![0; num_atom];
    let primitive = spa_transform_to_primitive(
        &mut mapping_table, &cell, &dataset.transformation_matrix, centering, symprec,
    ).ok_or(SpglibError::CellStandardizationFailed)?;

    // Validate mapping
    for i in 0..num_atom {
        if mapping_table[i] != dataset.mapping_to_primitive[i] as usize {
            debug::warning_print(format_args!(
                "spglib: spa_transform_to_primitive failed ({} != {})\n",
                mapping_table[i], dataset.mapping_to_primitive[i]
            ));
            return Err(SpglibError::CellStandardizationFailed);
        }
    }

    if to_primitive || matches!(centering, Centering::Primitive) {
        return Ok(primitive);
    }

    let std_cell = spa_transform_from_primitive(&primitive, centering, symprec)
        .ok_or(SpglibError::CellStandardizationFailed)?;
    Ok(std_cell)
}

/// 获取国际符号。
fn get_international(
    lattice: &Mat3,
    position: &[Vec3],
    types: &[i32],
    symprec: f64,
    angle_tolerance: f64,
) -> Result<(i32, String), SpglibError> {
    let dataset = get_dataset(lattice, position, types, None, 0, symprec, angle_tolerance)?;
    if dataset.spacegroup_number > 0 {
        Ok((dataset.spacegroup_number, dataset.international_symbol))
    } else {
        Err(SpglibError::SpacegroupSearchFailed)
    }
}

/// 获取 Schoenflies 符号。
fn get_schoenflies(
    lattice: &Mat3,
    position: &[Vec3],
    types: &[i32],
    symprec: f64,
    angle_tolerance: f64,
) -> Result<(i32, String), SpglibError> {
    let dataset = get_dataset(lattice, position, types, None, 0, symprec, angle_tolerance)?;
    if dataset.spacegroup_number > 0 {
        if let Ok(spgtype) = get_spacegroup_type(dataset.hall_number) {
            return Ok((dataset.spacegroup_number, spgtype.schoenflies));
        }
    }
    Err(SpglibError::SpacegroupSearchFailed)
}

/// 获取 Hall 编号对应的 Centering。
fn get_centering(hall_number: i32) -> Option<Centering> {
    Some(spgdb_get_spacegroup_type(hall_number).centering)
}

/// 从对称操作获取 Hall 编号。
fn get_hall_number_from_symmetry(
    rotations: &[Mat3I],
    translations: &[Vec3],
    lattice: &Mat3,
    transform_lattice_by_tmat: bool,
    symprec: f64,
) -> Result<i32, SpglibError> {
    let num_ops = rotations.len();
    let mut symmetry = Symmetry::new(num_ops);
    for i in 0..num_ops {
        symmetry.rot[i] = rotations[i];
        symmetry.trans[i] = translations[i];
    }

    let mut t_mat = [[0.0; 3]; 3];
    let prim_sym = prm_get_primitive_symmetry(&mut t_mat, &symmetry, symprec)
        .ok_or(SpglibError::SpacegroupSearchFailed)?;

    let prim_lat = if transform_lattice_by_tmat {
        let t_mat_inv = mat_inverse_matrix_d3(&t_mat, symprec)
            .ok_or(SpglibError::SpacegroupSearchFailed)?;
        mat_multiply_matrix_d3(lattice, &t_mat_inv)
    } else {
        *lattice
    };

    let spacegroup = spa_search_spacegroup_with_symmetry(&prim_sym, &prim_lat, symprec)
        .ok_or(SpglibError::SpacegroupSearchFailed)?;
    Ok(spacegroup.hall_number)
}

/// 获取 SpglibSpacegroupType。
fn get_spacegroup_type(hall_number: i32) -> Result<SpglibSpacegroupType, SpglibError> {
    if hall_number <= 0 || hall_number >= 531 {
        return Err(SpglibError::SpacegroupSearchFailed);
    }

    let spgtype = spgdb_get_spacegroup_type(hall_number);
    let pointgroup = ptg_get_pointgroup(spgtype.pointgroup_number);

    Ok(SpglibSpacegroupType {
        number: spgtype.number,
        hall_number,
        schoenflies: spgtype.schoenflies,
        hall_symbol: spgtype.hall_symbol,
        choice: spgtype.choice,
        international: spgtype.international,
        international_full: spgtype.international_full,
        international_short: spgtype.international_short,
        pointgroup_international: pointgroup.symbol.to_string(),
        pointgroup_schoenflies: pointgroup.schoenflies.to_string(),
        arithmetic_crystal_class_number: 0, // TODO: arth_get_symbol
        arithmetic_crystal_class_symbol: String::new(),
    })
}

/// 内部：不可约网格计算。
fn get_ir_reciprocal_mesh(
    grid_address: &mut [[i32; 3]],
    ir_mapping_table: &mut [usize],
    mesh: &[i32; 3],
    is_shift: &[i32; 3],
    is_time_reversal: bool,
    lattice: &Mat3,
    position: &[Vec3],
    types: &[i32],
    symprec: f64,
    angle_tolerance: f64,
) -> Result<usize, SpglibError> {
    let dataset = get_dataset(lattice, position, types, None, 0, symprec, angle_tolerance)?;

    use crate::mathfunc::{MatINT, mat_copy_matrix_i3};
    let mut rotations = MatINT::new(dataset.n_operations);
    for i in 0..dataset.n_operations {
        mat_copy_matrix_i3(&mut rotations.mat[i], &dataset.rotations[i]);
    }
    let rot_reciprocal = crate::kpoint::kpt_get_point_group_reciprocal(
        &rotations,
        if is_time_reversal { 1 } else { 0 },
    ).ok_or(SpglibError::SpacegroupSearchFailed)?;
    let num_ir = crate::kpoint::kpt_get_irreducible_reciprocal_mesh(
        grid_address, ir_mapping_table, mesh, is_shift, &rot_reciprocal,
    );
    Ok(num_ir)
}

/// 内部：密集不可约网格。
fn get_dense_ir_reciprocal_mesh(
    grid_address: &mut [[i32; 3]],
    ir_mapping_table: &mut [usize],
    mesh: &[i32; 3],
    is_shift: &[i32; 3],
    is_time_reversal: bool,
    lattice: &Mat3,
    position: &[Vec3],
    types: &[i32],
    symprec: f64,
    angle_tolerance: f64,
) -> Result<usize, SpglibError> {
    let dataset = get_dataset(lattice, position, types, None, 0, symprec, angle_tolerance)?;

    use crate::mathfunc::{MatINT, mat_copy_matrix_i3};
    let mut rotations = MatINT::new(dataset.n_operations);
    for i in 0..dataset.n_operations {
        mat_copy_matrix_i3(&mut rotations.mat[i], &dataset.rotations[i]);
    }
    let rot_reciprocal = crate::kpoint::kpt_get_point_group_reciprocal(
        &rotations,
        if is_time_reversal { 1 } else { 0 },
    ).ok_or(SpglibError::SpacegroupSearchFailed)?;
    let num_ir = crate::kpoint::kpt_get_dense_irreducible_reciprocal_mesh(
        grid_address, ir_mapping_table, mesh, is_shift, &rot_reciprocal,
    );
    Ok(num_ir)
}
