// spacegroup.rs

use crate::cell::{Cell, TensorRank, cel_trim_cell};
use crate::debug;
use crate::delaunay::del_layer_delaunay_reduce_2D;
use crate::delaunay::{del_delaunay_reduce, del_layer_delaunay_reduce};
use crate::hall_symbol::hal_match_hall_symbol_db;
use crate::mathfunc::{
    Mat3, Mat3I, Vec3, mat_cast_matrix_3d_to_3i, mat_cast_matrix_3i_to_3d,
    mat_check_identity_matrix_d3, mat_check_identity_matrix_id3, mat_copy_matrix_d3,
    mat_copy_vector_d3, mat_dabs, mat_get_determinant_d3, mat_get_determinant_i3, mat_get_metric,
    mat_get_similar_matrix_d3, mat_inverse_matrix_d3, mat_is_int_matrix, mat_multiply_matrix_d3,
    mat_multiply_matrix_di3, mat_multiply_matrix_id3, mat_multiply_matrix_vector_d3,
    mat_norm_squared_d3,
};
use crate::niggli::niggli_reduce;
use crate::pointgroup::{
    Holohedry, Laue, Pointgroup, ptg_get_pointsymmetry, ptg_get_transformation_matrix,
};
use crate::primitive::Primitive;
use crate::spg_database::{Centering, SpacegroupType, spgdb_get_spacegroup_type};
use crate::symmetry::{Symmetry, sym_get_operation, sym_reduce_operation};

const REDUCE_RATE: f64 = 0.95;
const NUM_ATTEMPT: i32 = 100;
const INT_PREC: f64 = 0.1;
const ZERO_PREC: f64 = 1e-10;

// --- Static Data (Matrices) ---

pub static M_BCC: Mat3I = [[0, 1, 1], [1, 0, 1], [1, 1, 0]];
pub static M_FCC: Mat3I = [[-1, 1, 1], [1, -1, 1], [1, 1, -1]];
pub static M_AC: Mat3I = [[1, 0, 0], [0, 1, 1], [0, -1, 1]];
pub static M_BC: Mat3I = [[1, 0, 1], [0, 1, 0], [-1, 0, 1]];
pub static M_CC: Mat3I = [[1, -1, 0], [1, 1, 0], [0, 0, 1]];
pub static M_RC: Mat3I = [[1, 0, 1], [-1, 1, 1], [0, -1, 1]];

static M_BCC_INV: Mat3 = [[-0.5, 0.5, 0.5], [0.5, -0.5, 0.5], [0.5, 0.5, -0.5]];
static M_FCC_INV: Mat3 = [[0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]];
static M_AC_INV: Mat3 = [[1.0, 0.0, 0.0], [0.0, 0.5, -0.5], [0.0, 0.5, 0.5]];
static M_BC_INV: Mat3 = [[0.5, 0.0, -0.5], [0.0, 1.0, 0.0], [0.5, 0.0, 0.5]];
static M_CC_INV: Mat3 = [[0.5, 0.5, 0.0], [-0.5, 0.5, 0.0], [0.0, 0.0, 1.0]];
static M_RC_INV: Mat3 = [
    [2. / 3., -1. / 3., -1. / 3.],
    [1. / 3., 1. / 3., -2. / 3.],
    [1. / 3., 1. / 3., 1. / 3.],
];

static IDENTITY: Mat3 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
static MONOCLI_I2C: Mat3 = [[1.0, 0.0, -1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]];
static MONOCLI_A2C: Mat3 = [[0.0, 0.0, 1.0], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0]];
static RHOMBO_OBVERSE: Mat3 = [
    [2. / 3., -1. / 3., -1. / 3.],
    [1. / 3., 1. / 3., -2. / 3.],
    [1. / 3., 1. / 3., 1. / 3.],
];
static RHOMBO_REVERSE: Mat3 = [
    [1. / 3., -2. / 3., 1. / 3.],
    [2. / 3., -1. / 3., -1. / 3.],
    [1. / 3., 1. / 3., 1. / 3.],
];
static A2C: Mat3 = [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
static B2C: Mat3 = [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]];

static A_MAT: Mat3 = [[1.0, 0.0, 0.0], [0.0, 0.5, -0.5], [0.0, 0.5, 0.5]];
static C_MAT: Mat3 = [[0.5, 0.5, 0.0], [-0.5, 0.5, 0.0], [0.0, 0.0, 1.0]];
static R_MAT: Mat3 = [
    [2. / 3., -1. / 3., -1. / 3.],
    [1. / 3., 1. / 3., -2. / 3.],
    [1. / 3., 1. / 3., 1. / 3.],
];
static I_MAT: Mat3 = [[-0.5, 0.5, 0.5], [0.5, -0.5, 0.5], [0.5, 0.5, -0.5]];
static F_MAT: Mat3 = [[0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]];

static CHANGE_OF_BASIS_MONOCLI: [Mat3; 48] = [
    [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
    [[0., 0., 1.], [0., -1., 0.], [1., 0., 0.]],
    [[0., 0., 1.], [1., 0., 0.], [0., 1., 0.]],
    [[1., 0., 0.], [0., 0., 1.], [0., -1., 0.]],
    [[0., 1., 0.], [0., 0., 1.], [1., 0., 0.]],
    [[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]],
    [[-1., 0., 1.], [0., 1., 0.], [-1., 0., 0.]],
    [[1., 0., -1.], [0., -1., 0.], [0., 0., -1.]],
    [[0., 1., -1.], [1., 0., 0.], [0., 0., -1.]],
    [[-1., -1., 0.], [0., 0., 1.], [-1., 0., 0.]],
    [[1., -1., 0.], [0., 0., 1.], [0., -1., 0.]],
    [[0., 1., 1.], [1., 0., 0.], [0., 1., 0.]],
    [[0., 0., -1.], [0., 1., 0.], [1., 0., -1.]],
    [[-1., 0., 0.], [0., -1., 0.], [-1., 0., 1.]],
    [[0., -1., 0.], [1., 0., 0.], [0., -1., 1.]],
    [[0., 1., 0.], [0., 0., 1.], [1., 1., 0.]],
    [[-1., 0., 0.], [0., 0., 1.], [-1., 1., 0.]],
    [[0., 0., -1.], [1., 0., 0.], [0., -1., -1.]],
    [[1., 0., 0.], [0., -1., 0.], [0., 0., -1.]],
    [[0., 0., -1.], [0., 1., 0.], [1., 0., 0.]],
    [[0., 0., 1.], [-1., 0., 0.], [0., -1., 0.]],
    [[-1., 0., 0.], [0., 0., -1.], [0., -1., 0.]],
    [[0., 1., 0.], [0., 0., -1.], [-1., 0., 0.]],
    [[0., 1., 0.], [-1., 0., 0.], [0., 0., 1.]],
    [[-1., 0., -1.], [0., -1., 0.], [-1., 0., 0.]],
    [[1., 0., 1.], [0., 1., 0.], [0., 0., 1.]],
    [[0., -1., -1.], [-1., 0., 0.], [0., 0., -1.]],
    [[1., -1., 0.], [0., 0., -1.], [1., 0., 0.]],
    [[-1., -1., 0.], [0., 0., -1.], [0., -1., 0.]],
    [[0., -1., 1.], [-1., 0., 0.], [0., -1., 0.]],
    [[0., 0., 1.], [0., -1., 0.], [1., 0., 1.]],
    [[-1., 0., 0.], [0., 1., 0.], [-1., 0., -1.]],
    [[0., 1., 0.], [-1., 0., 0.], [0., 1., 1.]],
    [[0., 1., 0.], [0., 0., -1.], [-1., 1., 0.]],
    [[1., 0., 0.], [0., 0., -1.], [1., 1., 0.]],
    [[0., 0., -1.], [-1., 0., 0.], [0., 1., -1.]],
    [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
    [[0., 1., 0.], [1., 0., 0.], [0., 0., -1.]],
    [[0., -1., 0.], [1., -1., 0.], [0., 0., 1.]],
    [[-1., 0., 0.], [-1., 1., 0.], [0., 0., -1.]],
    [[-1., 1., 0.], [-1., 0., 0.], [0., 0., 1.]],
    [[1., -1., 0.], [0., -1., 0.], [0., 0., -1.]],
    [[1., 0., 0.], [0., -1., 0.], [0., 0., -1.]],
    [[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]],
    [[0., 1., 0.], [1., 1., 0.], [0., 0., -1.]],
    [[-1., 0., 0.], [-1., -1., 0.], [0., 0., 1.]],
    [[-1., -1., 0.], [-1., 0., 0.], [0., 0., -1.]],
    [[1., 1., 0.], [0., 1., 0.], [0., 0., 1.]],
];

static SPACEGROUP_TO_HALL_NUMBER: [i32; 230] = [
    1, 2, 3, 6, 9, 18, 21, 30, 39, 57, 60, 63, 72, 81, 90, 108, 109, 112, 115, 116, 119, 122, 123,
    124, 125, 128, 134, 137, 143, 149, 155, 161, 164, 170, 173, 176, 182, 185, 191, 197, 203, 209,
    212, 215, 218, 221, 227, 228, 230, 233, 239, 245, 251, 257, 263, 266, 269, 275, 278, 284, 290,
    292, 298, 304, 310, 313, 316, 322, 334, 335, 337, 338, 341, 343, 349, 350, 351, 352, 353, 354,
    355, 356, 357, 358, 359, 361, 363, 364, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376,
    377, 378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395,
    396, 397, 398, 399, 400, 401, 402, 404, 406, 407, 408, 410, 412, 413, 414, 416, 418, 419, 420,
    422, 424, 425, 426, 428, 430, 431, 432, 433, 435, 436, 438, 439, 440, 441, 442, 443, 444, 446,
    447, 448, 449, 450, 452, 454, 455, 456, 457, 458, 460, 462, 463, 464, 465, 466, 467, 468, 469,
    470, 471, 472, 473, 474, 475, 476, 477, 478, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488,
    489, 490, 491, 492, 493, 494, 495, 497, 498, 500, 501, 502, 503, 504, 505, 506, 507, 508, 509,
    510, 511, 512, 513, 514, 515, 516, 517, 518, 520, 521, 523, 524, 525, 527, 529, 530,
];

static LAYER_GROUP_TO_HALL_NUMBER: [i32; 80] = [
    -1, -2, -3, -4, -5, -8, -9, -12, -14, -16, -18, -20, -22, -24, -26, -28, -30, -32, -34, -35,
    -37, -38, -39, -40, -42, -43, -44, -46, -48, -50, -52, -54, -56, -58, -60, -62, -64, -65, -67,
    -68, -70, -72, -74, -76, -77, -79, -80, -81, -82, -83, -84, -85, -87, -88, -89, -90, -91, -92,
    -93, -94, -95, -96, -98, -99, -101, -102, -103, -104, -105, -106, -107, -108, -109, -110, -111,
    -112, -113, -114, -115, -116,
];

// --- Structs ---

#[derive(Clone, Debug)]
pub struct Spacegroup {
    pub number: i32,
    pub hall_number: i32,
    pub pointgroup_number: i32,
    pub schoenflies: String,
    pub hall_symbol: String,
    pub international: String,
    pub international_long: String,
    pub international_short: String,
    pub choice: String,
    pub bravais_lattice: Mat3,
    pub origin_shift: Vec3,
}

impl Spacegroup {
    pub fn new() -> Self {
        Spacegroup {
            number: 0,
            hall_number: 0,
            pointgroup_number: 0,
            schoenflies: String::new(),
            hall_symbol: String::new(),
            international: String::new(),
            international_long: String::new(),
            international_short: String::new(),
            choice: String::new(),
            bravais_lattice: [[0.0; 3]; 3],
            origin_shift: [0.0; 3],
        }
    }
}

// --- Public Functions ---

pub fn spa_search_spacegroup(
    primitive: &Primitive,
    hall_number: i32,
    symprec: f64,
    angle_tolerance: f64,
) -> Option<Spacegroup> {
    debug::debug_print(format_args!(
        "search_spacegroup (tolerance = {}):\n",
        symprec
    ));

    let cell = primitive.cell.as_ref()?;
    let symmetry = sym_get_operation(cell, symprec, angle_tolerance)?;

    let candidates = if hall_number != 0 {
        vec![hall_number]
    } else if cell.aperiodic_axis == -1 {
        SPACEGROUP_TO_HALL_NUMBER.to_vec()
    } else {
        LAYER_GROUP_TO_HALL_NUMBER.to_vec()
    };

    search_spacegroup_with_symmetry(primitive, &candidates, &symmetry, symprec, angle_tolerance)
}

pub fn spa_search_spacegroup_with_symmetry(
    symmetry: &Symmetry,
    prim_lat: &Mat3,
    symprec: f64,
) -> Option<Spacegroup> {
    let mut primitive = Primitive::new(1);
    let mut cell = Cell::new(1, TensorRank::NoSpin);
    mat_copy_matrix_d3(&mut cell.lattice, prim_lat);
    cell.position[0] = [0.0; 3];
    primitive.cell = Some(cell);

    // Use all spacegroups (0-229)
    let candidates = SPACEGROUP_TO_HALL_NUMBER.to_vec();

    search_spacegroup_with_symmetry(&primitive, &candidates, symmetry, symprec, -1.0)
}

pub fn spa_transform_to_primitive(
    mapping_table: &mut [usize],
    cell: &Cell,
    trans_mat: &Mat3,
    centering: Centering,
    symprec: f64,
) -> Option<Cell> {
    let tmat_inv = match mat_inverse_matrix_d3(trans_mat, symprec) {
        Some(m) => m,
        None => return None,
    };

    let mut tmat = [[0.0; 3]; 3];
    match centering {
        Centering::Primitive => mat_copy_matrix_d3(&mut tmat, &tmat_inv),
        Centering::AFace => tmat = mat_multiply_matrix_d3(&tmat_inv, &A_MAT),
        Centering::CFace => tmat = mat_multiply_matrix_d3(&tmat_inv, &C_MAT),
        Centering::Face => tmat = mat_multiply_matrix_d3(&tmat_inv, &F_MAT),
        Centering::Body => tmat = mat_multiply_matrix_d3(&tmat_inv, &I_MAT),
        Centering::RCenter => tmat = mat_multiply_matrix_d3(&tmat_inv, &R_MAT),
        _ => return None,
    }

    let prim_lat = mat_multiply_matrix_d3(&cell.lattice, &tmat);

    let primitive = cel_trim_cell(mapping_table, &prim_lat, cell, symprec);
    if primitive.is_none() {
        debug::warning_print(format_args!("spglib: cel_trim_cell failed.\n"));
    }
    primitive
}

pub fn spa_transform_from_primitive(
    primitive: &Cell,
    centering: Centering,
    symprec: f64,
) -> Option<Cell> {
    let mut tmat = [[0.0; 3]; 3];
    let mut inv_tmat = [[0.0; 3]; 3];

    match centering {
        Centering::Primitive => {}
        Centering::AFace => {
            mat_copy_matrix_d3(&mut tmat, &A_MAT);
            inv_tmat = mat_inverse_matrix_d3(&A_MAT, 0.0).unwrap();
        }
        Centering::CFace => {
            mat_copy_matrix_d3(&mut tmat, &C_MAT);
            inv_tmat = mat_inverse_matrix_d3(&C_MAT, 0.0).unwrap();
        }
        Centering::Face => {
            mat_copy_matrix_d3(&mut tmat, &F_MAT);
            inv_tmat = mat_inverse_matrix_d3(&F_MAT, 0.0).unwrap();
        }
        Centering::Body => {
            mat_copy_matrix_d3(&mut tmat, &I_MAT);
            inv_tmat = mat_inverse_matrix_d3(&I_MAT, 0.0).unwrap();
        }
        Centering::RCenter => {
            mat_copy_matrix_d3(&mut tmat, &R_MAT);
            inv_tmat = mat_inverse_matrix_d3(&R_MAT, 0.0).unwrap();
        }
        _ => return None,
    }

    let mut shift = [[0.0; 3]; 3];
    let multi = get_centering_shifts(&mut shift, centering);

    let mut mapping_table = vec![0; primitive.size * multi];
    let mut std_cell = Cell::new(primitive.size * multi, primitive.tensor_rank);

    std_cell.lattice = mat_multiply_matrix_d3(&primitive.lattice, &inv_tmat);

    let mut num_atom = 0;
    for i in 0..primitive.size {
        std_cell.position[num_atom] = mat_multiply_matrix_vector_d3(&tmat, &primitive.position[i]);
        std_cell.types[num_atom] = primitive.types[i];
        num_atom += 1;
    }

    for i in 0..(multi - 1) {
        for j in 0..primitive.size {
            let src_pos = std_cell.position[j];
            std_cell.position[num_atom] = src_pos;
            for k in 0..3 {
                std_cell.position[num_atom][k] += shift[i][k];
            }
            std_cell.types[num_atom] = std_cell.types[j];
            num_atom += 1;
        }
    }

    let trimmed_cell = cel_trim_cell(&mut mapping_table, &std_cell.lattice, &std_cell, symprec);
    trimmed_cell
}

pub fn spa_copy_spacegroup(dst: &mut Spacegroup, src: &Spacegroup) {
    dst.number = src.number;
    dst.hall_number = src.hall_number;
    dst.pointgroup_number = src.pointgroup_number;
    dst.schoenflies = src.schoenflies.clone();
    dst.hall_symbol = src.hall_symbol.clone();
    dst.international = src.international.clone();
    dst.international_long = src.international_long.clone();
    dst.international_short = src.international_short.clone();
    dst.choice = src.choice.clone();
    mat_copy_matrix_d3(&mut dst.bravais_lattice, &src.bravais_lattice);
    mat_copy_vector_d3(&mut dst.origin_shift, &src.origin_shift);
}

// --- Internal Functions ---

fn search_spacegroup_with_symmetry(
    primitive: &Primitive,
    candidates: &[i32],
    symmetry: &Symmetry,
    symprec: f64,
    angle_tolerance: f64,
) -> Option<Spacegroup> {
    debug::debug_print(format_args!(
        "search_spacegroup (tolerance = {}):\n",
        symprec
    ));

    let mut origin_shift = [0.0; 3];
    let mut conv_lattice = [[0.0; 3]; 3];

    let pointsym = ptg_get_pointsymmetry(&symmetry.rot);
    if pointsym.size < symmetry.size {
        debug::info_print(format_args!(
            "spglib: Point symmetry of primitive cell is broken.\n"
        ));
        return None;
    }

    let hall_number = iterative_search_hall_number(
        &mut origin_shift,
        &mut conv_lattice,
        candidates,
        primitive,
        symmetry,
        symprec,
        angle_tolerance,
    );

    if hall_number == 0 {
        return None;
    }

    get_spacegroup(hall_number, &origin_shift, &conv_lattice)
}

fn get_spacegroup(
    hall_number: i32,
    origin_shift: &Vec3,
    conv_lattice: &Mat3,
) -> Option<Spacegroup> {
    let spg_type = spgdb_get_spacegroup_type(hall_number)?;

    let mut spacegroup = Spacegroup::new();
    mat_copy_matrix_d3(&mut spacegroup.bravais_lattice, conv_lattice);
    mat_copy_vector_d3(&mut spacegroup.origin_shift, origin_shift);
    spacegroup.number = spg_type.number;
    spacegroup.hall_number = hall_number;
    spacegroup.pointgroup_number = spg_type.pointgroup_number;
    spacegroup.schoenflies = spg_type.schoenflies;
    spacegroup.hall_symbol = spg_type.hall_symbol;
    spacegroup.international = spg_type.international;
    spacegroup.international_long = spg_type.international_full;
    spacegroup.international_short = spg_type.international_short;
    spacegroup.choice = spg_type.choice;

    Some(spacegroup)
}

fn iterative_search_hall_number(
    origin_shift: &mut Vec3,
    conv_lattice: &mut Mat3,
    candidates: &[i32],
    primitive: &Primitive,
    symmetry: &Symmetry,
    symprec: f64,
    angle_tolerance: f64,
) -> i32 {
    debug::debug_print(format_args!("iterative_search_hall_number:\n"));

    let mut hall_number = search_hall_number(
        origin_shift,
        conv_lattice,
        candidates,
        primitive,
        symmetry,
        symprec,
    );

    if hall_number != 0 {
        return hall_number;
    }

    let mut tolerance = symprec;
    let mut current_symmetry = symmetry.clone();

    for attempt in 0..NUM_ATTEMPT {
        debug::debug_print(format_args!(
            "spglib: Attempt {} tolerance = {:e} failed",
            attempt, tolerance
        ));

        tolerance *= REDUCE_RATE;
        if let Some(sym_reduced) = sym_reduce_operation(
            primitive.cell.as_ref().unwrap(),
            &current_symmetry,
            tolerance,
            angle_tolerance,
        ) {
            hall_number = search_hall_number(
                origin_shift,
                conv_lattice,
                candidates,
                primitive,
                &sym_reduced,
                symprec,
            );
            current_symmetry = sym_reduced;
        }

        if hall_number != 0 {
            break;
        }
    }

    hall_number
}

fn search_hall_number(
    origin_shift: &mut Vec3,
    conv_lattice: &mut Mat3,
    candidates: &[i32],
    primitive: &Primitive,
    symmetry: &Symmetry,
    symprec: f64,
) -> i32 {
    debug::debug_print(format_args!("search_hall_number:\n"));

    let aperiodic_axis = primitive.cell.as_ref().unwrap().aperiodic_axis;
    let mut tmat_int = [[0; 3]; 3];

    let pointgroup = ptg_get_transformation_matrix(&mut tmat_int, &symmetry.rot, aperiodic_axis);

    // debug_print("initial transformation matrix\n");
    // debug_print_matrix_i3(tmat_int);

    if pointgroup.number == 0 {
        return 0;
    }

    let mut conv_lattice_tmp = [[0.0; 3]; 3];

    if pointgroup.laue == Laue::Laue1 || pointgroup.laue == Laue::Laue2M {
        conv_lattice_tmp =
            mat_multiply_matrix_di3(&primitive.cell.as_ref().unwrap().lattice, &tmat_int);

        if pointgroup.laue == Laue::Laue1 {
            if !change_basis_tricli(
                &mut tmat_int,
                &conv_lattice_tmp,
                &primitive.cell.as_ref().unwrap().lattice,
                symprec,
                aperiodic_axis,
            ) {
                return 0;
            }
        }

        if pointgroup.laue == Laue::Laue2M {
            if !change_basis_monocli(
                &mut tmat_int,
                &conv_lattice_tmp,
                &primitive.cell.as_ref().unwrap().lattice,
                symprec,
                aperiodic_axis,
            ) {
                return 0;
            }
        }
    }

    let mut correction_mat = [[0.0; 3]; 3];
    let centering = get_centering(&mut correction_mat, &tmat_int, pointgroup.laue);
    if centering == Centering::Error {
        return 0;
    }

    let tmat = mat_multiply_matrix_id3(&tmat_int, &correction_mat);
    *conv_lattice = mat_multiply_matrix_d3(&primitive.cell.as_ref().unwrap().lattice, &tmat);

    // debug_print("transformation matrix\n");
    // debug_print_matrix_d3(tmat);

    let conv_symmetry = get_initial_conventional_symmetry(centering, &tmat, symmetry);
    if conv_symmetry.is_none() {
        return 0;
    }
    let conv_symmetry = conv_symmetry.unwrap();

    for &cand in candidates {
        if hal_match_hall_symbol_db(
            origin_shift,
            conv_lattice,
            cand,
            centering,
            &conv_symmetry,
            symprec,
        ) {
            debug::debug_print(format_args!("origin shift\n"));
            // debug_print_vector_d3(origin_shift);
            return cand;
        }
    }

    0
}

fn change_basis_tricli(
    tmat_int: &mut Mat3I,
    conv_lattice: &Mat3,
    primitive_lattice: &Mat3,
    symprec: f64,
    aperiodic_axis: i32,
) -> bool {
    let mut niggli_cell = *conv_lattice;

    if !niggli_reduce(&mut niggli_cell, symprec * symprec, aperiodic_axis) {
        return false;
    }

    let mut smallest_lattice = niggli_cell;
    if mat_get_determinant_d3(&smallest_lattice) < 0.0 {
        for i in 0..3 {
            for j in 0..3 {
                smallest_lattice[i][j] = -smallest_lattice[i][j];
            }
        }
    }

    let inv_lattice = mat_inverse_matrix_d3(primitive_lattice, 0.0).unwrap();
    let tmat = mat_multiply_matrix_d3(&inv_lattice, &smallest_lattice);
    *tmat_int = mat_cast_matrix_3d_to_3i(&tmat);

    true
}

fn change_basis_monocli(
    tmat_int: &mut Mat3I,
    conv_lattice: &Mat3,
    primitive_lattice: &Mat3,
    symprec: f64,
    aperiodic_axis_prim: i32,
) -> bool {
    let mut smallest_lattice = [[0.0; 3]; 3];
    let mut aperiodic_axis_conv = -1;
    let unique_axis;

    if aperiodic_axis_prim == -1 {
        unique_axis = 1;
    } else {
        for i in 0..3 {
            if tmat_int[aperiodic_axis_prim as usize][i] != 0 {
                aperiodic_axis_conv = i as i32;
            }
        }
        unique_axis = 0;
    }

    if !del_layer_delaunay_reduce_2D(
        &mut smallest_lattice,
        conv_lattice,
        unique_axis,
        aperiodic_axis_conv,
        symprec,
    ) {
        return false;
    }

    if aperiodic_axis_conv == 0 {
        // CHANGE_OF_BASIS_MONOCLI[2] corresponds to index 2 in C array
        smallest_lattice = mat_multiply_matrix_d3(&smallest_lattice, &CHANGE_OF_BASIS_MONOCLI[2]);
    }

    let inv_lattice = mat_inverse_matrix_d3(primitive_lattice, 0.0).unwrap();
    let tmat = mat_multiply_matrix_d3(&inv_lattice, &smallest_lattice);
    *tmat_int = mat_cast_matrix_3d_to_3i(&tmat);

    true
}

fn get_initial_conventional_symmetry(
    centering: Centering,
    tmat: &Mat3,
    symmetry: &Symmetry,
) -> Option<Symmetry> {
    debug::debug_print(format_args!("get_initial_conventional_symmetry\n"));

    if centering == Centering::RCenter {
        get_conventional_symmetry(tmat, Centering::Primitive, symmetry)
    } else {
        get_conventional_symmetry(tmat, centering, symmetry)
    }
}

fn get_conventional_symmetry(
    tmat: &Mat3,
    centering: Centering,
    primitive_sym: &Symmetry,
) -> Option<Symmetry> {
    let size = primitive_sym.size;
    let mut symmetry;

    match centering {
        Centering::Face => symmetry = Symmetry::new(size * 4),
        Centering::RCenter => symmetry = Symmetry::new(size * 3),
        Centering::Body | Centering::AFace | Centering::BFace | Centering::CFace => {
            symmetry = Symmetry::new(size * 2)
        }
        _ => symmetry = Symmetry::new(size),
    }

    let inv_tmat = mat_inverse_matrix_d3(tmat, 0.0).unwrap_or([[0.0; 3]; 3]);

    for i in 0..size {
        let mut primitive_sym_rot_d3 = [[0.0; 3]; 3];
        mat_cast_matrix_3i_to_3d(&primitive_sym.rot[i]);

        // C*S*C^-1
        let mut symmetry_rot_d3 = [[0.0; 3]; 3];
        if let Some(res) = mat_get_similar_matrix_d3(&primitive_sym_rot_d3, tmat, 0.0) {
            symmetry_rot_d3 = res;
        }
        symmetry.rot[i] = mat_cast_matrix_3d_to_3i(&symmetry_rot_d3);

        // translation: C = B^-1*P
        symmetry.trans[i] = mat_multiply_matrix_vector_d3(&inv_tmat, &primitive_sym.trans[i]);
    }

    let mut shift = [[0.0; 3]; 3];
    let multi = get_centering_shifts(&mut shift, centering);

    if centering != Centering::Primitive {
        for i in 0..(multi - 1) {
            for j in 0..size {
                //mat_copy_matrix_i3(&mut symmetry.rot[(i + 1) * size + j], &symmetry.rot[j]);
                let src_rot = symmetry.rot[j];
                symmetry.rot[(i + 1) * size + j] = src_rot;
                for k in 0..3 {
                    symmetry.trans[(i + 1) * size + j][k] = symmetry.trans[j][k] + shift[i][k];
                }
            }
        }
    }

    Some(symmetry)
}

fn get_centering(correction_mat: &mut Mat3, tmat: &Mat3I, laue: Laue) -> Centering {
    mat_copy_matrix_d3(correction_mat, &IDENTITY);
    let det = mat_get_determinant_i3(tmat).abs();
    debug::debug_print(format_args!("laue class: {:?}\n", laue));
    debug::debug_print(format_args!("multiplicity: {}\n", det));

    match det {
        1 => Centering::Primitive,
        2 => {
            let mut centering = get_base_center(tmat);
            if centering == Centering::AFace {
                if laue == Laue::Laue2M {
                    debug::debug_print(format_args!("Monocli A to C\n"));
                    mat_copy_matrix_d3(correction_mat, &MONOCLI_A2C);
                } else {
                    mat_copy_matrix_d3(correction_mat, &A2C);
                }
                centering = Centering::CFace;
            }
            if centering == Centering::BFace {
                mat_copy_matrix_d3(correction_mat, &B2C);
                centering = Centering::CFace;
            }
            if laue == Laue::Laue2M && centering == Centering::Body {
                debug::debug_print(format_args!("Monocli I to C\n"));
                mat_copy_matrix_d3(correction_mat, &MONOCLI_I2C);
                centering = Centering::CFace;
            }
            centering
        }
        3 => {
            let mut trans_corr_mat = [[0.0; 3]; 3];
            trans_corr_mat = mat_multiply_matrix_id3(tmat, &RHOMBO_OBVERSE);
            if mat_is_int_matrix(&trans_corr_mat, INT_PREC) {
                mat_copy_matrix_d3(correction_mat, &RHOMBO_OBVERSE);
                debug::debug_print(format_args!("R-center observe setting\n"));
                return Centering::RCenter;
            }
            trans_corr_mat = mat_multiply_matrix_id3(tmat, &RHOMBO_REVERSE);
            if mat_is_int_matrix(&trans_corr_mat, INT_PREC) {
                mat_copy_matrix_d3(correction_mat, &RHOMBO_REVERSE);
                debug::debug_print(format_args!("R-center reverse setting\n"));
                return Centering::RCenter;
            }
            Centering::RCenter
        }
        4 => Centering::Face,
        _ => Centering::Error,
    }
}

fn get_base_center(tmat: &Mat3I) -> Centering {
    debug::debug_print(format_args!("lat_get_base_center\n"));

    // C center
    for i in 0..3 {
        if tmat[i][0] == 0 && tmat[i][1] == 0 && tmat[i][2].abs() == 1 {
            return Centering::CFace;
        }
    }

    // A center
    for i in 0..3 {
        if tmat[i][0].abs() == 1 && tmat[i][1] == 0 && tmat[i][2] == 0 {
            return Centering::AFace;
        }
    }

    // B center
    for i in 0..3 {
        if tmat[i][0] == 0 && tmat[i][1].abs() == 1 && tmat[i][2] == 0 {
            return Centering::BFace;
        }
    }

    // Body center
    if (tmat[0][0].abs() + tmat[0][1].abs() + tmat[0][2].abs() == 2)
        && (tmat[1][0].abs() + tmat[1][1].abs() + tmat[1][2].abs() == 2)
        && (tmat[2][0].abs() + tmat[2][1].abs() + tmat[2][2].abs() == 2)
    {
        return Centering::Body;
    }

    debug::warning_print(format_args!("spglib: No centring was found.\n"));
    Centering::Primitive
}

fn get_centering_shifts(shift: &mut [[f64; 3]; 3], centering: Centering) -> usize {
    let mut multi = 1;
    for i in 0..3 {
        shift[i] = [0.0; 3];
    }

    if centering != Centering::Primitive {
        if centering != Centering::Face && centering != Centering::RCenter {
            for i in 0..3 {
                shift[0][i] = 0.5;
            } // BODY
            if centering == Centering::AFace {
                shift[0][0] = 0.0;
            }
            if centering == Centering::BFace {
                shift[0][1] = 0.0;
            }
            if centering == Centering::CFace {
                shift[0][2] = 0.0;
            }
            multi = 2;
        }

        if centering == Centering::RCenter {
            shift[0][0] = 2.0 / 3.0;
            shift[0][1] = 1.0 / 3.0;
            shift[0][2] = 1.0 / 3.0;
            shift[1][0] = 1.0 / 3.0;
            shift[1][1] = 2.0 / 3.0;
            shift[1][2] = 2.0 / 3.0;
            multi = 3;
        }

        if centering == Centering::Face {
            shift[0][0] = 0.0;
            shift[0][1] = 0.5;
            shift[0][2] = 0.5;
            shift[1][0] = 0.5;
            shift[1][1] = 0.0;
            shift[1][2] = 0.5;
            shift[2][0] = 0.5;
            shift[2][1] = 0.5;
            shift[2][2] = 0.0;
            multi = 4;
        }
    }
    multi
}

fn is_equivalent_lattice(
    tmat: &mut Mat3,
    mode: i32,
    lattice: &Mat3,
    orig_lattice: &Mat3,
    symprec: f64,
) -> bool {
    if (mat_get_determinant_d3(lattice) - mat_get_determinant_d3(orig_lattice)).abs() > symprec {
        return false;
    }

    let mut inv_lat = [[0.0; 3]; 3];
    if mat_inverse_matrix_d3(lattice, symprec).is_none() {
        return false;
    }

    // orig_lattice == lattice @ tmat
    *tmat = mat_multiply_matrix_d3(&inv_lat, orig_lattice);

    match mode {
        0 => {
            if mat_check_identity_matrix_d3(&IDENTITY, tmat, symprec) {
                return true;
            }
        }
        1 => {
            let mut tmat_abs = [[0.0; 3]; 3];
            for i in 0..3 {
                for j in 0..3 {
                    tmat_abs[i][j] = mat_dabs(tmat[i][j]);
                }
            }
            if mat_check_identity_matrix_d3(&IDENTITY, &tmat_abs, symprec) {
                return true;
            }
        }
        2 => {
            let tmat_int = mat_cast_matrix_3d_to_3i(tmat);
            if mat_check_identity_matrix_id3(&tmat_int, tmat, symprec) {
                if mat_get_determinant_i3(&tmat_int) == 1 {
                    let orig_metric = mat_get_metric(orig_lattice);
                    let metric = mat_get_metric(lattice);
                    if mat_check_identity_matrix_d3(&orig_metric, &metric, symprec) {
                        return true;
                    }
                }
            }
        }
        _ => {}
    }

    false
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cell::Cell;
    use crate::primitive::Primitive;
    use crate::symmetry::Symmetry;

    #[test]
    fn test_search_spacegroup_simple() {
        // Construct a simple P1 cell
        let lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let position = [[0.0, 0.0, 0.0]];
        let types = [1];
        let mut cell = Cell::new(1, crate::cell::TensorRank::NoSpin);
        cell.set_cell(&lattice, &position, &types);

        let mut primitive = Primitive::new(1);
        primitive.cell = Some(cell);
        mat_copy_matrix_d3(&mut primitive.orig_lattice, &lattice);

        // Should find P1 (Hall number 1)
        // Note: A single atom at origin in cubic lattice is actually Pm-3m (221).
        // If we want P1, we should break symmetry, e.g. triclinic lattice.
        // But let's see what it finds. It should find something valid.
        let spg = spa_search_spacegroup(&primitive, 0, 1e-5, -1.0);
        assert!(spg.is_some());
        let s = spg.unwrap();
        // For cubic lattice with 1 atom, it is Pm-3m (221)
        assert_eq!(s.number, 221);
    }
}
