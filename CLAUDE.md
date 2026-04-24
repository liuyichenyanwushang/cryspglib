# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build & Test Commands

```bash
cargo build              # Build the library
cargo test               # Run all 69 tests (magnetic: ~24 end-to-end)
cargo test test_fcc -- --nocapture  # FCC FM 测试（含已知 bug 的原始输出）
cargo test -- --nocapture  # Run tests with debug output (eprintln visible)
cargo test <test_name>   # Run a single test (e.g., `cargo test test_graphene_p6mmm`)
cargo check              # Check compilation without building
cargo clippy             # Run lints
```

The project uses `edition = "2024"` (requires nightly Rust).

## Critical Convention: Matrix Layout `lattice[cart][vec]`

All 3x3 matrices (`Mat3 = [[f64; 3]; 3]`) use the convention **rows = Cartesian components, columns = lattice vectors**:

```
lattice = [[a_x, b_x, c_x],
           [a_y, b_y, c_y],
           [a_z, b_z, c_z]]
```

- `lattice[i][j]` — the i-th Cartesian component (0=x, 1=y, 2=z) of the j-th lattice vector (0=a, 1=b, 2=c)
- Each **column** is a lattice vector; each **row** is a Cartesian direction

**Why this matters:**

- `mat_get_metric(lattice) = L^T * L` computes the metric tensor G where G[i][j] = dot(vec_i, vec_j). This only works when columns are vectors.
- `mat_multiply_matrix_vector_d3(lattice, frac)` transforms fractional coordinates to Cartesian: `cart = L * frac`
- Cubic matrices are symmetric (a=b=c, all angles 90°), so both conventions look identical — this is why cubic tests passed before the convention was corrected
- Hexagonal matrices are asymmetric (`a != b` direction), so swapping conventions silently produces wrong results (e.g., graphene identified as space group #10 instead of #191)

**When constructing test lattices**, always use `lattice[cart][vec]`:
```rust
// Hexagonal: a along x, b at 120° in xy-plane, c along z
fn hexagonal_lattice(a: f64, c: f64) -> Mat3 {
    let sqrt3 = 3.0_f64.sqrt();
    [
        [a, -a / 2.0, 0.0],      // x components of a, b, c
        [0.0, a * sqrt3 / 2.0, 0.0],  // y components
        [0.0, 0.0, c],            // z components
    ]
}
```

**When reading C reference code**, note that C uses the same memory layout (`double lattice[3][3]` with first index as row/cartesian), but the C code's array initialization may be transposed relative to Rust's `[[f64; 3]; 3]` depending on how it's written.

## Architecture: Space Group Identification Pipeline

This is a pure Rust port of the C library `spglib` with zero external dependencies.

The pipeline for identifying a crystal's space group:

1. **Cell** (`cell.rs`) — crystal structure input: lattice vectors, atomic positions, types, overlap checks
2. **Symmetry detection** (`symmetry.rs`) — find symmetry operations (rotations + translations) of the cell
3. **Point group classification** (`pointgroup.rs`) — classify the rotation part into a crystallographic point group, determine Laue class and holohedry
4. **Lattice reduction** (`delaunay.rs`, `niggli.rs`) — reduce lattice to canonical form (Delaunay for general, Niggli for triclinic/monoclinic)
5. **Primitive cell** (`primitive.rs`) — find primitive cell and its symmetry operations
6. **Hall symbol matching** (`hall_symbol.rs`) — match primitive cell symmetry against the Hall symbol database (`spg_database.rs`) to identify the Hall number (1–530)
7. **Space group search** (`spacegroup.rs`) — iterate Hall number candidates, perform origin shift and conventional cell transformation, select best match; returns the space group number (1–230)
8. **Structure refinement** (`refinement.rs`) — rotate Bravais lattice to standard orientation, determine exact Wyckoff positions via site symmetry database, recover symmetry operations in the original cell. Produces `ExactStructure` with standardized positions, Wyckoff labels, equivalent atom mappings.
9. **Determination** (`determination.rs`) — top-level retry loop that combines primitive + spacegroup + refinement steps with decreasing tolerance
10. **Public API** (`lib.rs`) — wraps the pipeline behind a clean Rust API: `SpglibDataset`, `spg_get_dataset`, `spg_standardize_cell`, etc.

Key supporting modules:
- `mathfunc.rs` — 3x3 matrix/vector operations with type aliases `Mat3`, `Mat3I`, `Vec3`, `Vec3I` (no external BLAS/LAPACK)
- `arithmetic.rs` — integer arithmetic for symmetry checks
- `overlap.rs` — overlap detection for atomic positions
- `site_symmetry.rs` — Wyckoff position exact location and equivalent atom assignment using site symmetry database
- `sitesym_database.rs` — site symmetry database (Wyckoff positions, site symmetry symbols)
- `spin.rs` — magnetic symmetry types
- `kgrid.rs`, `kpoint.rs` — k-point grid generation (Monkhorst-Pack)
- `debug.rs` — debug/warning print macros

## Critical Porting Pattern: C Output Parameters → Rust Return Values

In the C source (`src_c/*.txt`), many functions write results through pointer parameters:
```c
void mat_cast_matrix_3i_to_3d(double out[3][3], int const a[3][3]);
```

The Rust port converts these to return values:
```rust
pub fn mat_cast_matrix_3i_to_3d(a: &Mat3I) -> Mat3;
```

**When comparing Rust against C code, always check that the return value is actually captured.** The most severe bug found in this port was ignoring the return value of `mat_cast_matrix_3i_to_3d` in `get_conventional_symmetry`, which caused all rotation matrices to be zero and silently broke space group identification.

## Rust API Design Conventions (C → Rust mapping)

The `lib.rs` public API follows these conventions when porting from C:

| C Pattern | Rust Pattern |
|-----------|-------------|
| `malloc`/`free` | `Vec` (ownership, auto-drop) |
| `NULL` return on failure | `Option<T>` (None = failure) |
| `thread_local` error code + return 0 | `Option<T>` return |
| Output parameter `double out[3][3]` | Return `Mat3` |
| Output parameter `int *types` + `int *num_atom` | Return `Cell` (owns lattice, positions, types) |
| `spg_free_dataset(dataset)` | Not needed — Rust drops `SpglibDataset` automatically |
| `char[7]` / `char[11]` buffers | `String` |
| `int *mapping_table` pre-allocated by caller | `Vec<i32>` returned |
| K-point pre-allocated grid arrays | `&mut [i32]` / `&mut [usize]` slices |

## Reference C Code

`src_c/` contains the original spglib C source (`*.txt` files, renamed from `*.c` for tooling compatibility). Most `.rs` files correspond to a `.txt` file of the same name. The magnetic modules (`msg_database`, `magnetic_spacegroup`) are partially ported.

| C file | Rust module | Purpose |
|--------|------------|---------|
| `spglib.txt` | `lib.rs` | Public API, dataset, standardization |
| `determination.txt` | `determination.rs` | Top-level retry loop |
| `refinement.txt` | `refinement.rs` | Bravais lattice refinement, Wyckoff positions |
| `site_symmetry.txt` | `site_symmetry.rs` | Site symmetry exact positions |
| `spacegroup.txt` | `spacegroup.rs` | Space group search and identification |
| `primitive.txt` | `primitive.rs` | Primitive cell finding |
| `hall_symbol.txt` | `hall_symbol.rs` | Hall symbol matching |
| `pointgroup.txt` | `pointgroup.rs` | Point group classification |
| `spg_database.txt` | `spg_database.rs` | Space group database |
| `sitesym_database.txt` | `sitesym_database.rs` | Site symmetry database |
| `mathfunc.txt` | `mathfunc.rs` | Matrix/vector math |
| `delaunay.txt` | `delaunay.rs` | Delaunay lattice reduction |
| `niggli.txt` | `niggli.rs` | Niggli lattice reduction |
| `kgrid.txt` | `kgrid.rs` | k-point grid generation |
| `kpoint.txt` | `kpoint.rs` | k-point address mapping |
| `symmetry.txt` | `symmetry.rs` | Symmetry operation search |
| `cell.txt` | `cell.rs` | Crystal cell and overlap checks |
| `msg_database.txt` | `msg_database.rs` + `msg_database_gen.rs` | Magnetic space group database (auto-converted from C) |
| `magnetic_spacegroup.txt` | `magnetic_spacegroup.rs` | Magnetic space group identification |
| `spin.txt` | `spin.rs` | Magnetic symmetry operations with site tensors |

When debugging, compare Rust implementation against the C original line-by-line, paying special attention to:
- Output parameter → return value conversions
- Array indexing differences (C is row-major, Rust uses same layout)
- Floating point tolerance constants

## Test Coverage

Tests in `spacegroup.rs` verify the full pipeline end-to-end:
- Cubic: simple (#221), bcc (#229), fcc (#225), diamond (#227), rocksalt (#225)
- Hexagonal: graphene (#191), hcp (#194), AB-buckled silicene (#164)
- Supercells: 2x2x2 simple cubic (#221), 2x2x2 CsCl (#221), 2x2x1 graphene (#191)

Magnetic space group tests in `magnetic_spacegroup_test.rs`:
- Type-1/2/3 identification via Pm-3m (#221) symmetry operations from database
- Public API `spg_get_magnetic_spacegroup_type_from_symmetry`
- Edge cases: empty symmetry, missing identity

## Magnetic Space Group Status

Magnetic space group identification is fully functional for direct magnetic symmetry analysis.

### Implemented Components

| Component | Lines | Description |
|-----------|-------|-------------|
| `convert_msgdb.py` | ~200 | Python script to auto-convert C data arrays to Rust static data |
| `msg_database_gen.rs` | **28,319** | Generated static arrays (6 arrays: 1652 MSG types, 531 Hall mappings, 76683 encoded operations, etc.) |
| `msg_database.rs` | ~120 | Database query functions |
| `magnetic_spacegroup.rs` | ~500 | Magnetic identification pipeline: FSG/XSG extraction, type classification (1-4), database matching, auto-detect lower-symmetry parent |
| `magnetic_spacegroup_test.rs` | ~600 | 24 tests: Type-1/2/3, Fe/Cu AFM, FCC FM, end-to-end, edge cases |
| `spin.rs` | ~650 | Magnetic symmetry from site tensors (`spn_get_operations_with_site_tensors`) |
| `lib.rs` public API | ~400 | `SpglibMagneticSymmetry`, `spg_get_magnetic_dataset`, `spg_read_structure`, `spg_format_magnetic_symmetry` |

### Public API: End-to-End Pipeline

```rust
/// 从晶格 + 原子位置 + 磁矩分析磁空间群和对称操作。
pub fn spg_get_magnetic_dataset(
    lattice: &Mat3,
    positions: &[Vec3],
    types: &[i32],
    magnetic_moments: Option<&[[f64; 3]]>,  // None = 非磁
    symprec: f64,
) -> Option<SpglibMagneticSymmetry>;

/// 从类似 POSCAR 的格式解析结构。
pub fn spg_read_structure(data: &str) -> Option<(Mat3, Vec<Vec3>, Vec<i32>, Option<Vec<[f64; 3]>>)>;

/// 格式化输出（类似 phonopy --symmetry 风格）。
pub fn spg_format_magnetic_symmetry(result: &SpglibMagneticSymmetry) -> String;
```

管线流程:
1. `prm_get_primitive` + `spa_search_spacegroup` → 非磁空间群
2. `sym_get_operation` → 非磁对称操作
3. `spn_get_operations_with_site_tensors` → 磁对称操作（含 timerev 标记）
4. `msg_identify_magnetic_space_group_type` → 磁空间群 UNI + BNS 符号

### Database Architecture

```
HALL_MAPPING[hall_number] → [min_uni, max_uni]        # Hall → UNI range
UNI_MAPPING[uni_number]   → [num_halls, first_hall]    # UNI → Hall mapping
OPERATION_INDEX[uni][offset] → [count, start]           # UNI → operation index
```

Each operation in `MAGNETIC_SYMMETRY_OPERATIONS` is a single i32 encoding:
- Rotation matrix (base-3) + translation vector (base-12) + time-reversal flag

### Verified Test Cases (69 tests)

| Test | Structure | Non-mag Spg | Magnetic Result |
|------|-----------|------------|----------------|
| Type-1/2/3 (database ops) | Pm-3m hall 517 | #221 | UNI 1594-1598 |
| Fe [1,0,0] body center | 1 atom @ SC body center | Pm-3m (#221) | P4/mmm (#123), Type-3 |
| Fe [1,1,1] body center | 1 atom @ SC body center | Pm-3m (#221) | R-3m (#166), UNI=1331 |
| Co AFM [111] BCC | 2 atoms BCC | Im-3m (#229) | R-3m (#166), UNI=1331 |
| 8-atom AFM [111] | 8 atoms SC | Pm-3m (#221) | 24 ops (collinear AFM) |
| 8-atom AFM [001] | 8 atoms SC | Pm-3m (#221) | 32 ops (z-axis collinear AFM) |
| FCC FM [001] | 4 atoms FCC | Fm-3m (#225) | Type-4, 4 ops, UNI=0 |
| FCC FM [111] | 4 atoms FCC | Fm-3m (#225) | Type-4, 4 ops, UNI=0 |
| End-to-end POSCAR | Fe SC body center | Pm-3m (#221) | 16 ops, Type-3 |
| End-to-end POSCAR AFM | Fe BCC AFM | Im-3m (#229) | R-3m, UNI=1331 |

### Known Issue: `spn_get_operations_with_site_tensors` Base Mismatch

**问题描述**: `spn_get_operations_with_site_tensors` 接收的「非磁对称操作」来自**原胞**
（`sym_get_operation` on primitive cell），但这些操作是在原胞的晶格基（如 FCC 的菱面体基）
下的整数旋转矩阵。而函数内部的原子位置匹配 (`apply_symmetry_to_position`) 直接用这些
旋转矩阵去映射**常规晶胞**（如 FCC 的立方基）的原子坐标分数——基矢不同导致位置映射失败。

**影响**: 对面心/体心等非常规原胞（BCC、FCC 等），`spn_get_operations_with_site_tensors`
只能找到少量磁对称操作（如 FCC 铁磁只有 4 个），因为大多数操作的位置映射不正确。

**具体表现 (FCC FM [001])**:
- 非磁: Fm-3m (#225)，常规晶胞 4 原子，原胞 1 原子
- `sym_get_operation` → 48 个原胞对称操作（菱面体基）
- `spn_get_operations_with_site_tensors`: 用这 48 个操作去匹配 4 个立方基原子的磁矩
- 位置匹配无法正确对应 → 只保留 4 个有效操作 → Type-4
- 实际上 FCC 铁磁应有更多对称操作

**根因**: `spin.rs` 中 `get_operations` 函数的原子位置循环:
```rust
// rot 来自原胞（如菱面体基），position 来自常规晶胞（立方基）
apply_symmetry_to_position(&mut pos, &cell.position[j],
    &sym_nonspin.rot[i], &sym_nonspin.trans[i]);  // 基矢不匹配
```

**修复方向**: 将原胞旋转矩阵转换到常规晶胞基后，用转换后的矩阵进行位置匹配：
```rust
// 先用 lattice 将旋转从原胞基转到笛卡尔，再转到常规晶胞基
let rot_conventional = convert_to_conventional_basis(sym_nonspin.rot[i], lattice, prim_lattice);
apply_symmetry_to_position(&mut pos, &cell.position[j], &rot_conventional, &trans);
```

**当前状态**: ❌ 未修复 → 影响 BCC/FCC 等非常规原胞的磁对称操作数

## Notes

- `debug_test.rs` and similar ad-hoc files in the repo root are gitignored/untracked — use them for quick experiments but don't commit.
- `periodic_axes` in `cell.rs` layer functions uses `[usize; 2]` through the `AperiodicAxis` enum (not raw `i32`).
