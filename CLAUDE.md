# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build & Test Commands

```bash
cargo build              # Build the library
cargo test               # Run all 69 tests (magnetic: 24+ end-to-end)
cargo test test_fcc -- --nocapture  # FCC FM 测试
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

## Magnetic Space Group: Complete Pipeline

End-to-end: **POSCAR-like input → non-magnetic space group + magnetic space group (BNS) + symmetry operations**.

### API

```rust
pub fn spg_get_magnetic_dataset(
    lattice, positions, types, magnetic_moments, symprec
) -> Option<SpglibMagneticSymmetry>;

pub fn spg_read_structure(data) -> Option<(Mat3, Vec<Vec3>, Vec<i32>, Option<Vec<[f64;3]>>)>;

pub fn spg_format_magnetic_symmetry(result) -> String;
```

### Pipeline Flow

```
POSCAR → Cell (TensorRank::NonCollinear)
  → prm_get_primitive + spa_search_spacegroup → 非磁空间群
  → sym_get_operation(&cell)                  → 非磁对称操作（常规晶胞基）
  → spn_get_operations_with_site_tensors      → 磁对称操作（含 timerev 标记）
  → reduce_to_primitive_magsym                → 纯平移约化到原胞表示
  → msg_identify_magnetic_space_group_type    → 磁空间群 UNI + BNS
```

### Key Fixes Made

| # | Problem | Fix | Impact |
|---|---------|-----|--------|
| 1 | `spn_get_operations_with_site_tensors` 用原胞基旋转映射常规晶胞位置 | 改用 `sym_get_operation(&cell)` 常规晶胞对称操作 | FCC ops: 4→64, BCC ops: 正确 |
| 2 | 磁匹配要求 `changed_symmetry` 与 DB 操作数严格相等 | `is_subset` 替换 `is_equal`，允许子集匹配 | 磁矩破缺对称性时可匹配 DB |
| 3 | `changed_symmetry` 是常规晶胞表示，DB 用原胞表示 | `reduce_to_primitive_magsym` 用纯平移约化到原胞 | FCC 64→16 ops 正确匹配 |

### `reduce_to_primitive_magsym` (核心修复)

检测磁对称操作中的纯平移（identity + timerev=0），通过它们将操作约化到最短 translation
并去重，实现从常规晶胞到原胞的表示转换。例如 FCC 4 原子晶胞的 64 个操作含 4 个中心化
平移，约化后得到 16 个原胞操作。

### Verified Test Results (69 tests)

| Test | Structure | Non-mag | Magnetic Space Group | BNS |
|------|-----------|---------|--------------------|-----|
| Type-1/2/3 (DB ops) | Pm-3m hall 517 | #221 | UNI 1594-1598 | 221.92-221.96 |
| Fe [1,0,0] | 1 atom @ SC body center | Pm-3m (#221) | P4/mmm (#123) Type-3 | - |
| Fe [1,1,1] | 1 atom @ SC body center | Pm-3m (#221) | R-3m (#166) UNI=1331 | 166.101 |
| Co AFM [111] | 2 atoms BCC | Im-3m (#229) | R-3m (#166) UNI=1331 | 166.101 |
| FCC FM [001] | 4 atoms FCC | Fm-3m (#225) | P4/mmm (#123) UNI=1005 | 123.345 |
| FCC FM [111] | 4 atoms FCC | Fm-3m (#225) | R-3m (#166) UNI=1331 | 166.101 |
| End-to-end POSCAR | Fe SC body center [001] | Pm-3m (#221) | P4/mmm (#123) UNI=1005 | 123.345 |
| End-to-end POSCAR AFM | Fe BCC AFM [111] | Im-3m (#229) | R-3m (#166) UNI=1331 | 166.101 |

### C Code Comparison

The original C spglib has the SAME limitations and is WORSE:

| Issue | C code | Rust |
|-------|--------|------|
| Space-group search from magnetic symmetry | NULL crash (no check) | Graceful `Option::None` |
| `parent_hall_number` fallback | TODO comment only | Implemented |
| Pure translation reduction | Not implemented | `reduce_to_primitive_magsym` |
| Subset matching vs exact | Exact only | `is_subset` + `is_equal` |

C code key lines in `src_c/magnetic_spacegroup.txt`:
- Line 112-113: `/* TODO(shinohara): add option to specify hall_number ... */`
- Line 593-596: NULL dereference on failed search
- Line 615-617: Dead NULL check (never reached)

## Notes

- `debug_test.rs` and similar ad-hoc files in the repo root are gitignored/untracked — use them for quick experiments but don't commit.
- `periodic_axes` in `cell.rs` layer functions uses `[usize; 2]` through the `AperiodicAxis` enum (not raw `i32`).
