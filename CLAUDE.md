# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build & Test Commands

```bash
cargo build              # Build the library
cargo test               # Run all 67 tests
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

Magnetic space group identification is fully functional for configurations where the magnetic symmetry retains all parent space group operations (with time-reversal flags).

### Implemented Components

| Component | Lines | Description |
|-----------|-------|-------------|
| `convert_msgdb.py` | ~200 | Python script to auto-convert C data arrays to Rust static data |
| `msg_database_gen.rs` | **28,319** | Generated static arrays (6 arrays: 1652 MSG types, 531 Hall mappings, 76683 encoded operations, etc.) |
| `msg_database.rs` | ~120 | Database query functions |
| `magnetic_spacegroup.rs` | ~500 | Magnetic identification pipeline: FSG/XSG extraction, type classification (1-4), database matching |
| `magnetic_spacegroup_test.rs` | ~200 | 10 tests covering Type-1/2/3 identification via Pm-3m |
| `lib.rs` public API | ~60 | `SpglibMagneticSpacegroupType`, `spg_get_magnetic_spacegroup_type`, `spg_get_magnetic_spacegroup_type_from_symmetry` |

### Database Architecture

The magnetic space group database stores operations using a three-level indexing scheme:

```
HALL_MAPPING[hall_number] → [min_uni, max_uni]        # Hall → UNI range
UNI_MAPPING[uni_number]   → [num_halls, first_hall]    # UNI → Hall mapping
OPERATION_INDEX[uni][offset] → [count, start]           # UNI → operation index
```

Each operation in `MAGNETIC_SYMMETRY_OPERATIONS` is a single i32 encoding:
- Rotation matrix (base-3, 3^9 = 19683 values)
- Translation vector (base-12, 12^3 = 1728 values)
- Time-reversal flag (34012224 = 3^9 × 12^3, timerev = encoded / 34012224)

**Critical: every MSG entry stores ALL parent space group operations** with timerev flags.
No entry ever stores a reduced operation set. For example, all Type-3 variants of Pm-3m (#221) under Hall 517 (UNI 1596-1598) have 96 encoded operations
(= 48 rotations × 2 for centering), with 48 timerev=0 and 48 timerev=1.

### Verified Magnetic Types

| Type | Description | Example | UNI |
|------|-------------|---------|-----|
| 1 | Ordinary (no time reversal) | Pm-3m (#221) | 1524 / 1594 |
| 2 | Grey / paramagnetic (pure 1') | Pm-3m (#221) | 1525 / 1595 |
| 3 | Black-white (proper/improper) | Pm-3m (#221) | 1526 / 1596-1598 |

### Known Limitation: Symmetry-Reduced Magnetic Configurations

The standard magnetic identification pipeline (`msg_identify_magnetic_space_group_type`) only works for
magnetic configurations where ALL parent symmetry operations are either ordinary (timerev=0) or
anti (timerev=1). It cannot handle configurations where the magnetic moment reduces the number
of valid symmetry operations below the parent space group's full set.

**Concrete example**: Fe at simple cubic body center [0.5,0.5,0.5] with moment along [1,0,0]:
- Parent (non-magnetic): Pm-3m (#221), 48 symmetry operations
- Moment [1,0,0]: only 16 of the 48 ops preserve or reverse the moment (8 ordinary, 8 anti)
- The remaining 32 ops are not symmetries of the magnetic structure at all
- The pipeline receives only 16 ops → FSG=16, XSG=8 → Type-3 classification
- But no database entry under Hall 517 has 16 ops (all entries use 48 or 96)
- Result: `None` (no matching magnetic space group)

This is NOT a code bug — it is a fundamental limitation of the Litvin magnetic space group
database and the spglib identification algorithm.

#### C Code Comparison

The original C spglib (`src_c/magnetic_spacegroup.txt`) has the **same limitation** and is
**strictly worse**:

1. **C line 112-113**: Explicit TODO that was never implemented:
   ```c
   /* TODO(shinohara): add option to specify hall_number in searching
    * space-group type */
   ```

2. **C line 593-596**: NULL pointer dereference (crash) on reduced symmetry:
   ```c
   *spacegroup = spa_search_spacegroup_with_symmetry(prim_sym, unit_lat, symprec);
   ref_find_similar_bravais_lattice(*spacegroup, symprec);  // CRASH when NULL
   ```

3. **C line 615-617**: NULL check is too late, never reached:
   ```c
   if (spacegroup == NULL) { return NULL; }  // dead code
   ```

The Rust implementation improves on C by:
- Returning `Option` (graceful `None` instead of crash)
- Providing `msg_identify_with_parent_hall` fallback that uses non-magnetic parent Hall number
- Proper error propagation throughout the pipeline

#### To Fully Support Arbitrary Magnetic Moments

Handling [1,0,0] / [1,1,1] Fe correctly requires re-indexing the crystal in the lower crystal
system (tetragonal or trigonal) determined by the magnetic symmetry, then looking up magnetic
variants of that lower-symmetry parent. This is an extension beyond what C spglib provides,
and would require:

1. Detecting the effective point group from the reduced magnetic symmetry
2. Computing a transformation matrix to the standard orientation of that crystal system
3. Transforming the lattice and atomic positions accordingly
4. Re-running space group identification on the transformed structure
5. Looking up magnetic variants of the identified lower-symmetry space group

None of the existing tests cover this case — `test_fe_center_cubic_magnetic` verifies the
symmetry analysis counts (8 ordinary + 8 anti for [1,0,0], 6+6 for [1,1,1]) but
gracefully accepts that no database match exists.

### API: `msg_identify_with_parent_hall`

```rust
pub fn msg_identify_with_parent_hall(
    lattice: &Mat3,
    magnetic_symmetry: &MagneticSymmetry,
    parent_hall_number: Option<i32>,  // Non-magnetic parent Hall number
    symprec: f64,
) -> Option<MagneticDataset>
```

When `parent_hall_number` is `Some(hall)`, the pipeline builds a reference spacegroup from
that Hall number directly (bypassing the space group search from magnetic symmetry).
This is useful when the non-magnetic parent is known but the magnetic symmetry is reduced.
The fallback still uses FSG/XSG extraction from the magnetic symmetry for type classification.

## Notes

- `debug_test.rs` and similar ad-hoc files in the repo root are gitignored/untracked — use them for quick experiments but don't commit.
- `periodic_axes` in `cell.rs` layer functions uses `[usize; 2]` through the `AperiodicAxis` enum (not raw `i32`).
