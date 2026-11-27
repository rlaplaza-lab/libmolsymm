# C API

Shared conventions (units, frames, indices) are on the
[home page](../index.md).

```c
#include "molsymm/molsymm.h"
```

Opaque handles: `molsymm_analysis*`, `molsymm_displacements*`,
`molsymm_fragmentation*`, `molsymm_assembly*`. Status `0` is success;
`molsymm_last_error()` returns the last message. Free with the matching
`molsymm_free*` / `molsymm_free`.

## Analyze and query

| Function | Description |
| --- | --- |
| `molsymm_analyze(z, xyz, n, primary, &out)` | Analyze Cartesian geometry |
| `molsymm_analyze_zmatrix(rows, n, primary, &out)` | Analyze Z-matrix |
| `molsymm_zmatrix_to_xyz(rows, n, xyz_out)` | Z-matrix → flat `xyz` |
| `molsymm_point_group` / `molsymm_order` / `molsymm_symmetry_number` / `molsymm_rmsd` | Metadata |
| `molsymm_n_atoms` / `molsymm_n_operations` / `molsymm_operation_matrix` | Geometry / ops |
| `molsymm_n_irreps` / `molsymm_irrep_name` | Character table labels |
| `molsymm_n_subgroups` / `molsymm_subgroup_name` / `molsymm_select_subgroup` | Subgroups |
| `molsymm_symmetrize(a, &rmsd, xyz_out)` | Symmetrize analysis geometry |
| `molsymm_symmetrize_data(...)` | One-shot geometry ± forces ± Hessian |

`molsymm_symmetrize_data`: `xyz` / `forces` length `3*n`; `hessian` row-major
`3n×3n`. Pass `NULL` to skip optional arrays.

## Displacements

| Function | Description |
| --- | --- |
| `molsymm_unique_displacements(a, step, &out)` | Build set |
| `molsymm_n_displacements` / `molsymm_displacement_atom` / `molsymm_displacement_direction` / `molsymm_displacement_preserved` | Metadata |
| `molsymm_displacement_geometry(d, index, plus, xyz_out)` | `plus≠0` → +, else − |
| `molsymm_gradient_from_energies` / `molsymm_hessian_from_gradients` | FD assembly |

## Fragmentation and assembly

| Function | Description |
| --- | --- |
| `molsymm_make_fragmentation` / representatives, copies, mask, free axes / `molsymm_reconstruct` | Asymmetric unit |
| `molsymm_subgroup_assembly` / pieces / free axes / `molsymm_assemble` | Subgroup rebuild |

`rep_xyz` for reconstruct is `3 * n_rep`. `piece_xyz` for assemble is the
concatenated piece coordinates.
