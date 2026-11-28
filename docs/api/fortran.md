# Fortran API

Shared conventions (units, frames, indices) are on the
[home page](../index.md).

```fortran
use molsymm
```

Cartesians are `(3, n)`. Check `status` after each call; use `ms_last_error()`
on failure. Type `ms_zatom` mirrors the C Z-matrix row.

| Procedure | Description |
| --- | --- |
| `ms_analyze(z, xyz, handle, status [, primary])` | Analyze XYZ |
| `ms_analyze_zmatrix(rows, handle, status [, primary])` | Analyze Z-matrix |
| `ms_zmatrix_to_xyz(rows, xyz, status)` | Z-matrix → XYZ |
| `ms_release(handle)` | Free analysis |
| `ms_point_group_name` / `ms_symmetry_number` / `ms_order` / `ms_rmsd` | Metadata |
| `ms_symmetrize(handle, xyz_out, rmsd, status)` | Symmetrize |
| `ms_symmetrize_data(..., status [, forces, forces_out, hessian, hessian_out, primary])` | One-shot |
| `ms_unique_displacements` / `ms_n_displacements` / `ms_displacement_geometry` / `ms_hessian_from_gradients` / `ms_release_displacements` | Displacements |
| `ms_fragmentation` / `ms_n_representatives` / `ms_representative_atom` / `ms_n_copies` / `ms_copy_atoms` / `ms_fragment_mask` / `ms_reconstruct` / `ms_release_fragmentation` | Fragmentation |
| `ms_subgroup_assembly` / `ms_assembly_n_pieces` / `ms_assembly_piece_atoms` / `ms_assemble` / `ms_release_assembly` | Assembly |
| `ms_wyckoff_positions` / `ms_wyckoff_n_sites` / `ms_wyckoff_site_label` / multiplicity / dof / `ms_release_wyckoff_table` | Special positions |
| `ms_wyckoff_plans_nsites` / `ms_wyckoff_n_plans` / orbit label accessors / `ms_release_wyckoff_plans` | Orbit plans (`n_sites`) |
| `ms_assemble_from_au(pg, z, au_xyz, labels, handle, status)` | Group + AU → analysis |

`ms_displacement_geometry(..., plus, ...)`: `1` = plus geometry, `0` = minus.
`rep_xyz` is `(3, nrep)`; `piece_xyz` is `(3, nasym)`.
`au_xyz` for `ms_assemble_from_au` is `(3, n_orbits)`; Wyckoff site indices in
table/plan accessors are **1-based**.
Representative, copy, and piece atom indices are **1-based** (matching Z-matrix
indices). C and Python use 0-based indices.
