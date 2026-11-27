# API reference

## Python

```python
import molsymm
import numpy as np
```

### Module functions

| Function | Description |
| --- | --- |
| `analyze(numbers, positions, primary=1e-3)` | Detect point group → `Analysis`. `numbers` is atomic numbers; `positions` is `(n,3)` or `(3,n)`. |
| `symmetrize(numbers, positions, point_group, forces=None, hessian=None, primary=1e-3)` | Project geometry and optional forces / Hessian onto a named group → `Symmetrized`. |
| `analyze_zmatrix(entries, primary=1e-3)` | Build molecule from Z-matrix rows and analyze → `Analysis`. |
| `zmatrix_to_xyz(entries)` | Z-matrix → `(n,3)` Cartesian positions. |

Z-matrix `entries` are tuples of length 1, 3, 5, or 7:

`(Z)`, `(Z, bond_to, R)`, `(Z, bond_to, R, angle_to, angle)`,
`(Z, bond_to, R, angle_to, angle, dihedral_to, dihedral)`.

### `Analysis`

| Member | Description |
| --- | --- |
| `point_group`, `order`, `symmetry_number`, `rmsd`, `subgroups` | Group metadata (`order` is `0` for infinite groups) |
| `atomic_numbers` | Atomic numbers in the analysis frame |
| `frame_xyz` | Coordinates in the analysis frame |
| `operations` | List of 3×3 operation matrices |
| `irrep_names` | Irreducible representation labels |
| `select_subgroup(name)` | Restrict to a named subgroup |
| `symmetrize()` | → `(rms, xyz)` in the input frame |
| `displacements(step=1e-3)` | → `DisplacementSet` |
| `fragmentation()` | → `Fragmentation` |
| `subgroup_assembly(name)` | → `SubgroupAssembly` |

### `DisplacementSet`

| Member | Description |
| --- | --- |
| `step`, `n` | Finite-difference step and unique count |
| `plus`, `minus` | Lists of `(n,3)` geometries in the input frame |
| `gradient_from_energies(ep, em)` | Central-difference gradient from ± energies |
| `hessian_from_gradients(gp, gm)` | Hessian from ± Cartesian gradients (each length `3n`) |

### `Fragmentation` / `SubgroupAssembly` / `Symmetrized`

| Type | Members |
| --- | --- |
| `Fragmentation` | `representatives`, `copies`, `mask`, `free_axes`; `reconstruct(positions)` |
| `SubgroupAssembly` | `subgroup`, `pieces`, `free_axes`; `assemble(pieces)` |
| `Symmetrized` | `rmsd`, `positions`, `forces`, `hessian` (`forces` / `hessian` are `None` when omitted) |

```python
a = molsymm.analyze([8, 1, 1], positions)
a.select_subgroup("Cs")
rms, xyz = a.symmetrize()
d = a.displacements(step=1e-3)
frag = a.fragmentation()
xyz = frag.reconstruct(positions[frag.representatives])
```

---

## C++

```cpp
#include "molsymm/molsymm.hpp"
```

Umbrella header re-exports `analysis`, `displacement`, `fragment`, `molecule`,
`types`, `zmatrix`, `masses`, and `error`.

### Core types

| Type | Role |
| --- | --- |
| `Molecule` | Atoms with atomic number, position, mass; ctors from flat `xyz` or `3×N` Eigen matrix |
| `Tolerances` | `primary` (Å, default `1e-3`), `angle`, `inertia`, `max_axis_order`, `max_group_order` |
| `Analysis` | Result of `analyze`: `point_group`, `order`, `symmetry_number`, frame/original molecules, `operations`, `character_table`, `subgroups`, `rmsd`; `select_subgroup`, `symmetrize`, `frame_to_input` |
| `ZMatrixAtom` | One Z-matrix row (1-based indices; Å / radians) |

### Workflow functions

| Function | Description |
| --- | --- |
| `analyze(mol, tol={})` | Detect point group → `Analysis` |
| `molecule_from_zmatrix(rows)` / `zmatrix_positions(rows)` | Z-matrix → `Molecule` / positions |
| `symmetrize(mol, point_group, forces*, hessian*, tol={})` | One-shot project → `Symmetrized` |
| `unique_displacements(analysis, step=1e-3)` | → `DisplacementSet` with `plus` / `minus` |
| `fragmentation(analysis)` / `reconstruct(frag, reps)` | Asymmetric unit ↔ full molecule |
| `subgroup_assembly(analysis, name)` / `assemble(asm, pieces)` | Lower-symmetry assembly |

`DisplacementSet` provides `gradient_from_energies` and `hessian_from_gradients`.
`Symmetrized` holds `rmsd`, `positions` (input frame), optional `forces` / `hessian`.

C++-only helpers (not wrapped in other bindings): `symmetrize_gradient`,
`symmetrize_hessian`, `project_external`, `center_of_mass`, `inertia_tensor`,
`atomic_mass`.

```cpp
molsymm::Analysis a = molsymm::analyze(molsymm::Molecule(/*Z*/, /*xyz flat*/));
a.select_subgroup("C2v");
a.symmetrize();
auto d = molsymm::unique_displacements(a, 1e-3);
auto frag = molsymm::fragmentation(a);
molsymm::Molecule full = molsymm::reconstruct(frag, reps);
```

---

## C

```c
#include "molsymm/molsymm.h"
```

Opaque handles: `molsymm_analysis*`, `molsymm_displacements*`,
`molsymm_fragmentation*`, `molsymm_assembly*`. Status `0` is success;
`molsymm_last_error()` returns the last message. Free with the matching
`molsymm_free*` / `molsymm_free`.

### Analyze and query

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

### Displacements

| Function | Description |
| --- | --- |
| `molsymm_unique_displacements(a, step, &out)` | Build set |
| `molsymm_n_displacements` / `molsymm_displacement_atom` / `molsymm_displacement_direction` / `molsymm_displacement_preserved` | Metadata |
| `molsymm_displacement_geometry(d, index, plus, xyz_out)` | `plus≠0` → +, else − |
| `molsymm_gradient_from_energies` / `molsymm_hessian_from_gradients` | FD assembly |

### Fragmentation and assembly

| Function | Description |
| --- | --- |
| `molsymm_make_fragmentation` / representatives, copies, mask, free axes / `molsymm_reconstruct` | Asymmetric unit |
| `molsymm_subgroup_assembly` / pieces / free axes / `molsymm_assemble` | Subgroup rebuild |

`rep_xyz` for reconstruct is `3 * n_rep`. `piece_xyz` for assemble is the
concatenated piece coordinates.

---

## Fortran

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
| `ms_point_group_name` / `ms_symmetry_number` / `ms_order` | Metadata |
| `ms_symmetrize(handle, xyz_out, rmsd, status)` | Symmetrize |
| `ms_symmetrize_data(..., status [, forces, forces_out, hessian, hessian_out, primary])` | One-shot |
| `ms_unique_displacements` / `ms_n_displacements` / `ms_displacement_geometry` / `ms_hessian_from_gradients` / `ms_release_displacements` | Displacements |
| `ms_fragmentation` / `ms_n_representatives` / `ms_representative_atom` / `ms_n_copies` / `ms_copy_atoms` / `ms_fragment_mask` / `ms_reconstruct` / `ms_release_fragmentation` | Fragmentation |
| `ms_subgroup_assembly` / `ms_assembly_n_pieces` / `ms_assembly_piece_atoms` / `ms_assemble` / `ms_release_assembly` | Assembly |

`ms_displacement_geometry(..., plus, ...)`: `1` = plus geometry, `0` = minus.
`rep_xyz` is `(3, nrep)`; `piece_xyz` is `(3, nasym)`.
