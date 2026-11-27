# C++ API

Shared conventions (units, frames, indices) are on the
[home page](../index.md).

```cpp
#include "molsymm/molsymm.hpp"
```

Umbrella header re-exports `analysis`, `displacement`, `fragment`, `molecule`,
`types`, `zmatrix`, `masses`, and `error`.

## Core types

| Type | Role |
| --- | --- |
| `Molecule` | Atoms with atomic number, position, mass; ctors from flat `xyz` or `3×N` Eigen matrix |
| `Tolerances` | `primary` (Å, default `1e-3`), `angle`, `inertia`, `max_axis_order`, `max_group_order` |
| `Analysis` | Result of `analyze`: `point_group`, `order`, `symmetry_number`, frame/original molecules, `operations`, `character_table`, `subgroups`, `rmsd`; `select_subgroup`, `symmetrize`, `frame_to_input` |
| `ZMatrixAtom` | One Z-matrix row (1-based indices; Å / radians) |

## Workflow functions

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
