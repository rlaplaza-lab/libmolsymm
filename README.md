# libmolsymm

[![PyPI](https://img.shields.io/pypi/v/libmolsymm.svg)](https://pypi.org/project/libmolsymm/)

Point-group detection, character tables, subgroups, symmetrization,
asymmetric-unit reconstruction, and symmetry-reduced finite-difference
displacements. C++20 core with C, Python, and Fortran bindings. Coordinates
are ångströms; Z-matrix angles are radians (1-based indices, `0` unused).

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DMOLSYMM_BUILD_PYTHON=ON -DMOLSYMM_BUILD_FORTRAN=ON
cmake --build build -j && ctest --test-dir build --output-on-failure
# or: pip install libmolsymm
```

## C++

```cpp
#include "molsymm/molsymm.hpp"

molsymm::Analysis a = molsymm::analyze(molsymm::Molecule(/*Z*/, /*xyz flat*/));
// a.point_group, a.order, a.symmetry_number, a.subgroups, a.operations, …
a.select_subgroup("C2v");
a.symmetrize();  // updates frame geometry; RMS in a.rmsd

auto d = molsymm::unique_displacements(a, 1e-3);
// d.plus / d.minus → evaluate energies or gradients
// d.gradient_from_energies(ep, em);
// d.hessian_from_gradients(gp, gm);

auto frag = molsymm::fragmentation(a);
// frag.representatives, frag.mask, frag.free_axes (original input frame)
molsymm::Molecule full = molsymm::reconstruct(frag, reps);
auto asmbl = molsymm::subgroup_assembly(a, "C2");
molsymm::Molecule lower = molsymm::assemble(asmbl, pieces);

molsymm::Analysis z = molsymm::analyze(
    molsymm::molecule_from_zmatrix({
        {8, 0, 0, 0, 0.0, 0.0, 0.0},
        {1, 1, 0, 0, 0.957, 0.0, 0.0},
        {1, 1, 2, 0, 0.957, 104.5 * 0.017453292519943295, 0.0},
    }));
```

## C

`#include "molsymm/molsymm.h"` — opaque `molsymm_analysis*` / `molsymm_displacements*` /
`molsymm_fragmentation*` / `molsymm_assembly*`. `xyz` is length `3*n` (x₀,y₀,z₀,…).
Status `0` is success; see `molsymm_last_error()`.

```c
molsymm_analysis* a = NULL;
molsymm_analyze(z, xyz, n, 1e-3, &a);
char pg[32]; int sigma, order;
molsymm_point_group(a, pg, 32);
molsymm_symmetry_number(a, &sigma);
molsymm_order(a, &order);
molsymm_select_subgroup(a, "C2v");
double rmsd, xyz_out[/*3*n*/];
molsymm_symmetrize(a, &rmsd, xyz_out);

molsymm_displacements* d = NULL;
molsymm_unique_displacements(a, 1e-3, &d);
molsymm_hessian_from_gradients(d, g_plus, g_minus, hess_out);
molsymm_free_displacements(d);

molsymm_fragmentation* f = NULL;
molsymm_make_fragmentation(a, &f);
molsymm_reconstruct(f, rep_xyz, xyz_out);  /* rep_xyz is 3*n_rep */
molsymm_free_fragmentation(f);
molsymm_assembly* asm_ = NULL;
molsymm_subgroup_assembly(a, "C2", &asm_);
molsymm_assemble(asm_, piece_xyz, xyz_out);  /* piece_xyz concatenated */
molsymm_free_assembly(asm_);
molsymm_free(a);

molsymm_zatom rows[] = {
    {8, 0, 0, 0, 0.0, 0.0, 0.0},
    {1, 1, 0, 0, 0.957, 0.0, 0.0},
    {1, 1, 2, 0, 0.957, 1.8325957, 0.0},
};
molsymm_analyze_zmatrix(rows, 3, 1e-3, &a);
molsymm_zmatrix_to_xyz(rows, 3, xyz_out);
```

## Python

```python
import molsymm
import numpy as np

positions = np.array([[0.0, 0.0, 0.1173],
                      [0.0, 0.7572, -0.4692],
                      [0.0, -0.7572, -0.4692]])
a = molsymm.analyze([8, 1, 1], positions)
# a.point_group, a.order, a.symmetry_number, a.subgroups
# a.operations, a.irrep_names, a.atomic_numbers, a.frame_xyz
a.select_subgroup("Cs")
rms, xyz = a.symmetrize()
d = a.displacements(step=1e-3)  # d.n, d.plus, d.minus
# d.gradient_from_energies(ep, em); d.hessian_from_gradients(gp, gm)

frag = a.fragmentation()  # .representatives, .mask, .copies, .free_axes
xyz = frag.reconstruct(positions[frag.representatives])
asmbl = a.subgroup_assembly("C2")
xyz = asmbl.assemble([positions[atoms] for atoms in asmbl.pieces])

ang = 104.5 * np.pi / 180.0
a = molsymm.analyze_zmatrix([(8,), (1, 1, 0.957), (1, 1, 0.957, 2, ang)])
xyz = molsymm.zmatrix_to_xyz([(8,), (1, 1, 0.957), (1, 1, 0.957, 2, ang)])
```

## Fortran

`use molsymm`. Cartesians are `(3, n)`. Check `status` / `ms_last_error()`.

```fortran
call ms_analyze(z, xyz, handle, status)           ! optional primary=
call ms_point_group_name(handle, name, status)
call ms_symmetry_number(handle, sigma, status)
call ms_order(handle, order, status)
call ms_symmetrize(handle, xyz_out, rmsd, status)
call ms_unique_displacements(handle, 1.0d-3, dhandle, status)
call ms_n_displacements(dhandle, nd, status)
call ms_displacement_geometry(dhandle, i, 1, xyz_p, status)  ! 1=+ / 0=-
call ms_hessian_from_gradients(dhandle, gp, gm, hess, status)
call ms_release_displacements(dhandle)
call ms_fragmentation(handle, fhandle, status)
call ms_n_representatives(fhandle, nrep, status)
call ms_reconstruct(fhandle, rep_xyz, xyz_out, status)  ! rep_xyz is (3,nrep)
call ms_release_fragmentation(fhandle)
call ms_subgroup_assembly(handle, "C2", ahandle, status)
call ms_assemble(ahandle, piece_xyz, xyz_out, status)   ! piece_xyz is (3,nasym)
call ms_release_assembly(ahandle)
call ms_release(handle)

type(ms_zatom) :: rows(3)
! fill atomic_number, bond_to, distance, angle_to, angle, …
call ms_analyze_zmatrix(rows, handle, status)
call ms_zmatrix_to_xyz(rows, xyz, status)
```

## License

MIT
