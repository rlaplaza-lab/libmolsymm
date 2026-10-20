# libmolsymm

C++20 library for molecular point-group symmetry: detection, character tables,
subgroups, geometry symmetrization, and **symmetry-reduced finite-difference
displacements** that rebuild Cartesian gradients and Hessians.

Independent implementation. Detection ideas are related to
Johansson, M., Veryazov, V. Automatic procedure for generating symmetry
adapted wavefunctions. *J Cheminform* **9**, 8 (2017).
https://doi.org/10.1186/s13321-017-0193-3
Orbital SALCs are out of scope.

## Features

- Point groups: `C1`, `Cs`, `Ci`, `Cn`/`Cnv`/`Cnh`, `Sn`, `Dn`/`Dnh`/`Dnd`,
  `T`/`Td`/`Th`, `O`/`Oh`, `I`/`Ih`, linear `C∞v`/`D∞h`, single atom `Kh`
- Real Mulliken character tables
- Constructive subgroup lists and `select_subgroup`
- Geometry symmetrization with RMS residual
- Unique Cartesian displacement orbits for numerical derivatives
- Z-matrix geometries in the Priroda / briling placement, converted to Cartesians
- Native APIs: C++, C, Python (nanobind), Fortran 2003

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Optional:

```bash
cmake -S . -B build -DMOLSYMM_BUILD_PYTHON=ON -DMOLSYMM_BUILD_FORTRAN=ON
cmake --build build -j
# Python module is build/molsymm*.so — add build/ to PYTHONPATH
```

## C++ quick start

```cpp
#include "molsymm/molsymm.hpp"

molsymm::Molecule mol(/* Z */, /* xyz flat, Å */);
molsymm::Analysis a = molsymm::analyze(mol);
std::cout << a.point_group << " σ=" << a.symmetry_number << "\n";

auto disp = molsymm::unique_displacements(a, 1e-3);
// evaluate energies at disp.plus[i] / disp.minus[i]
// auto g = disp.gradient_from_energies(e_plus, e_minus);
// auto H = disp.hessian_from_gradients(g_plus, g_minus);
```

## Python

```bash
pip install libmolsymm
```

```python
import molsymm
import numpy as np

a = molsymm.analyze([8, 1, 1], np.array([[0, 0, 0.1173],
                                         [0, 0.7572, -0.4692],
                                         [0, -0.7572, -0.4692]]))
print(a.point_group, a.symmetry_number)
d = a.displacements(step=1e-3)
print(d.n)

angle = 104.5 * np.pi / 180.0
z = molsymm.analyze_zmatrix([
    (8,),
    (1, 1, 0.957),
    (1, 1, 0.957, 2, angle),
])
print(z.point_group)
```

From a source checkout (needs CMake ≥ 3.20 and a C++20 compiler):

```bash
pip install .
```

## Fortran

`xyz` is `(3, n)`.

```fortran
use molsymm
call ms_analyze(z, xyz, handle, status)
call ms_point_group_name(handle, name, status)
call ms_unique_displacements(handle, 1.0d-3, dhandle, status)
```

## Z-matrix

A Z-matrix is a second way to build the Cartesian geometry. Detection,
symmetrization, and displacements are unchanged and still return Cartesians.
Indices are **1-based** (0 means unused). Distances are **ångströms** and angles
are **radians**, the units `zmat2cart` uses after a Priroda file's degree values
are converted.

The first three atoms follow the [briling/v](https://github.com/briling/v/blob/master/src/math/zmat.c)
special cases: atom 1 at the origin, atom 2 at `reference + (0, R, 0)`, atom 3 in
the plane `z = reference_z`. Later atoms use a bond, a bend, and a dihedral.
References must point at earlier atoms. Rows with `atomic_number <= 0` are kept
and receive the usual mass fallback of 1, so dummy centers should be omitted
before analysis.

```cpp
constexpr double deg = 0.017453292519943295; // π/180
std::vector<molsymm::ZMatrixAtom> z = {
    {8, 0, 0, 0, 0.0, 0.0, 0.0},
    {1, 1, 0, 0, 0.957, 0.0, 0.0},
    {1, 1, 2, 0, 0.957, 104.5 * deg, 0.0},
};
molsymm::Analysis a = molsymm::analyze(molsymm::molecule_from_zmatrix(z));
```

`zmatrix_positions` returns the Cartesian coordinates alone. The C row has the
same field order: `atomic_number`, `bond_to`, `angle_to`, `dihedral_to`,
`distance`, `angle`, `dihedral`.

```c
molsymm_zatom z[] = {
    {8, 0, 0, 0, 0.0, 0.0, 0.0},
    {1, 1, 0, 0, 0.957, 0.0, 0.0},
    {1, 1, 2, 0, 0.957, 104.5 * 0.017453292519943295, 0.0},
};
molsymm_analysis* a = NULL;
molsymm_analyze_zmatrix(z, 3, 1e-3, &a);
double xyz[9];
molsymm_zmatrix_to_xyz(z, 3, xyz);
```

Python entries are `(Z,)`, `(Z, bond_to, R)`, `(Z, bond_to, R, angle_to, angle)`,
or the full 7-tuple. `analyze_zmatrix` is shown above; `zmatrix_to_xyz(entries)`
returns `(n, 3)`.

```fortran
type(ms_zatom) :: z(3)
real(c_double) :: ang, xyz(3, 3)
ang = 104.5d0 * acos(-1.0d0) / 180.0d0
z = ms_zatom(0, 0, 0, 0, 0.0d0, 0.0d0, 0.0d0)
z(1)%atomic_number = 8
z(2)%atomic_number = 1
z(2)%bond_to = 1
z(2)%distance = 0.957d0
z(3)%atomic_number = 1
z(3)%bond_to = 1
z(3)%angle_to = 2
z(3)%distance = 0.957d0
z(3)%angle = ang
call ms_analyze_zmatrix(z, handle, status)
call ms_zmatrix_to_xyz(z, xyz, status)  ! xyz is (3, n)
```

## Finite-difference contract

`unique_displacements` returns one representative per orbit of
`(atom, Cartesian direction)` under the point group. Water (`C2v`) has **6**
unique coordinates of 9; methane (`Td`) has **3** of 15. Each entry includes
the preserved subgroup of the displaced geometry (for QC symmetry on the
displaced point).

- `gradient_from_energies` — central differences on unique coordinates, scattered with `R`
- `hessian_from_gradients` — one Hessian column per unique step, completed by group action
- `symmetrize_gradient` / `symmetrize_hessian` — project a noisy full tensor
- `project_external` — remove translation/rotation after assembly

## License

MIT
