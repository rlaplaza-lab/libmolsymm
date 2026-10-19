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
```

From a source checkout (needs CMake ≥ 3.20 and a C++20 compiler):

```bash
pip install .
```

## Fortran

```fortran
use molsymm
call ms_analyze(z, xyz, handle, status)
call ms_point_group_name(handle, name, status)
call ms_unique_displacements(handle, 1.0d-3, dhandle, status)
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
