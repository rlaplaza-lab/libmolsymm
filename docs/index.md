# libmolsymm

Point-group detection, character tables, subgroups, symmetrization,
asymmetric-unit reconstruction, and symmetry-reduced finite-difference
displacements. C++20 core with C, Python, and Fortran bindings.

API reference: [Python](api/python.md), [C++](api/cpp.md), [C](api/c.md),
[Fortran](api/fortran.md). Quickstarts live in the
[README](https://github.com/rlaplaza-lab/libmolsymm#readme).

## Install

```bash
pip install libmolsymm
```

From source (CMake):

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DMOLSYMM_BUILD_PYTHON=ON -DMOLSYMM_BUILD_FORTRAN=ON
cmake --build build -j && ctest --test-dir build --output-on-failure
```

The PyPI package is `libmolsymm`; the Python import name is `molsymm`.
The CMake / pkg-config package name is `molsymm`.

## Conventions

| Topic | Rule |
| --- | --- |
| Coordinates | Ångströms |
| Z-matrix angles | Radians |
| Z-matrix indices | 1-based; `0` means unused |
| Fortran atom indices | 1-based (representatives, copies, pieces) |
| Python positions | Prefer `(n, 3)`; `(3, n)` is also accepted |
| C `xyz` | Flat length `3*n` (`x0,y0,z0,…`) |
| Fortran `xyz` | Shape `(3, n)` |
| Primary tolerance | Default `1e-3` Å (geometry / equivalence) |

Analysis places the molecule in a Cotton / principal-axis frame with the
center of mass at the origin. Displacement and reconstruct outputs are
returned in the **original input frame**.

## Errors

| Binding | Behavior |
| --- | --- |
| C++ / Python | Throws `molsymm::Error` / a Python exception |
| C | Status `0` is success; message via `molsymm_last_error()` |
| Fortran | Check `status`; message via `ms_last_error()` |
