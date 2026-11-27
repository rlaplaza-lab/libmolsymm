# Python API

Shared conventions (units, frames, indices) are on the
[home page](../index.md). Import name is `molsymm` (PyPI package
`libmolsymm`). Failures raise `RuntimeError`. Numeric results are
C-contiguous `float64` ndarrays.

```python
import molsymm
import numpy as np

a = molsymm.analyze(numbers, positions)   # -> Analysis
print(a.point_group, a.order)
rms, xyz = a.symmetrize()                 # input-frame coordinates
out = molsymm.symmetrize(numbers, positions, "C2v", forces=forces)
```

Analysis-frame helpers (`frame_xyz`, `operations`) use the Cotton /
principal-axis frame. `Analysis.symmetrize`, displacements, reconstruct,
assemble, and module-level `symmetrize` return the **original input frame**.
Fragmentation indices are 0-based.

## Functions

::: molsymm.analyze
    options:
      heading_level: 3

::: molsymm.symmetrize
    options:
      heading_level: 3

::: molsymm.analyze_zmatrix
    options:
      heading_level: 3

::: molsymm.zmatrix_to_xyz
    options:
      heading_level: 3

## Types

::: molsymm.ZMatrixRow
    options:
      heading_level: 3

## Classes

::: molsymm.Analysis
    options:
      heading_level: 3

::: molsymm.Symmetrized
    options:
      heading_level: 3

::: molsymm.DisplacementSet
    options:
      heading_level: 3

::: molsymm.Fragmentation
    options:
      heading_level: 3

::: molsymm.SubgroupAssembly
    options:
      heading_level: 3

## Examples

### Symmetrize geometry, forces, and Hessian

```python
--8<-- "tests/test_python_types.py:symmetrize_forces"
```

### Symmetry-reduced finite differences

```python
--8<-- "tests/test_python_types.py:displacements"
```

### Asymmetric unit

```python
--8<-- "tests/test_python_types.py:fragmentation"
```

### Subgroup assembly

```python
--8<-- "tests/test_python_types.py:subgroup_assembly"
```
