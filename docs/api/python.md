# Python API

Import name is `molsymm` (PyPI package `libmolsymm`). Units, frames, and index
bases are on the [home page](../index.md). Analysis helpers return coordinates
in the analysis frame; displacement / reconstruct / assemble / one-shot
`symmetrize` return the original input frame. Fragmentation indices are
0-based. Failures raise `RuntimeError`. Numeric results are C-contiguous
`float64` ndarrays.

```python
import molsymm
import numpy as np
```

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

## Reference

::: molsymm
