"""Runtime shape checks and mypy-checked call patterns for the Python API."""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt

import molsymm

if TYPE_CHECKING:
    from molsymm import ZMatrixRow

Float64Array = npt.NDArray[np.float64]


def _water() -> tuple[list[int], Float64Array]:
    numbers = [8, 1, 1]
    positions = np.array(
        [
            [0.0, 0.0, 0.1173],
            [0.0, 0.7572, -0.4692],
            [0.0, -0.7572, -0.4692],
        ],
        dtype=np.float64,
    )
    return numbers, positions


def _assert_f64(arr: object, shape: tuple[int, ...]) -> Float64Array:
    assert isinstance(arr, np.ndarray)
    assert arr.dtype == np.float64
    assert arr.shape == shape
    assert arr.flags.c_contiguous
    return arr


def test_analyze_array_returns() -> None:
    numbers, positions = _water()
    a = molsymm.analyze(numbers, positions)
    assert a.point_group == "C2v"
    assert a.order == 4
    assert a.atomic_numbers == [8, 1, 1]
    assert isinstance(a.subgroups, list)
    assert isinstance(a.irrep_names, list)

    frame = _assert_f64(a.frame_xyz, (3, 3))
    # Analysis frame splits the H atoms along x (guards column-order bugs).
    assert abs(float(frame[1, 0] - frame[2, 0])) > 0.5
    assert len(a.operations) == a.order
    for op in a.operations:
        _assert_f64(op, (3, 3))


def test_symmetrize_method_and_module() -> None:
    numbers, positions = _water()
    a = molsymm.analyze(numbers, positions)
    rms, xyz = a.symmetrize()
    assert isinstance(rms, float)
    _assert_f64(xyz, (3, 3))

    # --8<-- [start:symmetrize_forces]
    forces = np.array(
        [
            [0.0, 0.0, 0.1],
            [0.05, 0.02, -0.05],
            [0.0, -0.02, -0.05],
        ],
        dtype=np.float64,
    )
    hessian = np.eye(9, dtype=np.float64)
    out = molsymm.symmetrize(
        numbers, positions, "C2v", forces=forces, hessian=hessian, primary=1e-2
    )
    pos = out.positions
    fout = out.forces
    hout = out.hessian
    # --8<-- [end:symmetrize_forces]
    _assert_f64(pos, (3, 3))
    assert fout is not None
    _assert_f64(fout, (3, 3))
    assert hout is not None
    _assert_f64(hout, (9, 9))

    geom = molsymm.symmetrize(numbers, positions, "C2v", primary=1e-2)
    assert geom.forces is None and geom.hessian is None
    _assert_f64(geom.positions, (3, 3))

    flat_out = molsymm.symmetrize(
        numbers, positions, "C2v", forces=forces.reshape(-1), primary=1e-2
    )
    assert flat_out.forces is not None
    _assert_f64(flat_out.forces, (3, 3))


def test_displacements() -> None:
    numbers, positions = _water()
    a = molsymm.analyze(numbers, positions)

    # --8<-- [start:displacements]
    d = a.displacements(step=1e-3)
    plus: list[Float64Array] = d.plus
    minus: list[Float64Array] = d.minus
    assert d.n == len(plus) == len(minus)
    g = d.gradient_from_energies([0.0] * d.n, [0.0] * d.n)
    h = d.hessian_from_gradients(
        [np.zeros(9, dtype=np.float64) for _ in range(d.n)],
        [np.zeros(9, dtype=np.float64) for _ in range(d.n)],
    )
    # --8<-- [end:displacements]
    assert d.step == 1e-3
    for p, m in zip(plus, minus, strict=True):
        _assert_f64(p, (3, 3))
        _assert_f64(m, (3, 3))
    _assert_f64(g, (9,))
    _assert_f64(h, (9, 9))


def test_fragmentation_and_assembly() -> None:
    numbers, positions = _water()
    a = molsymm.analyze(numbers, positions)

    # --8<-- [start:fragmentation]
    frag = a.fragmentation()
    reps: list[int] = frag.representatives
    full = frag.reconstruct(positions[reps])
    # --8<-- [end:fragmentation]
    assert len(frag.free_axes) == len(reps)
    for site_axes in frag.free_axes:
        for axis in site_axes:
            assert len(axis) == 3
    _assert_f64(full, (3, 3))
    assert np.allclose(full, positions, atol=1e-6)

    # --8<-- [start:subgroup_assembly]
    asmbl = a.subgroup_assembly("C2v")
    pieces: list[Float64Array] = [positions[atoms] for atoms in asmbl.pieces]
    rebuilt = asmbl.assemble(pieces)
    # --8<-- [end:subgroup_assembly]
    assert asmbl.subgroup == "C2v"
    site = 0
    for atoms in asmbl.pieces:
        for _ in atoms:
            axes = asmbl.free_axes[site]
            if axes:
                assert len(axes[0]) == 3
            site += 1
    _assert_f64(rebuilt, (3, 3))
    assert np.allclose(rebuilt, positions, atol=1e-6)


def test_zmatrix() -> None:
    angle = 104.5 * np.pi / 180.0
    entries: Sequence[ZMatrixRow] = [
        (8,),
        (1, 1, 0.957),
        (1, 1, 0.957, 2, angle),
    ]
    _assert_f64(molsymm.zmatrix_to_xyz(entries), (3, 3))
    assert molsymm.analyze_zmatrix(entries).point_group == "C2v"


def test_select_subgroup() -> None:
    numbers, positions = _water()
    a = molsymm.analyze(numbers, positions)
    a.select_subgroup("Cs")
    assert a.point_group == "Cs" and a.order == 2


if __name__ == "__main__":
    test_analyze_array_returns()
    test_symmetrize_method_and_module()
    test_displacements()
    test_fragmentation_and_assembly()
    test_zmatrix()
    test_select_subgroup()
    print("python types ok")
