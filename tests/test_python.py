"""Python tests for molsymm against the IQmol manifest corpus."""

from __future__ import annotations

import os
from collections.abc import Sequence
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

import molsymm

if TYPE_CHECKING:
    from molsymm import ZMatrixRow

DATA_DIR = Path(os.environ.get("MOLSYMM_TEST_DATA_DIR", Path(__file__).parent / "data"))
MANIFEST = DATA_DIR / "manifest.tsv"

_SYMBOLS = {
    "H": 1,
    "He": 2,
    "Li": 3,
    "Be": 4,
    "B": 5,
    "C": 6,
    "N": 7,
    "O": 8,
    "F": 9,
    "Ne": 10,
    "Na": 11,
    "Mg": 12,
    "Al": 13,
    "Si": 14,
    "P": 15,
    "S": 16,
    "Cl": 17,
    "Ar": 18,
    "K": 19,
    "Ca": 20,
    "Sc": 21,
    "Ti": 22,
    "V": 23,
    "Cr": 24,
    "Mn": 25,
    "Fe": 26,
    "Co": 27,
    "Ni": 28,
    "Cu": 29,
    "Zn": 30,
    "Ga": 31,
    "Ge": 32,
    "As": 33,
    "Se": 34,
    "Br": 35,
    "Kr": 36,
    "Rb": 37,
    "Sr": 38,
    "Y": 39,
    "Zr": 40,
    "Nb": 41,
    "Mo": 42,
    "Tc": 43,
    "Ru": 44,
    "Rh": 45,
    "Pd": 46,
    "Ag": 47,
    "Cd": 48,
    "In": 49,
    "Sn": 50,
    "Sb": 51,
    "Te": 52,
    "I": 53,
    "Xe": 54,
    "Cs": 55,
    "Ba": 56,
    "La": 57,
    "Hf": 72,
    "Ta": 73,
    "W": 74,
    "Re": 75,
    "Os": 76,
    "Ir": 77,
    "Pt": 78,
    "Au": 79,
    "Hg": 80,
}


def _atomic_number(symbol: str) -> int:
    if symbol not in _SYMBOLS:
        raise ValueError(f"unknown element symbol: {symbol}")
    return _SYMBOLS[symbol]


def load_xyz(path: Path) -> tuple[list[int], np.ndarray]:
    lines = path.read_text().splitlines()
    nat = int(lines[0].strip())
    numbers: list[int] = []
    coords: list[list[float]] = []
    for line in lines[2 : 2 + nat]:
        parts = line.split()
        numbers.append(_atomic_number(parts[0]))
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return numbers, np.asarray(coords, dtype=float)


def load_manifest() -> list[tuple[str, str, int, int, float, float, str]]:
    rows = []
    for line in MANIFEST.read_text().splitlines()[1:]:
        if not line.strip():
            continue
        path, group, order, sigma, primary, max_rmsd, source = line.split("\t")
        rows.append((path, group, int(order), int(sigma), float(primary), float(max_rmsd), source))
    return rows


def test_synthetic_water() -> None:
    numbers = [8, 1, 1]
    positions = np.array(
        [
            [0.0, 0.0, 0.1173],
            [0.0, 0.7572, -0.4692],
            [0.0, -0.7572, -0.4692],
        ]
    )
    a = molsymm.analyze(numbers, positions)
    assert a.point_group == "C2v"
    assert a.symmetry_number == 2
    assert a.order == 4
    assert set(a.irrep_names) >= {"A1", "A2", "B1", "B2"}
    d = a.displacements()
    assert d.n == 6
    g_rec = np.asarray(d.gradient_from_energies(np.zeros(d.n), np.zeros(d.n)), dtype=np.float64)
    assert g_rec.shape == (9,) and float(np.linalg.norm(g_rec)) < 1e-8
    zeros = [np.zeros(9) for _ in range(d.n)]
    assert np.asarray(d.hessian_from_gradients(zeros, zeros), dtype=float).shape == (9, 9)

    rattled = positions.copy()
    rattled[0, 0] += 0.02
    rmsd, sym_xyz = molsymm.analyze(numbers, rattled, primary=1e-2).symmetrize()
    assert rmsd >= 0.0
    assert molsymm.analyze(numbers, sym_xyz, primary=1e-2).point_group == "C2v"

    frag = a.fragmentation()
    assert frag.representatives == [0, 1] and frag.mask == [1, 1, 0]
    assert np.allclose(frag.reconstruct(positions[frag.representatives]), positions, atol=1e-6)
    asmbl = a.subgroup_assembly("C2v")
    assert len(asmbl.pieces) == 1
    pieces = [positions[atoms] for atoms in asmbl.pieces]
    assert np.allclose(asmbl.assemble(pieces), positions, atol=1e-6)


def _synthetic_benzene() -> tuple[list[int], np.ndarray]:
    zs: list[int] = []
    xyz: list[list[float]] = []
    for i in range(6):
        ang = 2.0 * np.pi * i / 6.0
        zs.append(6)
        xyz.append([1.39 * np.cos(ang), 1.39 * np.sin(ang), 0.0])
    for i in range(6):
        ang = 2.0 * np.pi * i / 6.0
        zs.append(1)
        xyz.append([2.47 * np.cos(ang), 2.47 * np.sin(ang), 0.0])
    return zs, np.asarray(xyz, dtype=float)


def test_synthetic_benzene_subgroups() -> None:
    zs, xyz = _synthetic_benzene()
    a = molsymm.analyze(zs, xyz)
    assert a.point_group == "D6h"
    assert molsymm.analyze(zs, xyz.T).point_group == "D6h"
    assert "D6" in a.subgroups and "C6v" in a.subgroups and "C2v" in a.subgroups
    for name, order in (("C2", 2), ("C3", 3), ("C6v", 12), ("D6", 12)):
        b = molsymm.analyze(zs, xyz)
        b.select_subgroup(name)
        assert b.point_group == name and b.order == order

    asmbl = a.subgroup_assembly("C2v")
    assert len(asmbl.pieces) > 1
    pieces = [xyz[atoms].copy() for atoms in asmbl.pieces]
    assert np.allclose(asmbl.assemble(pieces), xyz, atol=1e-6)
    site = 0
    for p, atoms in enumerate(asmbl.pieces):
        for i, _ in enumerate(atoms):
            if asmbl.free_axes[site]:
                pieces[p][i] = pieces[p][i] + (0.04 if p == 0 else -0.03) * np.asarray(
                    asmbl.free_axes[site][0]
                )
            site += 1
    assert not np.allclose(asmbl.assemble(pieces), xyz, atol=1e-3)


def test_synthetic_s4() -> None:
    z = [6, 6, 6, 1, 1, 1, 1]
    xyz = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.31],
            [0.0, 0.0, -1.31],
            [0.95, 0.0, 1.89],
            [-0.95, 0.0, 1.89],
            [0.0, 0.95, -1.89],
            [0.0, -0.95, -1.89],
        ]
    )
    a = molsymm.analyze(z, xyz)
    assert a.point_group == "D2d"
    assert "E" in a.irrep_names
    assert not any(n.endswith(("g", "u")) for n in a.irrep_names)
    a.select_subgroup("S4")
    assert a.point_group == "S4"
    assert a.order == 4
    assert a.symmetry_number == 2


def test_zmatrix_water() -> None:
    angle = 104.5 * np.pi / 180.0
    entries: Sequence[ZMatrixRow] = [
        (8,),
        (1, 1, 0.957),
        (1, 1, 0.957, 2, angle),
    ]
    xyz = molsymm.zmatrix_to_xyz(entries)
    assert xyz.shape == (3, 3)
    a = molsymm.analyze_zmatrix(entries)
    assert a.point_group == "C2v"
    assert a.symmetry_number == 2


def test_python_errors() -> None:
    numbers = [8, 1, 1]
    positions = np.array([[0.0, 0.0, 0.1173], [0.0, 0.7572, -0.4692], [0.0, -0.7572, -0.4692]])
    try:
        molsymm.analyze(numbers, positions[:2])
        raise AssertionError("expected Error for shape mismatch")
    except Exception:
        pass
    try:
        molsymm.analyze(numbers, positions).subgroup_assembly("C2v").assemble([])
        raise AssertionError("expected Error for piece count mismatch")
    except Exception:
        pass
    try:
        molsymm.zmatrix_to_xyz([(8, 1)])  # type: ignore[list-item]
        raise AssertionError("expected Error for bad z-matrix row")
    except Exception:
        pass


def test_symmetrize_water() -> None:
    numbers = [8, 1, 1]
    positions = np.array(
        [
            [0.0, 0.0, 0.1173],
            [0.02, 0.7572, -0.4692],
            [0.0, -0.7572, -0.4692],
        ]
    )
    forces = np.array(
        [
            [0.0, 0.0, 0.1],
            [0.05, 0.02, -0.05],
            [0.0, -0.02, -0.05],
        ]
    )
    hessian = np.eye(9)
    hessian[3, 3] += 0.3
    out = molsymm.symmetrize(
        numbers, positions, "C2v", forces=forces, hessian=hessian, primary=1e-2
    )
    pos = np.asarray(out.positions, dtype=float)
    fout = np.asarray(out.forces, dtype=float)
    hout = np.asarray(out.hessian, dtype=float)
    assert pos.shape == (3, 3)
    assert fout.shape == (3, 3)
    assert hout.shape == (9, 9)
    a = molsymm.analyze(numbers, pos, primary=1e-2)
    assert a.point_group == "C2v" and a.rmsd < 1e-8
    again = molsymm.symmetrize(numbers, pos, "C2v", forces=fout, hessian=hout)
    assert np.allclose(np.asarray(again.forces, dtype=float), fout, atol=1e-8)
    assert np.allclose(np.asarray(again.hessian, dtype=float), hout, atol=1e-8)
    assert abs(fout[1, 0]) < abs(forces[1, 0]) - 1e-4

    only_f = molsymm.symmetrize(numbers, positions, "C2v", forces=forces, primary=1e-2)
    assert only_f.forces is not None and only_f.hessian is None
    only_h = molsymm.symmetrize(numbers, positions, "C2v", hessian=hessian, primary=1e-2)
    assert only_h.forces is None and only_h.hessian is not None
    geom = molsymm.symmetrize(numbers, positions, "C2v", primary=1e-2)
    assert geom.forces is None and geom.hessian is None
    assert np.asarray(geom.positions, dtype=float).shape == (3, 3)

    flat = forces.reshape(-1)
    flat_out = molsymm.symmetrize(numbers, positions, "C2v", forces=flat, primary=1e-2)
    assert flat_out.forces is not None
    assert np.asarray(flat_out.forces, dtype=float).shape == (3, 3)

    try:
        molsymm.symmetrize(numbers, positions, "Td")
        raise AssertionError("expected RuntimeError for Td on water")
    except RuntimeError:
        pass


def test_symmetrize_benzene_subgroup() -> None:
    zs, xyz = _synthetic_benzene()
    a = molsymm.analyze(zs, xyz)
    asmbl = a.subgroup_assembly("C2v")
    pieces = [xyz[atoms].copy() for atoms in asmbl.pieces]
    site = 0
    for p, atoms in enumerate(asmbl.pieces):
        for i, _ in enumerate(atoms):
            if asmbl.free_axes[site]:
                pieces[p][i] = pieces[p][i] + 0.005 * np.asarray(asmbl.free_axes[site][0])
            site += 1
    rattled = asmbl.assemble(pieces)
    keep = np.asarray(molsymm.symmetrize(zs, rattled, "C2v", primary=1e-2).positions, dtype=float)
    full = np.asarray(molsymm.symmetrize(zs, rattled, "D6h", primary=1e-2).positions, dtype=float)
    keep_err = float(np.linalg.norm(keep - rattled))
    full_err = float(np.linalg.norm(full - rattled))
    assert keep_err < 1e-3
    assert full_err > keep_err + 1e-3


def test_symmetrize_co2() -> None:
    numbers = [6, 8, 8]
    positions = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.16], [0.0, 0.0, -1.16]])
    forces = np.zeros((3, 3))
    forces[0, 0] = 0.1
    forces[1, 2] = 0.2
    forces[2, 2] = -0.15
    out = molsymm.symmetrize(numbers, positions, "Dinfh", forces=forces)
    assert out.forces is not None
    fout = np.asarray(out.forces, dtype=float)
    assert abs(fout[0, 0]) < 1e-10
    assert abs(fout[1, 2] + fout[2, 2]) < 1e-10
    assert abs(abs(fout[1, 2]) - 0.175) < 1e-10


def _ops_permute(a: molsymm.Analysis, tol: float) -> bool:
    ops = a.operations
    if not ops:
        return False
    if a.order > 0 and len(ops) != a.order:
        return False
    numbers = list(a.atomic_numbers)
    xyz = np.asarray(a.frame_xyz, dtype=float)
    n = len(numbers)
    match = max(tol, 1e-4)
    for R in ops:
        M = np.asarray(R, dtype=float)
        claimed = [False] * n
        for i in range(n):
            rp = M @ xyz[i]
            best_j = -1
            best_d = match
            for j in range(n):
                if numbers[j] != numbers[i] or claimed[j]:
                    continue
                d = float(np.linalg.norm(rp - xyz[j]))
                if d <= best_d:
                    best_d = d
                    best_j = j
            if best_j < 0:
                return False
            claimed[best_j] = True
    return True


def test_manifest_corpus() -> None:
    rows = load_manifest()
    assert rows, "manifest must be non-empty"
    for rel, group, order, sigma, primary, max_rmsd, _source in rows:
        numbers, positions = load_xyz(DATA_DIR / rel)
        a = molsymm.analyze(numbers, positions, primary=primary)
        assert a.point_group == group, f"{rel}: expect {group} got {a.point_group}"
        assert a.order == order, f"{rel}: order {order} got {a.order}"
        assert a.symmetry_number == sigma, f"{rel}: sigma {sigma} got {a.symmetry_number}"
        assert a.rmsd < max_rmsd, f"{rel}: rmsd {a.rmsd} max {max_rmsd}"
        assert _ops_permute(a, primary), f"{rel}: operations must permute atoms"


if __name__ == "__main__":
    test_synthetic_water()
    test_synthetic_benzene_subgroups()
    test_synthetic_s4()
    test_zmatrix_water()
    test_python_errors()
    test_symmetrize_water()
    test_symmetrize_benzene_subgroup()
    test_symmetrize_co2()
    test_manifest_corpus()
    print("python corpus ok")
