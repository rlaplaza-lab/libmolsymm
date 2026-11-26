"""Python tests for molsymm against the IQmol manifest corpus."""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np

import molsymm

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
    entries = [
        (8,),
        (1, 1, 0.957),
        (1, 1, 0.957, 2, angle),
    ]
    xyz = molsymm.zmatrix_to_xyz(entries)
    assert xyz.shape == (3, 3)
    a = molsymm.analyze_zmatrix(entries)
    assert a.point_group == "C2v"
    assert a.symmetry_number == 2


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
    test_manifest_corpus()
    print("python corpus ok")
