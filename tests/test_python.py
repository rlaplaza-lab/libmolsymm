"""Python smoke tests for molsymm, including real XYZ geometries."""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np

import molsymm

DATA_DIR = Path(os.environ.get("MOLSYMM_TEST_DATA_DIR", Path(__file__).parent / "data"))


def _atomic_number(symbol: str) -> int:
    table = {
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
        "Pt": 78,
    }
    if symbol not in table:
        raise ValueError(f"unknown element symbol: {symbol}")
    return table[symbol]


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
    d = a.displacements()
    assert d.n == 6


def test_xyz_water() -> None:
    numbers, positions = load_xyz(DATA_DIR / "water.xyz")
    a = molsymm.analyze(numbers, positions)
    assert a.point_group == "C2v"
    assert a.symmetry_number == 2
    assert a.displacements().n == 6


def test_xyz_anion_pi() -> None:
    numbers, positions = load_xyz(DATA_DIR / "anion-pi.xyz")
    a = molsymm.analyze(numbers, positions)
    assert a.point_group == "C6v"
    assert a.symmetry_number == 6
    assert a.order == 12


def test_xyz_c60() -> None:
    numbers, positions = load_xyz(DATA_DIR / "c60.xyz")
    a = molsymm.analyze(numbers, positions)
    assert a.point_group == "Ih"
    assert a.symmetry_number == 60
    assert a.order == 120


if __name__ == "__main__":
    test_synthetic_water()
    test_xyz_water()
    test_xyz_anion_pi()
    test_xyz_c60()
    print("python smoke ok")
