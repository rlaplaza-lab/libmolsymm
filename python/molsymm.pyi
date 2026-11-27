"""Molecular point-group symmetry, including Z-matrix input.

Coordinates are ångströms; Z-matrix angles are radians (1-based indices,
``0`` unused). Fragmentation indices are 0-based. Prefer ``(n, 3)`` positions
(``(3, n)`` also accepted). Failures raise :class:`RuntimeError`.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import TypeAlias

import numpy as np
import numpy.typing as npt

ArrayLike: TypeAlias = npt.ArrayLike
Float64Array: TypeAlias = npt.NDArray[np.float64]

ZMatrixRow: TypeAlias = (
    tuple[int]
    | tuple[int, int, float]
    | tuple[int, int, float, int, float]
    | tuple[int, int, float, int, float, int, float]
)

class Analysis:
    """Point-group analysis of a molecule."""

    point_group: str
    order: int
    """Group order; ``0`` for infinite groups."""
    symmetry_number: int
    rmsd: float
    """RMS change after the last :meth:`symmetrize`, else ``0.0``."""
    subgroups: list[str]

    @property
    def atomic_numbers(self) -> list[int]:
        """Atomic numbers in analysis-frame atom order."""

    @property
    def frame_xyz(self) -> Float64Array:
        """Analysis-frame coordinates ``(n, 3)`` (COM at origin)."""

    @property
    def operations(self) -> list[Float64Array]:
        """Symmetry operation matrices in the analysis frame, each ``(3, 3)``."""

    @property
    def irrep_names(self) -> list[str]:
        """Irreducible-representation labels."""

    def select_subgroup(self, name: str) -> None:
        """Restrict to a named subgroup in :attr:`subgroups` (in place)."""

    def symmetrize(self) -> tuple[float, Float64Array]:
        """Project onto the totally symmetric geometry; updates :attr:`rmsd`.

        Returns
        -------
        rms : float
            RMS change in the analysis frame.
        xyz : ndarray of float64, shape (n, 3)
            Symmetrized positions in the input frame.
        """

    def displacements(self, step: float = 1e-3) -> DisplacementSet:
        """Symmetry-reduced finite-difference displacements (input frame).

        Parameters
        ----------
        step : float, optional
            Cartesian step in Å (default ``1e-3``).
        """

    def fragmentation(self) -> Fragmentation:
        """Asymmetric-unit partition (finite groups only)."""

    def subgroup_assembly(self, name: str) -> SubgroupAssembly:
        """Pieces for rebuilding under a named subgroup in :attr:`subgroups`."""

class DisplacementSet:
    """Symmetry-reduced ± Cartesian displacements."""

    step: float

    @property
    def n(self) -> int:
        """Number of unique displacements."""

    @property
    def plus(self) -> list[Float64Array]:
        """``+step`` geometries in the input frame, each ``(n_atoms, 3)``."""

    @property
    def minus(self) -> list[Float64Array]:
        """``-step`` geometries in the input frame, each ``(n_atoms, 3)``."""

    def gradient_from_energies(
        self,
        ep: ArrayLike,
        em: ArrayLike,
    ) -> Float64Array:
        """Central-difference gradient from ± energies (length :attr:`n` each).

        Returns
        -------
        ndarray of float64, shape (3 * n_atoms,)
            Input frame.
        """

    def hessian_from_gradients(
        self,
        gp: Sequence[ArrayLike],
        gm: Sequence[ArrayLike],
    ) -> Float64Array:
        """Hessian from ± gradients (length :attr:`n`, each vector ``3 * n_atoms``).

        Returns
        -------
        ndarray of float64, shape (3 * n_atoms, 3 * n_atoms)
            Input frame.
        """

class Fragmentation:
    """Asymmetric unit of a finite point group."""

    representatives: list[int]
    """One 0-based atom index per orbit (lowest index)."""

    @property
    def copies(self) -> list[list[int]]:
        """Disjoint partitions of the molecule (0-based)."""

    @property
    def mask(self) -> list[int]:
        """Length-``n`` mask; ``1`` on representatives."""

    @property
    def free_axes(self) -> list[list[list[float]]]:
        """Input-frame free axes per representative; empty site lists are falsy."""

    def reconstruct(self, positions: ArrayLike) -> Float64Array:
        """Rebuild full geometry from representative coords ``(n_rep, 3)`` / ``(3, n_rep)``.

        Returns
        -------
        ndarray of float64, shape (n_atoms, 3)
            Input frame.
        """

class SubgroupAssembly:
    """Pieces for assembling under a lower-symmetry subgroup."""

    subgroup: str

    @property
    def pieces(self) -> list[list[int]]:
        """0-based atom-index lists, one per piece."""

    @property
    def free_axes(self) -> list[list[list[float]]]:
        """Input-frame free axes in piece order; empty site lists are falsy."""

    def assemble(self, pieces: Sequence[ArrayLike]) -> Float64Array:
        """Assemble full geometry from one array per :attr:`pieces` entry.

        Returns
        -------
        ndarray of float64, shape (n_atoms, 3)
            Input frame.
        """

class Symmetrized:
    """Result of one-shot projection onto a named point group."""

    rmsd: float

    @property
    def positions(self) -> Float64Array:
        """Symmetrized coordinates ``(n, 3)`` in the input frame."""

    @property
    def forces(self) -> Float64Array | None:
        """Symmetrized forces ``(n, 3)``, or ``None`` if omitted (also when input was flat)."""

    @property
    def hessian(self) -> Float64Array | None:
        """Symmetrized Hessian ``(3n, 3n)``, or ``None`` if omitted."""

def analyze(
    numbers: Sequence[int],
    positions: ArrayLike,
    primary: float = 1e-3,
) -> Analysis:
    """Detect the point group of a Cartesian geometry.

    Parameters
    ----------
    numbers : sequence of int
        Atomic numbers, length ``n``.
    positions : array_like
        Shape ``(n, 3)`` or ``(3, n)``.
    primary : float, optional
        Geometry tolerance in Å (default ``1e-3``).
    """

def symmetrize(
    numbers: Sequence[int],
    positions: ArrayLike,
    point_group: str,
    forces: ArrayLike | None = None,
    hessian: ArrayLike | None = None,
    primary: float = 1e-3,
) -> Symmetrized:
    """Project geometry and optional forces / Hessian onto a named group.

    Parameters
    ----------
    numbers : sequence of int
        Atomic numbers, length ``n``.
    positions : array_like
        Shape ``(n, 3)`` or ``(3, n)``.
    point_group : str
        Target point-group label.
    forces : array_like, optional
        ``(n, 3)``, ``(3, n)``, or flat ``(3n,)``. Omitted → ``forces`` is ``None``.
    hessian : array_like, optional
        ``(3n, 3n)``. Omitted → ``hessian`` is ``None``.
    primary : float, optional
        Geometry tolerance in Å (default ``1e-3``).
    """

def analyze_zmatrix(
    entries: Sequence[ZMatrixRow],
    primary: float = 1e-3,
) -> Analysis:
    """Build from a Z-matrix and analyze.

    Parameters
    ----------
    entries : sequence of Z-matrix rows
        Length 1, 3, 5, or 7: ``(Z,)``, ``(Z, bond_to, R)``,
        ``(Z, bond_to, R, angle_to, angle)``,
        ``(Z, bond_to, R, angle_to, angle, dihedral_to, dihedral)``.
        1-based indices; Å / radians.
    primary : float, optional
        Geometry tolerance in Å (default ``1e-3``).
    """

def zmatrix_to_xyz(entries: Sequence[ZMatrixRow]) -> Float64Array:
    """Z-matrix → Cartesian ``(n, 3)`` (same row forms as :func:`analyze_zmatrix`)."""
