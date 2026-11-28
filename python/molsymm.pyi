"""Molecular point-group symmetry for Cartesian and Z-matrix geometries.

Import as ``molsymm`` (PyPI package ``libmolsymm``). Coordinates are ångströms;
Z-matrix angles are radians with 1-based indices (``0`` unused). Fragmentation
indices are 0-based. Prefer positions shaped ``(n, 3)`` (``(3, n)`` is also
accepted). Numeric results are C-contiguous ``float64`` ndarrays. Failures
raise ``RuntimeError``.

Typical workflow:

1. ``a = analyze(numbers, positions)`` (or ``analyze_zmatrix``)
2. Read ``a.point_group``, ``a.order``, ``a.operations``, …
3. Optionally ``a.select_subgroup(...)``, then ``a.symmetrize()``,
   ``a.displacements()``, ``a.fragmentation()``, or ``a.subgroup_assembly(...)``
4. Or one-shot ``symmetrize(numbers, positions, point_group, ...)`` without a
   prior analysis object
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
"""One Z-matrix row (1-based bond/angle/dihedral indices; Å / radians).

Notes
-----
Accepted lengths are 1, 3, 5, or 7:

- ``(Z,)`` — first atom
- ``(Z, bond_to, R)`` — second atom
- ``(Z, bond_to, R, angle_to, angle)`` — third atom
- ``(Z, bond_to, R, angle_to, angle, dihedral_to, dihedral)`` — later atoms
"""

class Analysis:
    """Point-group analysis result for one molecule.

    Created by ``analyze`` or ``analyze_zmatrix``. Stores the detected (or
    selected) point group and the molecule in the analysis frame (COM at the
    origin, principal axes). Use the methods below to project the geometry,
    build symmetry-reduced finite differences, or recover an asymmetric unit.

    Attributes
    ----------
    point_group : str
        Detected or selected point-group label (e.g. ``C2v``).
    order : int
        Group order; ``0`` for infinite groups (``Cinfv``, ``Dinfh``).
    symmetry_number : int
        Rotational symmetry number for thermochemistry.
    rmsd : float
        RMS change from the last ``symmetrize`` call, else ``0.0``.
    subgroups : list of str
        Labels accepted by ``select_subgroup`` and ``subgroup_assembly``.

    Examples
    --------
    >>> a = molsymm.analyze(numbers, positions)
    >>> a.point_group, a.order
    ('C2v', 4)
    >>> rms, xyz = a.symmetrize()
    """

    point_group: str
    """Detected (or selected) point-group label."""
    order: int
    """Group order; ``0`` for infinite groups."""
    symmetry_number: int
    """Rotational symmetry number for thermochemistry."""
    rmsd: float
    """RMS change after the last ``symmetrize`` call, else ``0.0``."""
    subgroups: list[str]
    """Named subgroups for ``select_subgroup`` / ``subgroup_assembly``."""

    @property
    def atomic_numbers(self) -> list[int]:
        """Atomic numbers in analysis-frame atom order."""

    @property
    def frame_xyz(self) -> Float64Array:
        """Analysis-frame coordinates, shape ``(n, 3)`` (COM at origin)."""

    @property
    def operations(self) -> list[Float64Array]:
        """Symmetry operation matrices in the analysis frame, each ``(3, 3)``."""

    @property
    def irrep_names(self) -> list[str]:
        """Irreducible-representation labels for the current group."""

    def select_subgroup(self, name: str) -> None:
        """Replace the current group with a named subgroup (in place).

        After this call, ``point_group``, ``order``, ``operations``, and related
        fields describe ``name``. Subsequent ``symmetrize`` / ``displacements`` /
        ``fragmentation`` use the subgroup.

        Parameters
        ----------
        name : str
            A label from ``subgroups``.

        Raises
        ------
        RuntimeError
            If ``name`` is not in ``subgroups``.
        """

    def symmetrize(self) -> tuple[float, Float64Array]:
        """Project the geometry onto the totally symmetric subspace.

        Updates ``rmsd``. The returned coordinates are in the **input** frame
        (the frame passed to ``analyze``), not ``frame_xyz``.

        Returns
        -------
        rms : float
            RMS Cartesian change in the analysis frame (also stored in ``rmsd``).
        xyz : ndarray of float64, shape (n, 3)
            Symmetrized positions in the input frame.

        Examples
        --------
        >>> a = molsymm.analyze(numbers, positions)
        >>> rms, xyz = a.symmetrize()
        """

    def displacements(self, step: float = 1e-3) -> DisplacementSet:
        """Build symmetry-reduced ± Cartesian displacements.

        Returns unique displacements in the **input** frame for finite-difference
        gradients and Hessians.

        Parameters
        ----------
        step : float, optional
            Cartesian step in Å (default ``1e-3``).

        Returns
        -------
        DisplacementSet
            ``plus`` / ``minus`` geometries and FD assembly helpers.
        """

    def fragmentation(self) -> Fragmentation:
        """Partition the molecule into an asymmetric unit (finite groups).

        Returns
        -------
        Fragmentation
            Representative atom indices and ``reconstruct`` to rebuild the full
            geometry from those sites.

        Raises
        ------
        RuntimeError
            If the current group is infinite.
        """

    def subgroup_assembly(self, name: str) -> SubgroupAssembly:
        """Describe how to rebuild the molecule under a lower-symmetry subgroup.

        Parameters
        ----------
        name : str
            A label from ``subgroups``.

        Returns
        -------
        SubgroupAssembly
            Atom-index pieces and ``assemble`` to rebuild the full geometry.

        Raises
        ------
        RuntimeError
            If ``name`` is not in ``subgroups``.
        """

class DisplacementSet:
    """Symmetry-reduced ± Cartesian displacements from ``Analysis.displacements``.

    Evaluate energies or gradients at ``plus`` / ``minus``, then call
    ``gradient_from_energies`` or ``hessian_from_gradients``.

    Attributes
    ----------
    step : float
        Cartesian step size in Å.
    """

    step: float
    """Cartesian step size in Å used to build ``plus`` / ``minus``."""

    @property
    def n(self) -> int:
        """Number of unique displacements (``len(plus)`` / ``len(minus)``)."""

    @property
    def plus(self) -> list[Float64Array]:
        """``+step`` geometries in the input frame, each shape ``(n_atoms, 3)``."""

    @property
    def minus(self) -> list[Float64Array]:
        """``-step`` geometries in the input frame, each shape ``(n_atoms, 3)``."""

    def gradient_from_energies(
        self,
        ep: ArrayLike,
        em: ArrayLike,
    ) -> Float64Array:
        """Central-difference gradient from ± energies.

        Parameters
        ----------
        ep : array_like
            Energies at ``plus`` geometries, length ``n``.
        em : array_like
            Energies at ``minus`` geometries, length ``n``.

        Returns
        -------
        ndarray of float64, shape (3 * n_atoms,)
            Gradient in the input frame.
        """

    def hessian_from_gradients(
        self,
        gp: Sequence[ArrayLike],
        gm: Sequence[ArrayLike],
    ) -> Float64Array:
        """Assemble a Hessian from ± gradients.

        Parameters
        ----------
        gp : sequence of array_like
            Gradients at ``plus`` geometries (length ``n``, each length
            ``3 * n_atoms``).
        gm : sequence of array_like
            Gradients at ``minus`` geometries (same shapes as ``gp``).

        Returns
        -------
        ndarray of float64, shape (3 * n_atoms, 3 * n_atoms)
            Hessian in the input frame.
        """

class Fragmentation:
    """Asymmetric unit of a finite point group from ``Analysis.fragmentation``.

    ``representatives`` are the independent atoms. Pass their coordinates to
    ``reconstruct`` to generate the full molecule by symmetry.

    Attributes
    ----------
    representatives : list of int
        One 0-based atom index per orbit (lowest index).
    """

    representatives: list[int]
    """One 0-based atom index per orbit (lowest index)."""

    @property
    def copies(self) -> list[list[int]]:
        """Disjoint 0-based atom-index partitions of the molecule."""

    @property
    def mask(self) -> list[int]:
        """Length-``n`` mask with ``1`` on representatives and ``0`` elsewhere."""

    @property
    def free_axes(self) -> list[list[list[float]]]:
        """Input-frame free axes per representative; empty site lists are falsy."""

    def reconstruct(self, positions: ArrayLike) -> Float64Array:
        """Rebuild the full geometry from representative coordinates.

        Parameters
        ----------
        positions : array_like
            Representative coordinates, shape ``(n_rep, 3)`` or ``(3, n_rep)``,
            in the same order as ``representatives``.

        Returns
        -------
        ndarray of float64, shape (n_atoms, 3)
            Full geometry in the input frame.

        Raises
        ------
        RuntimeError
            If the shape does not match ``representatives``.
        """

class SubgroupAssembly:
    """Piece lists for assembling under a subgroup from ``subgroup_assembly``.

    Attributes
    ----------
    subgroup : str
        Subgroup label used to build this assembly.
    """

    subgroup: str
    """Subgroup label used to build this assembly."""

    @property
    def pieces(self) -> list[list[int]]:
        """0-based atom-index lists, one per independent piece."""

    @property
    def free_axes(self) -> list[list[list[float]]]:
        """Input-frame free axes in piece order; empty site lists are falsy."""

    def assemble(self, pieces: Sequence[ArrayLike]) -> Float64Array:
        """Assemble the full geometry from one coordinate array per piece.

        Parameters
        ----------
        pieces : sequence of array_like
            One array per entry in ``pieces``, each shaped ``(n_i, 3)`` or
            ``(3, n_i)``.

        Returns
        -------
        ndarray of float64, shape (n_atoms, 3)
            Full geometry in the input frame.

        Raises
        ------
        RuntimeError
            If the number or shapes of ``pieces`` do not match.
        """

class Symmetrized:
    """Result of the module-level ``symmetrize`` function.

    Holds the projected geometry in the **input** frame and optional projected
    forces / Hessian when those were passed in.

    Attributes
    ----------
    rmsd : float
        RMS change of the geometry projection in the analysis frame.
    """

    rmsd: float
    """RMS change of the geometry projection in the analysis frame."""

    @property
    def positions(self) -> Float64Array:
        """Symmetrized coordinates, shape ``(n, 3)``, in the input frame."""

    @property
    def forces(self) -> Float64Array | None:
        """Symmetrized forces ``(n, 3)``, or ``None`` if omitted (also if input was flat)."""

    @property
    def hessian(self) -> Float64Array | None:
        """Symmetrized Hessian ``(3n, 3n)``, or ``None`` if omitted."""

def analyze(
    numbers: Sequence[int],
    positions: ArrayLike,
    primary: float = 1e-3,
) -> Analysis:
    """Detect the point group of a Cartesian geometry.

    This is the main entry point. The returned ``Analysis`` object carries the
    group label and methods such as ``symmetrize``, ``displacements``, and
    ``fragmentation``.

    Parameters
    ----------
    numbers : sequence of int
        Atomic numbers, length ``n``.
    positions : array_like
        Cartesian coordinates, shape ``(n, 3)`` or ``(3, n)``, in ångströms.
    primary : float, optional
        Geometry tolerance in Å (default ``1e-3``).

    Returns
    -------
    Analysis
        Point-group analysis for this molecule.

    Raises
    ------
    RuntimeError
        If ``numbers`` and ``positions`` are incompatible or analysis fails.

    Examples
    --------
    >>> a = molsymm.analyze([8, 1, 1], positions)
    >>> a.point_group
    'C2v'
    """

def symmetrize(
    numbers: Sequence[int],
    positions: ArrayLike,
    point_group: str,
    forces: ArrayLike | None = None,
    hessian: ArrayLike | None = None,
    primary: float = 1e-3,
) -> Symmetrized:
    """Project geometry (and optional forces / Hessian) onto a named point group.

    One-shot alternative to ``analyze(...).symmetrize()`` when you already know
    the target group and may also want forces / Hessian projected. Returns a
    ``Symmetrized`` object; coordinates are in the **input** frame.

    Parameters
    ----------
    numbers : sequence of int
        Atomic numbers, length ``n``.
    positions : array_like
        Cartesian coordinates, shape ``(n, 3)`` or ``(3, n)``.
    point_group : str
        Target point-group label (e.g. ``C2v``, ``D6h``).
    forces : array_like, optional
        Forces shaped ``(n, 3)``, ``(3, n)``, or flat ``(3n,)``. If omitted,
        ``result.forces`` is ``None``.
    hessian : array_like, optional
        Hessian shaped ``(3n, 3n)``. If omitted, ``result.hessian`` is ``None``.
    primary : float, optional
        Geometry tolerance in Å (default ``1e-3``).

    Returns
    -------
    Symmetrized
        ``positions`` and optional ``forces`` / ``hessian`` after projection.

    Raises
    ------
    RuntimeError
        If shapes are incompatible or ``point_group`` cannot be applied.

    Examples
    --------
    >>> out = molsymm.symmetrize(numbers, positions, "C2v", forces=forces)
    >>> out.positions.shape
    (n, 3)
    """

def analyze_zmatrix(
    entries: Sequence[ZMatrixRow],
    primary: float = 1e-3,
) -> Analysis:
    """Build Cartesian coordinates from a Z-matrix and detect the point group.

    Same return type and methods as ``analyze``.

    Parameters
    ----------
    entries : sequence of ZMatrixRow
        Rows of length 1, 3, 5, or 7 (see ``ZMatrixRow``). Indices are 1-based;
        distances in Å, angles in radians.
    primary : float, optional
        Geometry tolerance in Å (default ``1e-3``).

    Returns
    -------
    Analysis
        Point-group analysis of the built geometry.

    Raises
    ------
    RuntimeError
        If a row is malformed or analysis fails.

    Examples
    --------
    >>> entries = [(8,), (1, 1, 0.957), (1, 1, 0.957, 2, angle)]
    >>> molsymm.analyze_zmatrix(entries).point_group
    'C2v'
    """

def zmatrix_to_xyz(entries: Sequence[ZMatrixRow]) -> Float64Array:
    """Convert a Z-matrix to Cartesian coordinates without symmetry analysis.

    Parameters
    ----------
    entries : sequence of ZMatrixRow
        Same row forms as ``analyze_zmatrix``.

    Returns
    -------
    ndarray of float64, shape (n, 3)
        Cartesian positions in ångströms.

    Raises
    ------
    RuntimeError
        If a row is malformed.
    """

class WyckoffSite:
    """One special position of a finite point group (Cotton frame)."""

    label: str
    site_symmetry: str
    multiplicity: int
    degrees_of_freedom: int
    projector: Float64Array
    chamber: Float64Array

class WyckoffOrbit:
    """One orbit assignment inside a ``WyckoffPlan``."""

    atomic_number: int
    label: str
    multiplicity: int
    degrees_of_freedom: int

class WyckoffPlan:
    """Orbit assignment for a composition under a point group."""

    point_group: str
    order: int
    orbits: list[WyckoffOrbit]

    def assemble(self, positions: ArrayLike) -> Analysis:
        """Expand one Cotton-frame position per orbit into an ``Analysis``."""

def wyckoff_positions(point_group: str) -> list[WyckoffSite]:
    """List special positions of a finite Schönflies group."""

def wyckoff_plans(
    point_group: str,
    composition: Sequence[int] | int,
) -> list[WyckoffPlan]:
    """Partition a composition (or ``n_sites``) into orbit plans."""

def wyckoff_plan(
    point_group: str,
    site_labels: Sequence[str],
    atomic_numbers: Sequence[int] = (),
) -> WyckoffPlan:
    """Build one plan from explicit site labels."""

def assemble_from_au(
    point_group: str,
    numbers: Sequence[int],
    positions: ArrayLike,
    site_labels: Sequence[str],
) -> Analysis:
    """Expand asymmetric-unit coordinates for the given sites into an ``Analysis``."""
