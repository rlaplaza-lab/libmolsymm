#pragma once

#include "molsymm/molecule.hpp"

#include <span>
#include <vector>

namespace molsymm {

/// One row of a Z-matrix in the Priroda / briling placement.
/// Indices are 1-based; 0 means the field is unused. Distances are Å and angles are radians.
struct ZMatrixAtom {
  int atomic_number = 0;
  int bond_to = 0;
  int angle_to = 0;
  int dihedral_to = 0;
  double distance = 0;
  double angle = 0;
  double dihedral = 0;
};

/// Place atoms with the same rules as briling/v `zmat2cart`.
[[nodiscard]] std::vector<Vec3> zmatrix_positions(std::span<const ZMatrixAtom> atoms);

/// Build a Cartesian molecule from a Z-matrix. Masses are resolved from atomic numbers.
[[nodiscard]] Molecule molecule_from_zmatrix(std::span<const ZMatrixAtom> atoms);

} // namespace molsymm
