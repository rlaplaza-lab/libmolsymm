#pragma once

#include "molsymm/molecule.hpp"
#include "molsymm/types.hpp"

#include <string>
#include <vector>

namespace molsymm {

struct Analysis {
  std::string point_group;
  int order = 1; // 0 for infinite groups
  int symmetry_number = 1;
  bool infinite = false;
  bool linear = false;

  Molecule molecule; // in Cotton/principal-axis frame, COM at origin
  Molecule original; // input molecule (unmodified)
  Vec3 center_of_mass = Vec3::Zero();
  Mat3 principal_axes = Mat3::Identity(); // columns are axes; R maps body→lab as R * r_body
  Mat3 orientation = Mat3::Identity();    // maps input (COM-centered) → frame used in molecule
  std::vector<EquivalenceSet> equivalence_sets;
  std::vector<Operation> operations;
  CharacterTable character_table;
  std::vector<std::string> subgroups;
  double rmsd = 0.0; // after symmetrization of the analysis geometry

  /// Map a point from the analysis frame back to the original input coordinates.
  [[nodiscard]] Vec3 frame_to_input(const Vec3& p) const {
    return orientation.transpose() * p + center_of_mass;
  }

  /// Restrict operations and rebuild the character table for a named subgroup.
  void select_subgroup(const std::string& name);

  /// Project coordinates onto the totally symmetric geometry. Returns RMS change.
  double symmetrize();
};

Analysis analyze(const Molecule& molecule, const Tolerances& tol = {});

} // namespace molsymm
