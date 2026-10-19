#pragma once

#include <Eigen/Dense>
#include <cstdint>
#include <string>
#include <vector>

namespace molsymm {

using Vec3 = Eigen::Vector3d;
using Mat3 = Eigen::Matrix3d;

enum class OperationKind : std::uint8_t {
  Identity,
  Inversion,
  ProperRotation,
  ImproperRotation,
  Reflection
};

enum class PlaneKind : std::uint8_t { None, Horizontal, Vertical, Dihedral };

struct Atom {
  int atomic_number = 0;
  Vec3 position = Vec3::Zero();
  double mass = 0.0; // 0 means use standard atomic weight
};

struct Tolerances {
  double primary = 1e-3; // Å — geometry / equivalence matching
  double angle = 1e-2;   // radians for axis/plane comparisons
  double inertia = 1e-2; // relative inertia degeneracy
  int max_axis_order = 24;
  int max_group_order = 240;
};

struct Operation {
  OperationKind kind = OperationKind::Identity;
  int order = 1; // n for Cn / Sn; 1 for E, i, σ
  int power = 0; // k in C_n^k or S_n^k
  Vec3 axis = Vec3::UnitZ();
  PlaneKind plane = PlaneKind::None;
  Mat3 matrix = Mat3::Identity();
  std::vector<int> image; // atom permutation: image[i] = j means g maps atom i → j
  int class_index = 0;
};

struct CharacterTable {
  std::vector<std::string> irrep_names;
  std::vector<int> irrep_dims;
  std::vector<std::string> class_labels;
  std::vector<int> class_sizes;
  /// characters[irrep][class]
  std::vector<std::vector<double>> characters;
};

struct EquivalenceSet {
  std::vector<int> atoms;
};

} // namespace molsymm
