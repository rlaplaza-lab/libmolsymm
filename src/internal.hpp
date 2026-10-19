#pragma once

#include "molsymm/molecule.hpp"
#include "molsymm/types.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <string>
#include <vector>

namespace molsymm::detail {

inline constexpr double kPi = 3.14159265358979323846;

inline int parse_n(const std::string& name) {
  std::size_t i = 0;
  while (i < name.size() && !std::isdigit(static_cast<unsigned char>(name[i]))) {
    ++i;
  }
  if (i >= name.size()) {
    return 0;
  }
  int n = 0;
  while (i < name.size() && std::isdigit(static_cast<unsigned char>(name[i]))) {
    n = n * 10 + (name[i] - '0');
    ++i;
  }
  return n;
}

inline bool parallel_axis(const Vec3& a, const Vec3& b, double angle_tol) {
  const double na = a.norm();
  const double nb = b.norm();
  if (na < 1e-12 || nb < 1e-12) {
    return false;
  }
  const double c = std::abs(a.dot(b) / (na * nb));
  return c >= std::cos(angle_tol);
}

inline Mat3 rotation_matrix(const Vec3& axis, double angle) {
  Vec3 u = axis.normalized();
  const double c = std::cos(angle);
  const double s = std::sin(angle);
  Mat3 K = Mat3::Zero();
  K(0, 1) = -u.z();
  K(0, 2) = u.y();
  K(1, 0) = u.z();
  K(1, 2) = -u.x();
  K(2, 0) = -u.y();
  K(2, 1) = u.x();
  return Mat3::Identity() + s * K + (1.0 - c) * (K * K);
}

inline Mat3 reflection_matrix(const Vec3& normal) {
  Vec3 n = normal.normalized();
  return Mat3::Identity() - 2.0 * n * n.transpose();
}

inline Mat3 improper_rotation_matrix(const Vec3& axis, double angle) {
  return reflection_matrix(axis) * rotation_matrix(axis, angle);
}

inline Mat3 inversion_matrix() {
  return -Mat3::Identity();
}

inline bool matrices_equal(const Mat3& a, const Mat3& b, double tol) {
  return (a - b).cwiseAbs().maxCoeff() <= tol;
}

inline void dedupe_operations(std::vector<Operation>& ops, double tol) {
  std::vector<Operation> uniq;
  uniq.reserve(ops.size());
  for (const Operation& op : ops) {
    bool found = false;
    for (const Operation& u : uniq) {
      if (matrices_equal(u.matrix, op.matrix, tol)) {
        found = true;
        break;
      }
    }
    if (!found) {
      uniq.push_back(op);
    }
  }
  ops.swap(uniq);
}

using ElementIndex = std::vector<std::vector<int>>;

inline ElementIndex build_element_index(const std::vector<Atom>& atoms) {
  int max_z = 0;
  for (const Atom& a : atoms) {
    max_z = std::max(max_z, a.atomic_number);
  }
  ElementIndex idx(static_cast<std::size_t>(std::max(max_z, 0)) + 1);
  for (int i = 0; i < static_cast<int>(atoms.size()); ++i) {
    const int z = atoms[static_cast<std::size_t>(i)].atomic_number;
    if (z >= 0) {
      idx[static_cast<std::size_t>(z)].push_back(i);
    }
  }
  return idx;
}

int find_image(
    const std::vector<Atom>& atoms, const ElementIndex& by_z, int from, const Mat3& R, double tol);

std::vector<int>
permutation(const std::vector<Atom>& atoms, const ElementIndex& by_z, const Mat3& R, double tol);

CharacterTable
build_character_table(const std::string& name, int n, const std::vector<Operation>& ops);

int symmetry_number_for(const std::string& name, int n);

std::vector<std::string> constructive_subgroups(const std::string& name, int n);

std::string classify_point_group(const std::vector<Operation>& ops,
                                 bool linear,
                                 bool has_inversion,
                                 bool has_sigma_h,
                                 bool has_sigma_v,
                                 int highest_cn,
                                 int n_perp_c2,
                                 int n_c3,
                                 int n_c4,
                                 int n_c5,
                                 bool spherical_top);

void close_group(std::vector<Operation>& ops,
                 const std::vector<Atom>& atoms,
                 double tol,
                 int max_order);

void assign_conjugacy_classes(std::vector<Operation>& ops, double tol);

void align_to_frame(Molecule& mol,
                    std::vector<Operation>& ops,
                    Mat3& orientation,
                    const Tolerances& tol);

} // namespace molsymm::detail

namespace molsymm {

std::vector<EquivalenceSet> build_equivalence_sets(const Molecule& mol, double tol);

} // namespace molsymm
