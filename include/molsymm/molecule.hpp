#pragma once

#include "molsymm/types.hpp"

#include <span>
#include <vector>

namespace molsymm {

struct Molecule {
  std::vector<Atom> atoms{};

  Molecule() = default;
  Molecule(std::span<const int> atomic_numbers, std::span<const double> xyz_flat);
  Molecule(std::span<const int> atomic_numbers,
           const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>& xyz);

  [[nodiscard]] int size() const { return static_cast<int>(atoms.size()); }
  void resolve_masses();
};

Vec3 center_of_mass(const Molecule& mol);
Mat3 inertia_tensor(const Molecule& mol, const Vec3& origin);

} // namespace molsymm
