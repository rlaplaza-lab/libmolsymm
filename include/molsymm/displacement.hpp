#pragma once

#include "molsymm/analysis.hpp"

#include <Eigen/Dense>
#include <span>
#include <string>
#include <vector>

namespace molsymm {

struct DisplacementImage {
  int operation_index = 0;
  int target_atom = 0;
  Mat3 rotation = Mat3::Identity();
};

struct UniqueDisplacement {
  int atom = 0;
  Vec3 direction = Vec3::UnitZ();
  std::vector<DisplacementImage> images;
  std::vector<std::string> irreps;
  std::string preserved_subgroup;
};

struct DisplacementSet {
  double step = 1e-3;
  Analysis analysis;
  std::vector<UniqueDisplacement> displacements;
  std::vector<Molecule> plus;  // original input frame
  std::vector<Molecule> minus; // original input frame

  /// Central-difference Cartesian gradient from energies on ± geometries.
  [[nodiscard]] Eigen::VectorXd gradient_from_energies(std::span<const double> e_plus,
                                                       std::span<const double> e_minus) const;

  /// Assemble Cartesian Hessian from gradients evaluated at ± geometries.
  /// Each gradient is length 3N in the original frame.
  [[nodiscard]] Eigen::MatrixXd
  hessian_from_gradients(std::span<const Eigen::VectorXd> g_plus,
                         std::span<const Eigen::VectorXd> g_minus) const;
};

DisplacementSet unique_displacements(const Analysis& analysis, double step = 1e-3);

Eigen::VectorXd symmetrize_gradient(const Analysis& analysis, const Eigen::VectorXd& gradient);
Eigen::MatrixXd symmetrize_hessian(const Analysis& analysis, const Eigen::MatrixXd& hessian);
Eigen::MatrixXd project_external(const Analysis& analysis, const Eigen::MatrixXd& hessian);

} // namespace molsymm
