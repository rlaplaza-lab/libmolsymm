#include "internal.hpp"

#include "molsymm/error.hpp"

#include <algorithm>
#include <cmath>
#include <utility>

namespace molsymm::detail {

int find_image(
    const std::vector<Atom>& atoms, const ElementIndex& by_z, int from, const Mat3& R, double tol) {
  const Vec3 target = R * atoms[static_cast<std::size_t>(from)].position;
  const int z = atoms[static_cast<std::size_t>(from)].atomic_number;
  if (z < 0 || static_cast<std::size_t>(z) >= by_z.size()) {
    return -1;
  }
  int best = -1;
  double best_d = tol + 1.0;
  for (int j : by_z[static_cast<std::size_t>(z)]) {
    const double d = (atoms[static_cast<std::size_t>(j)].position - target).norm();
    if (d < best_d) {
      best_d = d;
      best = j;
    }
  }
  return best_d <= tol ? best : -1;
}

std::vector<int>
permutation(const std::vector<Atom>& atoms, const ElementIndex& by_z, const Mat3& R, double tol) {
  std::vector<int> image(atoms.size(), -1);
  std::vector<char> taken(atoms.size(), 0);
  for (int i = 0; i < static_cast<int>(atoms.size()); ++i) {
    const int j = find_image(atoms, by_z, i, R, tol);
    if (j < 0 || taken[static_cast<std::size_t>(j)]) {
      return {};
    }
    taken[static_cast<std::size_t>(j)] = 1;
    image[static_cast<std::size_t>(i)] = j;
  }
  return image;
}

bool classify_proper_rotation(const Mat3& M, int& order, int& power, Vec3& axis) {
  const double tr = M.trace();
  const double c = std::clamp((tr - 1.0) * 0.5, -1.0, 1.0);
  Vec3 skew(M(2, 1) - M(1, 2), M(0, 2) - M(2, 0), M(1, 0) - M(0, 1));
  const double s = 0.5 * skew.norm();
  const double ang = std::atan2(s, c); // [0, π]
  if (ang < 1e-8) {
    order = 1;
    power = 0;
    axis = Vec3::UnitZ();
    return true;
  }
  if (skew.norm() > 1e-8) {
    axis = skew.normalized();
  } else {
    // 180°: axis is the +1 eigenvector.
    Eigen::SelfAdjointEigenSolver<Mat3> es(0.5 * (M + M.transpose()));
    axis = Vec3::UnitZ();
    for (int k = 0; k < 3; ++k) {
      if (std::abs(es.eigenvalues()(k) - 1.0) < 1e-4) {
        axis = es.eigenvectors().col(k);
        break;
      }
    }
  }

  constexpr double kAngTol = 1e-3;
  int best_n = 0;
  int best_k = 1;
  double best_err = 1e9;
  for (int n = 2; n <= 48; ++n) {
    for (int k = 1; k <= n / 2; ++k) {
      const double err = std::abs(ang - 2.0 * kPi * k / n);
      if (err < best_err - 1e-15 ||
          (std::abs(err - best_err) <= 1e-15 && (best_n == 0 || n < best_n))) {
        best_err = err;
        best_n = n;
        best_k = k;
      }
    }
  }
  if (best_n == 0 || best_err > kAngTol) {
    return false;
  }
  order = best_n;
  power = best_k;
  return true;
}

void close_group(std::vector<Operation>& ops,
                 const std::vector<Atom>& atoms,
                 double match_tol,
                 double matrix_tol,
                 int max_order) {
  const ElementIndex by_z = build_element_index(atoms);

  bool has_e = false;
  for (const Operation& op : ops) {
    if (op.kind == OperationKind::Identity) {
      has_e = true;
      break;
    }
  }
  if (!has_e) {
    Operation e;
    e.kind = OperationKind::Identity;
    e.matrix = Mat3::Identity();
    e.image.resize(atoms.size());
    for (int i = 0; i < static_cast<int>(atoms.size()); ++i) {
      e.image[static_cast<std::size_t>(i)] = i;
    }
    ops.push_back(e);
  }

  for (Operation& op : ops) {
    if (op.image.empty()) {
      op.image = permutation(atoms, by_z, op.matrix, match_tol);
    }
  }

  bool grew = true;
  while (grew) {
    grew = false;
    const int n = static_cast<int>(ops.size());
    if (n > max_order) {
      throw Error("group order exceeded max_group_order");
    }
    for (int i = 0; i < n && !grew; ++i) {
      for (int j = 0; j < n; ++j) {
        Mat3 M = ops[static_cast<std::size_t>(i)].matrix * ops[static_cast<std::size_t>(j)].matrix;
        bool found = false;
        for (const Operation& op : ops) {
          if (matrices_equal(op.matrix, M, matrix_tol)) {
            found = true;
            break;
          }
        }
        if (found) {
          continue;
        }
        std::vector<int> image = permutation(atoms, by_z, M, match_tol);
        if (image.empty()) {
          continue;
        }
        Operation neu;
        neu.matrix = M;
        neu.image = std::move(image);
        const double det = M.determinant();
        const double tr = M.trace();
        if (matrices_equal(M, Mat3::Identity(), matrix_tol)) {
          neu.kind = OperationKind::Identity;
          neu.order = 1;
        } else if (matrices_equal(M, -Mat3::Identity(), matrix_tol)) {
          neu.kind = OperationKind::Inversion;
          neu.order = 1;
        } else if (std::abs(det - 1.0) < 1e-6) {
          neu.kind = OperationKind::ProperRotation;
          int ord = 0;
          int pow = 0;
          Vec3 ax = Vec3::UnitZ();
          if (!classify_proper_rotation(M, ord, pow, ax)) {
            continue;
          }
          if (ord <= 1) {
            neu.kind = OperationKind::Identity;
            neu.order = 1;
          } else {
            neu.order = ord;
            neu.power = pow;
            neu.axis = ax;
          }
        } else if (std::abs(det + 1.0) < 1e-6) {
          // Improper: reflection if M^2 = I and trace = 1
          if (matrices_equal(M * M, Mat3::Identity(), matrix_tol) && std::abs(tr - 1.0) < 1e-4) {
            neu.kind = OperationKind::Reflection;
            neu.order = 1;
            Eigen::SelfAdjointEigenSolver<Mat3> es(0.5 * (M + M.transpose()));
            neu.axis = Vec3::UnitZ();
            for (int k = 0; k < 3; ++k) {
              if (std::abs(es.eigenvalues()(k) + 1.0) < 1e-4) {
                neu.axis = es.eigenvectors().col(k);
                break;
              }
            }
          } else {
            neu.kind = OperationKind::ImproperRotation;
            Vec3 ax(M(2, 1) - M(1, 2), M(0, 2) - M(2, 0), M(1, 0) - M(0, 1));
            if (ax.norm() < 1e-8) {
              Eigen::SelfAdjointEigenSolver<Mat3> es(0.5 * (M + M.transpose()));
              for (int k = 0; k < 3; ++k) {
                if (std::abs(es.eigenvalues()(k) + 1.0) < 0.5) {
                  ax = es.eigenvectors().col(k);
                  break;
                }
              }
            }
            if (ax.norm() > 1e-12) {
              neu.axis = ax.normalized();
            } else {
              neu.axis = Vec3::UnitZ();
            }
            // Matrix order of S, then Schönflies power from the proper part σ·S.
            neu.order = 2;
            neu.power = 1;
            Mat3 P = Mat3::Identity();
            for (int ord = 1; ord <= 48; ++ord) {
              P = M * P;
              if (matrices_equal(P, Mat3::Identity(), matrix_tol)) {
                neu.order = ord;
                break;
              }
            }
            Mat3 proper = reflection_matrix(neu.axis) * M;
            int prop_ord = 0;
            int prop_pow = 0;
            Vec3 prop_ax = neu.axis;
            if (classify_proper_rotation(proper, prop_ord, prop_pow, prop_ax) && prop_ord >= 2 &&
                neu.order % prop_ord == 0) {
              neu.power = prop_pow * (neu.order / prop_ord);
            }
          }
        } else {
          continue;
        }
        ops.push_back(neu);
        grew = true;
        break;
      }
    }
  }
}

void assign_conjugacy_classes(std::vector<Operation>& ops, double matrix_tol) {
  const int n = static_cast<int>(ops.size());
  std::vector<int> cls(n, -1);
  int next = 0;
  for (int i = 0; i < n; ++i) {
    if (cls[static_cast<std::size_t>(i)] >= 0) {
      continue;
    }
    cls[static_cast<std::size_t>(i)] = next;
    for (int g = 0; g < n; ++g) {
      Mat3 conj = ops[static_cast<std::size_t>(g)].matrix *
                  ops[static_cast<std::size_t>(i)].matrix *
                  ops[static_cast<std::size_t>(g)].matrix.transpose();
      // For orthogonal matrices, g^{-1} = g^T
      for (int j = 0; j < n; ++j) {
        if (cls[static_cast<std::size_t>(j)] >= 0) {
          continue;
        }
        if (matrices_equal(ops[static_cast<std::size_t>(j)].matrix, conj, matrix_tol)) {
          cls[static_cast<std::size_t>(j)] = next;
        }
      }
    }
    ++next;
  }
  for (int i = 0; i < n; ++i) {
    ops[static_cast<std::size_t>(i)].class_index = cls[static_cast<std::size_t>(i)];
  }
}

std::vector<EquivalenceSet> orbits_from_operations(const std::vector<Operation>& ops, int n_atoms) {
  std::vector<int> label(static_cast<std::size_t>(n_atoms), -1);
  int next = 0;
  for (int i = 0; i < n_atoms; ++i) {
    if (label[static_cast<std::size_t>(i)] >= 0) {
      continue;
    }
    label[static_cast<std::size_t>(i)] = next;
    bool grew = true;
    while (grew) {
      grew = false;
      for (int a = 0; a < n_atoms; ++a) {
        if (label[static_cast<std::size_t>(a)] != next) {
          continue;
        }
        for (const Operation& op : ops) {
          if (static_cast<int>(op.image.size()) != n_atoms) {
            continue;
          }
          const int j = op.image[static_cast<std::size_t>(a)];
          if (j >= 0 && j < n_atoms && label[static_cast<std::size_t>(j)] != next) {
            label[static_cast<std::size_t>(j)] = next;
            grew = true;
          }
        }
      }
    }
    ++next;
  }
  std::vector<EquivalenceSet> sets(static_cast<std::size_t>(next));
  for (int i = 0; i < n_atoms; ++i) {
    sets[static_cast<std::size_t>(label[static_cast<std::size_t>(i)])].atoms.push_back(i);
  }
  return sets;
}

} // namespace molsymm::detail
