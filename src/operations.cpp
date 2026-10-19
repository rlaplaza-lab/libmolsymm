#include "internal.hpp"

#include "molsymm/error.hpp"

#include <algorithm>
#include <cmath>

namespace molsymm::detail {

int find_image(const std::vector<Atom>& atoms, int from, const Mat3& R, double tol) {
  const Vec3 target = R * atoms[static_cast<std::size_t>(from)].position;
  const int z = atoms[static_cast<std::size_t>(from)].atomic_number;
  int best = -1;
  double best_d = tol + 1.0;
  for (int j = 0; j < static_cast<int>(atoms.size()); ++j) {
    if (atoms[static_cast<std::size_t>(j)].atomic_number != z) {
      continue;
    }
    const double d = (atoms[static_cast<std::size_t>(j)].position - target).norm();
    if (d < best_d) {
      best_d = d;
      best = j;
    }
  }
  return best_d <= tol ? best : -1;
}

std::vector<int> permutation(const std::vector<Atom>& atoms, const Mat3& R, double tol) {
  std::vector<int> image(atoms.size(), -1);
  std::vector<char> taken(atoms.size(), 0);
  for (int i = 0; i < static_cast<int>(atoms.size()); ++i) {
    const int j = find_image(atoms, i, R, tol);
    if (j < 0 || taken[static_cast<std::size_t>(j)]) {
      return {};
    }
    taken[static_cast<std::size_t>(j)] = 1;
    image[static_cast<std::size_t>(i)] = j;
  }
  return image;
}

bool is_symmetry(const std::vector<Atom>& atoms, const Mat3& R, double tol) {
  return !permutation(atoms, R, tol).empty();
}

void close_group(std::vector<Operation>& ops,
                 const std::vector<Atom>& atoms,
                 double tol,
                 int max_order) {
  // Ensure identity present
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
    e.image = permutation(atoms, e.matrix, tol);
    ops.push_back(e);
  }

  for (Operation& op : ops) {
    if (op.image.empty()) {
      op.image = permutation(atoms, op.matrix, tol);
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
          if (matrices_equal(op.matrix, M, tol * 10.0)) {
            found = true;
            break;
          }
        }
        if (found) {
          continue;
        }
        if (!is_symmetry(atoms, M, tol)) {
          continue;
        }
        Operation neu;
        neu.matrix = M;
        neu.image = permutation(atoms, M, tol);
        // Classify by matrix properties
        const double det = M.determinant();
        const double tr = M.trace();
        if (matrices_equal(M, Mat3::Identity(), tol * 10.0)) {
          neu.kind = OperationKind::Identity;
          neu.order = 1;
        } else if (matrices_equal(M, -Mat3::Identity(), tol * 10.0)) {
          neu.kind = OperationKind::Inversion;
          neu.order = 1;
        } else if (std::abs(det - 1.0) < 1e-6) {
          neu.kind = OperationKind::ProperRotation;
          // angle from trace: 1 + 2 cos θ
          double c = std::clamp((tr - 1.0) * 0.5, -1.0, 1.0);
          double ang = std::acos(c);
          if (ang < 1e-8) {
            neu.kind = OperationKind::Identity;
            neu.order = 1;
          } else {
            int best_n = 2;
            double best_err = 1e9;
            for (int ord = 2; ord <= 48; ++ord) {
              const double err = std::abs(ang - 2.0 * kPi / ord);
              if (err < best_err) {
                best_err = err;
                best_n = ord;
              }
            }
            neu.order = best_n;
            Eigen::SelfAdjointEigenSolver<Mat3> es(0.5 * (M + M.transpose()));
            // axis is eigenvector of eigenvalue 1
            neu.axis = Vec3::UnitZ();
            for (int k = 0; k < 3; ++k) {
              if (std::abs(es.eigenvalues()(k) - 1.0) < 1e-4) {
                neu.axis = es.eigenvectors().col(k);
                break;
              }
            }
            // Also try (M - M^T) axis extraction
            Vec3 ax(M(2, 1) - M(1, 2), M(0, 2) - M(2, 0), M(1, 0) - M(0, 1));
            if (ax.norm() > 1e-8) {
              neu.axis = ax.normalized();
            }
          }
        } else if (std::abs(det + 1.0) < 1e-6) {
          // Improper: reflection if M^2 = I and trace = 1
          if (matrices_equal(M * M, Mat3::Identity(), tol * 10.0) && std::abs(tr - 1.0) < 1e-4) {
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
            double c = std::clamp((-tr - 1.0) * 0.5, -1.0, 1.0);
            // For Sn: improper = σ * Cn, trace = -1 + 2 cos? Actually S_n = σh * Cn
            // character of improper rotation: -1 + 2 cos θ for the matrix with det -1
            // More carefully: reflection has tr=1; inversion tr=-3; S3 etc.
            // Use angle of the proper part.
            Mat3 proper = -M; // not always right
            // Extract axis from the -1 eigenspace of the reflection-like part
            Eigen::JacobiSVD<Mat3> svd(M - Mat3::Identity());
            (void)svd;
            Vec3 ax(M(2, 1) - M(1, 2), M(0, 2) - M(2, 0), M(1, 0) - M(0, 1));
            if (ax.norm() < 1e-8) {
              // pure reflection or inversion already handled; S2 = i
              Eigen::SelfAdjointEigenSolver<Mat3> es(0.5 * (M + M.transpose()));
              for (int k = 0; k < 3; ++k) {
                if (std::abs(es.eigenvalues()(k) + 1.0) < 0.5) {
                  ax = es.eigenvectors().col(k);
                  break;
                }
              }
            }
            neu.axis = ax.normalized();
            // Estimate order from M^{2k} = I
            neu.order = 2;
            Mat3 P = Mat3::Identity();
            for (int ord = 1; ord <= 48; ++ord) {
              P = M * P;
              if (matrices_equal(P, Mat3::Identity(), tol * 10.0)) {
                neu.order = ord;
                break;
              }
            }
            (void)c;
            (void)proper;
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

void assign_conjugacy_classes(std::vector<Operation>& ops, double tol) {
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
        if (matrices_equal(ops[static_cast<std::size_t>(j)].matrix, conj, tol * 10.0)) {
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

} // namespace molsymm::detail
