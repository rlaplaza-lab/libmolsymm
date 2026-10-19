#include "molsymm/analysis.hpp"

#include "internal.hpp"
#include "molsymm/error.hpp"

#include <algorithm>
#include <cmath>

namespace molsymm {

double Analysis::symmetrize() {
  if (molecule.size() == 0) {
    return 0.0;
  }
  if (point_group == "Kh") {
    return 0.0;
  }

  const int n = molecule.size();
  std::vector<Vec3> before(n);
  for (int i = 0; i < n; ++i) {
    before[static_cast<std::size_t>(i)] = molecule.atoms[static_cast<std::size_t>(i)].position;
  }

  if (linear || point_group == "Cinfv" || point_group == "Dinfh") {
    for (int i = 0; i < n; ++i) {
      Vec3& p = molecule.atoms[static_cast<std::size_t>(i)].position;
      p.x() = 0.0;
      p.y() = 0.0;
    }
    if (point_group == "Dinfh") {
      const Operation* inv = nullptr;
      for (const Operation& op : operations) {
        if (op.kind == OperationKind::Inversion && static_cast<int>(op.image.size()) == n) {
          inv = &op;
          break;
        }
      }
      std::vector<char> done(n, 0);
      for (int i = 0; i < n; ++i) {
        if (done[static_cast<std::size_t>(i)]) {
          continue;
        }
        const Vec3 pi = molecule.atoms[static_cast<std::size_t>(i)].position;
        int partner = -1;
        if (inv != nullptr) {
          partner = inv->image[static_cast<std::size_t>(i)];
          if (partner < 0 || partner >= n ||
              molecule.atoms[static_cast<std::size_t>(partner)].atomic_number !=
                  molecule.atoms[static_cast<std::size_t>(i)].atomic_number) {
            partner = -1;
          }
        }
        if (partner < 0) {
          double best = 1e9;
          for (int j = 0; j < n; ++j) {
            if (molecule.atoms[static_cast<std::size_t>(j)].atomic_number !=
                molecule.atoms[static_cast<std::size_t>(i)].atomic_number) {
              continue;
            }
            const double d = (molecule.atoms[static_cast<std::size_t>(j)].position + pi).norm();
            if (d < best) {
              best = d;
              partner = j;
            }
          }
        }
        if (partner >= 0) {
          const double z =
              0.5 * (std::abs(pi.z()) +
                     std::abs(molecule.atoms[static_cast<std::size_t>(partner)].position.z()));
          const double s = (pi.z() >= 0) ? 1.0 : -1.0;
          molecule.atoms[static_cast<std::size_t>(i)].position = Vec3(0, 0, s * z);
          molecule.atoms[static_cast<std::size_t>(partner)].position = Vec3(0, 0, -s * z);
          done[static_cast<std::size_t>(i)] = 1;
          done[static_cast<std::size_t>(partner)] = 1;
        }
      }
    }
  } else if (!operations.empty()) {
    // r'_i = (1/|G|) Σ_g R_g r_{g^{-1} i}
    const double inv_g = 1.0 / static_cast<double>(operations.size());
    std::vector<Vec3> avg(n, Vec3::Zero());
    std::vector<char> seen(static_cast<std::size_t>(n), 0);
    for (const Operation& op : operations) {
      if (static_cast<int>(op.image.size()) != n) {
        throw Error("operation image size mismatch");
      }
      std::fill(seen.begin(), seen.end(), 0);
      for (int j = 0; j < n; ++j) {
        const int i = op.image[static_cast<std::size_t>(j)];
        if (i < 0 || i >= n || seen[static_cast<std::size_t>(i)]) {
          throw Error("operation image is not a permutation");
        }
        seen[static_cast<std::size_t>(i)] = 1;
        avg[static_cast<std::size_t>(i)] += op.matrix * before[static_cast<std::size_t>(j)];
      }
    }
    for (int i = 0; i < n; ++i) {
      molecule.atoms[static_cast<std::size_t>(i)].position =
          inv_g * avg[static_cast<std::size_t>(i)];
    }
  }

  double sum = 0.0;
  for (int i = 0; i < n; ++i) {
    sum +=
        (molecule.atoms[static_cast<std::size_t>(i)].position - before[static_cast<std::size_t>(i)])
            .squaredNorm();
  }
  rmsd = std::sqrt(sum / std::max(n, 1));
  return rmsd;
}

} // namespace molsymm
