#include "molsymm/analysis.hpp"

#include "internal.hpp"
#include "molsymm/error.hpp"

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
    // Snap to z axis; for Dinfh also enforce inversion pairing
    for (int i = 0; i < n; ++i) {
      Vec3& p = molecule.atoms[static_cast<std::size_t>(i)].position;
      p.x() = 0.0;
      p.y() = 0.0;
    }
    if (point_group == "Dinfh") {
      // Average z with inversion partner of same Z
      std::vector<char> done(n, 0);
      for (int i = 0; i < n; ++i) {
        if (done[static_cast<std::size_t>(i)]) {
          continue;
        }
        const Vec3 pi = molecule.atoms[static_cast<std::size_t>(i)].position;
        int partner = -1;
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
    for (const Operation& op : operations) {
      if (static_cast<int>(op.image.size()) != n) {
        continue;
      }
      for (int i = 0; i < n; ++i) {
        // op maps j → i means j = g^{-1} i; image[j]=i
        int src = -1;
        for (int j = 0; j < n; ++j) {
          if (op.image[static_cast<std::size_t>(j)] == i) {
            src = j;
            break;
          }
        }
        if (src < 0) {
          throw Error("operation image is not a permutation");
        }
        avg[static_cast<std::size_t>(i)] += op.matrix * before[static_cast<std::size_t>(src)];
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
