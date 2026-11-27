#pragma once

#include "molsymm/molsymm.hpp"

#include <algorithm>
#include <array>
#include <random>
#include <utility>
#include <vector>

namespace molsymm::test {

inline Mat3 rot(const Vec3& axis, double ang) {
  return Eigen::AngleAxisd(ang, axis.normalized()).toRotationMatrix();
}

inline Mat3 refl(const Vec3& n) {
  const Vec3 u = n.normalized();
  return Mat3::Identity() - 2.0 * u * u.transpose();
}

inline bool mat_near(const Mat3& a, const Mat3& b, double tol = 1e-8) {
  return (a - b).cwiseAbs().maxCoeff() < tol;
}

inline std::vector<Mat3> close_group(std::vector<Mat3> g) {
  g.insert(g.begin(), Mat3::Identity());
  bool grew = true;
  while (grew) {
    grew = false;
    const std::size_t n = g.size();
    for (std::size_t i = 0; i < n; ++i) {
      for (std::size_t j = 0; j < n; ++j) {
        const Mat3 p = g[i] * g[j];
        if (std::none_of(g.begin(), g.end(), [&](const Mat3& m) { return mat_near(m, p); })) {
          g.push_back(p);
          grew = true;
        }
      }
    }
    if (g.size() > 240) {
      break;
    }
  }
  return g;
}

// Generic seeds with distinct Z; small |G| needs more orbits to avoid accidental symmetry.
inline Molecule molecule_of(const std::vector<Mat3>& g,
                            const Mat3& frame,
                            std::mt19937* shuffle, // NOLINT(readability-non-const-parameter)
                            bool axis_atoms = true) {
  const std::array<Vec3, 4> seeds{Vec3(1.13, 0.37, 0.71),
                                  Vec3(0.52, 1.41, -0.33),
                                  Vec3(-0.87, 0.29, 1.52),
                                  Vec3(0.41, -1.07, -1.23)};
  const std::array<int, 4> zs{6, 7, 8, 9};
  const std::size_t norb = g.size() <= 12 ? 4 : 2;
  std::vector<std::pair<int, Vec3>> atoms;
  for (std::size_t s = 0; s < norb; ++s) {
    for (const Mat3& m : g) {
      atoms.emplace_back(zs[s], frame * (m * seeds[s]));
    }
  }
  if (axis_atoms) {
    std::vector<Vec3> dirs;
    for (const Mat3& m : g) {
      const Mat3 r = (m.determinant() > 0) ? m : Mat3(-m);
      if (mat_near(r, Mat3::Identity())) {
        continue;
      }
      const Eigen::AngleAxisd aa(r);
      for (const Vec3& d : {Vec3(aa.axis()), Vec3(-aa.axis())}) {
        if (std::none_of(
                dirs.begin(), dirs.end(), [&](const Vec3& q) { return (q - d).norm() < 1e-6; })) {
          dirs.push_back(d);
        }
      }
    }
    for (const Vec3& d : dirs) {
      atoms.emplace_back(1, frame * (1.9 * d));
    }
  }
  if (shuffle != nullptr) {
    std::shuffle(atoms.begin(), atoms.end(), *shuffle);
  }
  std::vector<int> z;
  std::vector<double> xyz;
  for (const auto& [zi, p] : atoms) {
    z.push_back(zi);
    xyz.insert(xyz.end(), {p.x(), p.y(), p.z()});
  }
  return Molecule(z, xyz);
}

} // namespace molsymm::test
