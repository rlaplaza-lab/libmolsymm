#include "internal.hpp"

#include "molsymm/molecule.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>

namespace molsymm {

namespace {

struct Fingerprint {
  int z = 0;
  double r = 0.0;
  std::vector<double> dists;
};

bool fingerprints_match(const Fingerprint& a, const Fingerprint& b, double tol) {
  if (a.z != b.z) {
    return false;
  }
  if (std::abs(a.r - b.r) > tol) {
    return false;
  }
  if (a.dists.size() != b.dists.size()) {
    return false;
  }
  for (std::size_t i = 0; i < a.dists.size(); ++i) {
    if (std::abs(a.dists[i] - b.dists[i]) > tol) {
      return false;
    }
  }
  return true;
}

Fingerprint make_fp(const Molecule& mol, int i) {
  Fingerprint fp;
  fp.z = mol.atoms[static_cast<std::size_t>(i)].atomic_number;
  fp.r = mol.atoms[static_cast<std::size_t>(i)].position.norm();
  fp.dists.reserve(mol.atoms.size());
  for (int j = 0; j < mol.size(); ++j) {
    fp.dists.push_back((mol.atoms[static_cast<std::size_t>(i)].position -
                        mol.atoms[static_cast<std::size_t>(j)].position)
                           .norm());
  }
  std::sort(fp.dists.begin(), fp.dists.end());
  return fp;
}

} // namespace

std::vector<EquivalenceSet> build_equivalence_sets(const Molecule& mol, double tol) {
  const int n = mol.size();
  std::vector<Fingerprint> fps(n);
  for (int i = 0; i < n; ++i) {
    fps[static_cast<std::size_t>(i)] = make_fp(mol, i);
  }

  std::vector<int> label(n, -1);
  int next = 0;
  for (int i = 0; i < n; ++i) {
    if (label[static_cast<std::size_t>(i)] >= 0) {
      continue;
    }
    label[static_cast<std::size_t>(i)] = next;
    for (int j = i + 1; j < n; ++j) {
      if (label[static_cast<std::size_t>(j)] >= 0) {
        continue;
      }
      if (fingerprints_match(
              fps[static_cast<std::size_t>(i)], fps[static_cast<std::size_t>(j)], tol)) {
        label[static_cast<std::size_t>(j)] = next;
      }
    }
    ++next;
  }

  std::vector<EquivalenceSet> sets(next);
  for (int i = 0; i < n; ++i) {
    sets[static_cast<std::size_t>(label[static_cast<std::size_t>(i)])].atoms.push_back(i);
  }
  return sets;
}

} // namespace molsymm
