#include "molsymm/fragment.hpp"

#include "internal.hpp"
#include "molsymm/error.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace molsymm {
namespace {

void require_finite_group(const Analysis& analysis) {
  if (analysis.infinite) {
    throw Error("asymmetric unit is undefined for infinite point groups");
  }
}

int validate_operations(const std::vector<Operation>& ops, int n) {
  if (ops.empty()) {
    throw Error("analysis has no operations");
  }
  int identity = -1;
  for (std::size_t oi = 0; oi < ops.size(); ++oi) {
    const Operation& op = ops[oi];
    if (static_cast<int>(op.image.size()) != n) {
      throw Error("operation image size mismatch");
    }
    std::vector<char> seen(static_cast<std::size_t>(n), 0);
    for (int i = 0; i < n; ++i) {
      const int j = op.image[static_cast<std::size_t>(i)];
      if (j < 0 || j >= n || seen[static_cast<std::size_t>(j)]) {
        throw Error("operation image is not a permutation");
      }
      seen[static_cast<std::size_t>(j)] = 1;
    }
    if (op.kind == OperationKind::Identity) {
      identity = static_cast<int>(oi);
    }
  }
  if (identity < 0) {
    throw Error("group has no identity");
  }
  return identity;
}

void canonicalize(Vec3& d) {
  for (int k = 0; k < 3; ++k) {
    if (std::abs(d(k)) > 1e-10) {
      if (d(k) < 0) {
        d = -d;
      }
      break;
    }
  }
}

Mat3 site_projector(const std::vector<Operation>& ops, int atom) {
  Mat3 P = Mat3::Zero();
  int count = 0;
  for (const Operation& op : ops) {
    if (op.image[static_cast<std::size_t>(atom)] != atom) {
      continue;
    }
    P += op.matrix;
    ++count;
  }
  return P / static_cast<double>(count);
}

std::vector<Vec3> free_axes_of(const Mat3& projector, const Mat3& orientation) {
  Eigen::SelfAdjointEigenSolver<Mat3> es(0.5 * (projector + projector.transpose()));
  std::vector<Vec3> axes;
  for (int k = 0; k < 3; ++k) {
    if (es.eigenvalues()(k) <= 0.5) {
      continue;
    }
    Vec3 out = orientation.transpose() * es.eigenvectors().col(k).normalized();
    out.normalize();
    canonicalize(out);
    axes.push_back(out);
  }
  return axes;
}

FragmentSite make_site(const std::vector<Operation>& ops, int atom, const Mat3& orientation) {
  FragmentSite site;
  site.atom = atom;
  site.projector = site_projector(ops, atom);
  site.free_axes = free_axes_of(site.projector, orientation);
  return site;
}

int image_of(const Operation& op, int atom) {
  return op.image[static_cast<std::size_t>(atom)];
}

std::vector<int> representatives_of(const std::vector<EquivalenceSet>& orbits) {
  std::vector<int> reps;
  reps.reserve(orbits.size());
  for (const EquivalenceSet& es : orbits) {
    if (!es.atoms.empty()) {
      reps.push_back(*std::min_element(es.atoms.begin(), es.atoms.end()));
    }
  }
  std::sort(reps.begin(), reps.end());
  return reps;
}

int generic_orbit_index(const std::vector<EquivalenceSet>& orbits) {
  int best = -1;
  int best_size = -1;
  int best_rep = std::numeric_limits<int>::max();
  for (int i = 0; i < static_cast<int>(orbits.size()); ++i) {
    const auto& atoms = orbits[static_cast<std::size_t>(i)].atoms;
    if (atoms.empty()) {
      continue;
    }
    const int size = static_cast<int>(atoms.size());
    const int rep = *std::min_element(atoms.begin(), atoms.end());
    if (size > best_size || (size == best_size && rep < best_rep)) {
      best = i;
      best_size = size;
      best_rep = rep;
    }
  }
  return best;
}

Molecule molecule_shell(const std::vector<int>& atomic_numbers, const std::vector<double>& masses) {
  Molecule mol;
  mol.atoms.resize(atomic_numbers.size());
  for (std::size_t i = 0; i < atomic_numbers.size(); ++i) {
    mol.atoms[i].atomic_number = atomic_numbers[i];
    mol.atoms[i].mass = masses[i];
  }
  return mol;
}

void write_geometry(Molecule& mol,
                    const std::vector<Operation>& ops,
                    const std::vector<int>& reps,
                    const std::vector<Vec3>& frame_reps,
                    const Mat3& orientation,
                    const Vec3& com) {
  const int n = mol.size();
  std::vector<char> written(static_cast<std::size_t>(n), 0);
  for (std::size_t s = 0; s < reps.size(); ++s) {
    const Vec3& p = frame_reps[s];
    for (const Operation& op : ops) {
      const int img = image_of(op, reps[s]);
      if (written[static_cast<std::size_t>(img)]) {
        continue;
      }
      mol.atoms[static_cast<std::size_t>(img)].position =
          orientation.transpose() * (op.matrix * p) + com;
      written[static_cast<std::size_t>(img)] = 1;
    }
  }
  if (std::find(written.begin(), written.end(), 0) != written.end()) {
    throw Error("reconstruction did not cover all atoms");
  }
}

int cover_count(const Operation& op,
                const std::vector<int>& reps,
                const std::vector<char>& covered,
                int n) {
  int cover = 0;
  for (int rep : reps) {
    const int img = image_of(op, rep);
    if (img >= 0 && img < n && !covered[static_cast<std::size_t>(img)]) {
      ++cover;
    }
  }
  return cover;
}

template <typename Eligible>
int best_covering_op(const std::vector<Operation>& ops,
                     const std::vector<int>& reps,
                     const std::vector<char>& covered,
                     int n,
                     Eligible&& eligible) {
  int best_op = -1;
  int best_cover = -1;
  for (std::size_t oi = 0; oi < ops.size(); ++oi) {
    if (!eligible(ops[oi])) {
      continue;
    }
    const int cover = cover_count(ops[oi], reps, covered, n);
    if (cover > best_cover ||
        (cover == best_cover && (best_op < 0 || static_cast<int>(oi) < best_op))) {
      best_cover = cover;
      best_op = static_cast<int>(oi);
    }
  }
  return best_op;
}

std::vector<FragmentCopy> build_copies(const std::vector<Operation>& ops,
                                       const std::vector<int>& reps,
                                       const std::vector<EquivalenceSet>& orbits,
                                       int identity) {
  const int n = static_cast<int>(ops.front().image.size());
  const int gid = generic_orbit_index(orbits);
  if (gid < 0) {
    throw Error("no atom orbits for fragmentation");
  }
  const std::vector<int>& generic = orbits[static_cast<std::size_t>(gid)].atoms;
  const int generic_rep = *std::min_element(generic.begin(), generic.end());

  std::vector<char> covered(static_cast<std::size_t>(n), 0);
  std::vector<FragmentCopy> copies;

  auto append_copy = [&](int op_index) {
    const Operation& op = ops[static_cast<std::size_t>(op_index)];
    FragmentCopy copy;
    for (int rep : reps) {
      const int img = image_of(op, rep);
      if (covered[static_cast<std::size_t>(img)]) {
        continue;
      }
      covered[static_cast<std::size_t>(img)] = 1;
      copy.atoms.push_back(img);
    }
    if (!copy.atoms.empty()) {
      copies.push_back(std::move(copy));
    }
  };

  append_copy(identity);

  auto generic_uncovered = [&] {
    return std::any_of(generic.begin(), generic.end(), [&](int a) {
      return !covered[static_cast<std::size_t>(a)];
    });
  };

  while (generic_uncovered()) {
    const int best_op = best_covering_op(ops, reps, covered, n, [&](const Operation& op) {
      return !covered[static_cast<std::size_t>(image_of(op, generic_rep))];
    });
    if (best_op < 0) {
      throw Error("cannot build asymmetric-unit copies");
    }
    append_copy(best_op);
  }

  for (int a = 0; a < n; ++a) {
    if (covered[static_cast<std::size_t>(a)]) {
      continue;
    }
    const int best_op = best_covering_op(ops, reps, covered, n, [&](const Operation& op) {
      return std::any_of(reps.begin(), reps.end(), [&](int rep) { return image_of(op, rep) == a; });
    });
    if (best_op < 0) {
      throw Error("cannot cover atom in fragmentation");
    }
    append_copy(best_op);
  }

  return copies;
}

int parent_copy_of(const Fragmentation& frag, int atom) {
  for (std::size_t c = 0; c < frag.copies.size(); ++c) {
    const auto& atoms = frag.copies[c].atoms;
    if (std::find(atoms.begin(), atoms.end(), atom) != atoms.end()) {
      return static_cast<int>(c);
    }
  }
  return -1;
}

} // namespace

Fragmentation fragmentation(const Analysis& analysis) {
  require_finite_group(analysis);
  const int n = analysis.molecule.size();
  if (n == 0) {
    throw Error("empty molecule");
  }
  const int identity = validate_operations(analysis.operations, n);

  Fragmentation frag;
  frag.center_of_mass = analysis.center_of_mass;
  frag.orientation = analysis.orientation;
  frag.operations = analysis.operations;
  frag.atomic_numbers.resize(static_cast<std::size_t>(n));
  frag.masses.resize(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    frag.atomic_numbers[static_cast<std::size_t>(i)] =
        analysis.molecule.atoms[static_cast<std::size_t>(i)].atomic_number;
    frag.masses[static_cast<std::size_t>(i)] =
        analysis.molecule.atoms[static_cast<std::size_t>(i)].mass;
  }

  const auto orbits = detail::orbits_from_operations(frag.operations, n);
  frag.representatives = representatives_of(orbits);
  frag.sites.reserve(frag.representatives.size());
  for (int atom : frag.representatives) {
    frag.sites.push_back(make_site(frag.operations, atom, frag.orientation));
  }
  frag.mask.assign(static_cast<std::size_t>(n), 0);
  for (int atom : frag.representatives) {
    frag.mask[static_cast<std::size_t>(atom)] = 1;
  }
  frag.copies = build_copies(frag.operations, frag.representatives, orbits, identity);
  return frag;
}

Molecule reconstruct(const Fragmentation& frag, std::span<const Vec3> representative_positions) {
  if (representative_positions.size() != frag.representatives.size()) {
    throw Error("representative position count mismatch");
  }
  Molecule mol = molecule_shell(frag.atomic_numbers, frag.masses);
  std::vector<Vec3> frame_reps(frag.representatives.size());
  for (std::size_t s = 0; s < frag.representatives.size(); ++s) {
    frame_reps[s] = frag.sites[s].projector * frag.input_to_frame(representative_positions[s]);
  }
  write_geometry(mol,
                 frag.operations,
                 frag.representatives,
                 frame_reps,
                 frag.orientation,
                 frag.center_of_mass);
  return mol;
}

SubgroupAssembly subgroup_assembly(const Analysis& analysis, const std::string& subgroup) {
  require_finite_group(analysis);
  if (std::find(analysis.subgroups.begin(), analysis.subgroups.end(), subgroup) ==
      analysis.subgroups.end()) {
    throw Error("subgroup not listed: " + subgroup);
  }

  Fragmentation parent = fragmentation(analysis);
  Analysis sub = analysis;
  sub.select_subgroup(subgroup);

  const int n = sub.molecule.size();
  validate_operations(sub.operations, n);
  const auto h_reps = representatives_of(detail::orbits_from_operations(sub.operations, n));

  SubgroupAssembly assembly;
  assembly.subgroup = subgroup;
  assembly.center_of_mass = analysis.center_of_mass;
  assembly.orientation = analysis.orientation;
  assembly.operations = std::move(sub.operations);
  assembly.atomic_numbers = std::move(parent.atomic_numbers);
  assembly.masses = std::move(parent.masses);

  std::vector<std::vector<int>> by_copy(parent.copies.size());
  for (int atom : h_reps) {
    const int c = parent_copy_of(parent, atom);
    if (c < 0) {
      throw Error("H representative not owned by a parent copy");
    }
    by_copy[static_cast<std::size_t>(c)].push_back(atom);
  }

  for (std::size_t c = 0; c < by_copy.size(); ++c) {
    if (by_copy[c].empty()) {
      continue;
    }
    std::sort(by_copy[c].begin(), by_copy[c].end());
    AssemblyPiece piece;
    piece.parent_copy = static_cast<int>(c);
    piece.atoms = std::move(by_copy[c]);
    for (int atom : piece.atoms) {
      assembly.sites.push_back(make_site(assembly.operations, atom, assembly.orientation));
    }
    assembly.pieces.push_back(std::move(piece));
  }
  if (assembly.pieces.empty()) {
    throw Error("empty subgroup assembly");
  }
  return assembly;
}

Molecule assemble(const SubgroupAssembly& assembly,
                  const std::vector<std::vector<Vec3>>& piece_positions) {
  if (piece_positions.size() != assembly.pieces.size()) {
    throw Error("piece count mismatch");
  }
  std::vector<int> reps;
  std::vector<Vec3> frame_reps;
  std::size_t site = 0;
  for (std::size_t p = 0; p < assembly.pieces.size(); ++p) {
    if (piece_positions[p].size() != assembly.pieces[p].atoms.size()) {
      throw Error("piece position count mismatch");
    }
    for (std::size_t i = 0; i < assembly.pieces[p].atoms.size(); ++i) {
      frame_reps.emplace_back(assembly.sites[site].projector *
                              assembly.input_to_frame(piece_positions[p][i]));
      reps.push_back(assembly.pieces[p].atoms[i]);
      ++site;
    }
  }

  Molecule mol = molecule_shell(assembly.atomic_numbers, assembly.masses);
  write_geometry(
      mol, assembly.operations, reps, frame_reps, assembly.orientation, assembly.center_of_mass);
  return mol;
}

} // namespace molsymm
