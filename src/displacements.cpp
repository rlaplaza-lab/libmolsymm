#include "molsymm/displacement.hpp"

#include "internal.hpp"
#include "molsymm/error.hpp"

#include <algorithm>
#include <cmath>

namespace molsymm {
namespace {

bool dirs_same_orbit(const Vec3& a, const Vec3& b, double tol) {
  return (a - b).norm() <= tol || (a + b).norm() <= tol;
}

std::vector<Operation> stabilizer(const Analysis& analysis, int atom) {
  std::vector<Operation> stab;
  for (const Operation& op : analysis.operations) {
    if (static_cast<int>(op.image.size()) != analysis.molecule.size()) {
      continue;
    }
    if (op.image[static_cast<std::size_t>(atom)] == atom) {
      stab.push_back(op);
    }
  }
  return stab;
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

/// One representative direction per site-symmetry invariant subspace of R^3.
std::vector<Vec3> site_unique_directions(const Analysis& analysis, int atom) {
  const Vec3 pos = analysis.molecule.atoms[static_cast<std::size_t>(atom)].position;
  auto stab = stabilizer(analysis, atom);

  std::vector<Vec3> candidates;
  if (pos.norm() > 1e-8) {
    candidates.push_back(pos.normalized());
  }
  candidates.emplace_back(Vec3::UnitX());
  candidates.emplace_back(Vec3::UnitY());
  candidates.emplace_back(Vec3::UnitZ());
  if (pos.norm() > 1e-8) {
    Vec3 radial = pos.normalized();
    for (int c = 0; c < 3; ++c) {
      Vec3 axis = Vec3::Zero();
      axis(c) = 1.0;
      Vec3 perp = radial.cross(axis);
      if (perp.norm() > 1e-8) {
        candidates.push_back(perp.normalized());
      }
    }
  }

  auto in_span = [](const Vec3& v, const std::vector<Vec3>& basis) {
    Vec3 r = v;
    for (const Vec3& b : basis) {
      r -= b.dot(r) * b;
    }
    return r.norm() < 1e-8;
  };

  std::vector<Vec3> used_basis; // orthonormal basis of covered subspaces
  std::vector<Vec3> unique;
  for (Vec3 cand : candidates) {
    if (cand.norm() < 1e-10) {
      continue;
    }
    cand.normalize();
    if (in_span(cand, used_basis)) {
      continue;
    }

    // Invariant subspace = span{ R_g cand }
    std::vector<Vec3> subspace;
    for (const Operation& op : stab) {
      Vec3 v = op.matrix * cand;
      if (v.norm() < 1e-10) {
        continue;
      }
      v.normalize();
      Vec3 r = v;
      for (const Vec3& b : subspace) {
        r -= b.dot(r) * b;
      }
      if (r.norm() > 1e-8) {
        subspace.push_back(r.normalized());
      }
    }
    if (subspace.empty()) {
      subspace.push_back(cand);
    }

    Vec3 rep = subspace.front();
    canonicalize(rep);
    unique.push_back(rep);
    for (const Vec3& b : subspace) {
      Vec3 r = b;
      for (const Vec3& u : used_basis) {
        r -= u.dot(r) * u;
      }
      if (r.norm() > 1e-8) {
        used_basis.push_back(r.normalized());
      }
    }
    if (used_basis.size() >= 3) {
      break;
    }
  }
  return unique;
}

std::string preserved_subgroup_of(const Analysis& analysis, int atom, const Vec3& dir) {
  bool has_sigma = false;
  bool has_inv = false;
  int max_cn = 1;
  for (const Operation& op : analysis.operations) {
    if (op.kind == OperationKind::Identity) {
      continue;
    }
    if (static_cast<int>(op.image.size()) != analysis.molecule.size()) {
      continue;
    }
    if (op.image[static_cast<std::size_t>(atom)] != atom) {
      continue;
    }
    Vec3 Rd = op.matrix * dir;
    if (!dirs_same_orbit(Rd.normalized(), dir, 1e-6)) {
      continue;
    }
    if (op.kind == OperationKind::Reflection) {
      has_sigma = true;
    }
    if (op.kind == OperationKind::Inversion) {
      has_inv = true;
    }
    if (op.kind == OperationKind::ProperRotation) {
      max_cn = std::max(max_cn, op.order);
    }
  }
  if (analysis.point_group == "C1") {
    return "C1";
  }
  if (has_inv && has_sigma) {
    return "C2h";
  }
  if (has_inv) {
    return "Ci";
  }
  if (has_sigma && max_cn >= 2) {
    return "C" + std::to_string(max_cn) + "v";
  }
  if (has_sigma) {
    return "Cs";
  }
  if (max_cn >= 2) {
    return "C" + std::to_string(max_cn);
  }
  return "C1";
}

std::vector<std::string> project_irreps(const Analysis& analysis, int atom, const Vec3& dir) {
  std::vector<std::string> out;
  const CharacterTable& ct = analysis.character_table;
  if (ct.irrep_names.empty() || analysis.operations.empty()) {
    throw Error("character table missing for irrep projection");
  }
  const int g = static_cast<int>(analysis.operations.size());
  for (std::size_t ir = 0; ir < ct.irrep_names.size(); ++ir) {
    double proj = 0.0;
    for (const Operation& op : analysis.operations) {
      if (op.class_index < 0 || op.class_index >= static_cast<int>(ct.characters[ir].size())) {
        continue;
      }
      if (static_cast<int>(op.image.size()) != analysis.molecule.size()) {
        continue;
      }
      if (op.image[static_cast<std::size_t>(atom)] != atom) {
        continue;
      }
      const double chi_disp = dir.dot(op.matrix * dir);
      proj += ct.characters[ir][static_cast<std::size_t>(op.class_index)] * chi_disp;
    }
    proj *= static_cast<double>(ct.irrep_dims[ir]) / static_cast<double>(std::max(g, 1));
    if (std::abs(proj) > 1e-6) {
      out.push_back(ct.irrep_names[ir]);
    }
  }
  return out;
}

Molecule to_original_frame(const Analysis& analysis, const Molecule& frame_mol) {
  Molecule out = analysis.original;
  const Mat3 Rinv = analysis.orientation.transpose();
  for (int i = 0; i < out.size(); ++i) {
    out.atoms[static_cast<std::size_t>(i)].position =
        Rinv * frame_mol.atoms[static_cast<std::size_t>(i)].position + analysis.center_of_mass;
  }
  return out;
}

Eigen::VectorXd
rotate_pattern(const Analysis& analysis, const Operation& op, const Eigen::VectorXd& p) {
  const int n = analysis.molecule.size();
  Eigen::VectorXd out = Eigen::VectorXd::Zero(3 * n);
  for (int i = 0; i < n; ++i) {
    const int j = op.image[static_cast<std::size_t>(i)];
    out.segment<3>(3 * j) = op.matrix * p.segment<3>(3 * i);
  }
  return out;
}

bool pattern_seen(const Eigen::VectorXd& e, const std::vector<Eigen::VectorXd>& seen) {
  return std::any_of(seen.begin(), seen.end(), [&](const Eigen::VectorXd& s) {
    return (e - s).norm() < 1e-6 || (e + s).norm() < 1e-6;
  });
}

} // namespace

DisplacementSet unique_displacements(const Analysis& analysis, double step) {
  DisplacementSet set;
  set.step = step;
  set.analysis = analysis;

  if (analysis.point_group == "Kh" || analysis.molecule.size() == 0) {
    return set;
  }

  const int n = analysis.molecule.size();

  if (analysis.linear || analysis.point_group == "Cinfv" || analysis.point_group == "Dinfh") {
    std::vector<char> seen(n, 0);
    for (int i = 0; i < n; ++i) {
      if (seen[static_cast<std::size_t>(i)]) {
        continue;
      }
      const Vec3 pi = analysis.molecule.atoms[static_cast<std::size_t>(i)].position;
      for (int j = 0; j < n; ++j) {
        if (analysis.molecule.atoms[static_cast<std::size_t>(j)].atomic_number !=
            analysis.molecule.atoms[static_cast<std::size_t>(i)].atomic_number) {
          continue;
        }
        const Vec3 pj = analysis.molecule.atoms[static_cast<std::size_t>(j)].position;
        if (std::abs(std::abs(pj.z()) - std::abs(pi.z())) < 1e-6) {
          seen[static_cast<std::size_t>(j)] = 1;
        }
      }
      for (const Vec3& dir : {Vec3(Vec3::UnitZ()), Vec3(Vec3::UnitX())}) {
        UniqueDisplacement ud;
        ud.atom = i;
        ud.direction = dir;
        ud.preserved_subgroup = (dir.z() > 0.5) ? analysis.point_group : "Cs";
        if (analysis.point_group == "Dinfh") {
          ud.irreps = {(dir.z() > 0.5) ? "Σg+" : "Πu"};
        } else {
          ud.irreps = {(dir.z() > 0.5) ? "Σ+" : "Π"};
        }
        for (std::size_t oi = 0; oi < analysis.operations.size(); ++oi) {
          const Operation& op = analysis.operations[oi];
          if (static_cast<int>(op.image.size()) != n) {
            continue;
          }
          DisplacementImage im;
          im.operation_index = static_cast<int>(oi);
          im.target_atom = op.image[static_cast<std::size_t>(i)];
          im.rotation = op.matrix;
          ud.images.push_back(im);
        }
        set.displacements.push_back(ud);
      }
    }
  } else if (analysis.point_group == "C1") {
    for (int i = 0; i < n; ++i) {
      for (int c = 0; c < 3; ++c) {
        UniqueDisplacement ud;
        ud.atom = i;
        ud.direction = Vec3::Zero();
        ud.direction(c) = 1.0;
        ud.preserved_subgroup = "C1";
        ud.irreps = {"A"};
        DisplacementImage im;
        im.operation_index = 0;
        im.target_atom = i;
        im.rotation = Mat3::Identity();
        ud.images.push_back(im);
        set.displacements.push_back(ud);
      }
    }
  } else {
    for (const EquivalenceSet& es : analysis.equivalence_sets) {
      const int a = *std::min_element(es.atoms.begin(), es.atoms.end());
      for (const Vec3& dir : site_unique_directions(analysis, a)) {
        UniqueDisplacement ud;
        ud.atom = a;
        ud.direction = dir;
        ud.irreps = project_irreps(analysis, a, dir);
        ud.preserved_subgroup = preserved_subgroup_of(analysis, a, dir);
        for (std::size_t oi = 0; oi < analysis.operations.size(); ++oi) {
          const Operation& op = analysis.operations[oi];
          if (static_cast<int>(op.image.size()) != n) {
            continue;
          }
          DisplacementImage im;
          im.operation_index = static_cast<int>(oi);
          im.target_atom = op.image[static_cast<std::size_t>(a)];
          im.rotation = op.matrix;
          ud.images.push_back(im);
        }
        set.displacements.push_back(ud);
      }
    }
  }

  // ± geometries: displace ONLY the representative atom (unique Cartesian),
  // returned in the original input frame.
  for (const UniqueDisplacement& ud : set.displacements) {
    auto displaced = [&](double signed_step) {
      Molecule m = analysis.molecule;
      m.atoms[static_cast<std::size_t>(ud.atom)].position += signed_step * ud.direction;
      return to_original_frame(analysis, m);
    };
    set.plus.push_back(displaced(step));
    set.minus.push_back(displaced(-step));
  }

  return set;
}

Eigen::VectorXd DisplacementSet::gradient_from_energies(std::span<const double> e_plus,
                                                        std::span<const double> e_minus) const {
  const int n = analysis.original.size();
  if (static_cast<int>(e_plus.size()) != static_cast<int>(displacements.size()) ||
      static_cast<int>(e_minus.size()) != static_cast<int>(displacements.size())) {
    throw Error("energy array size mismatch");
  }

  // Build gradient in analysis frame, scatter each unique derivative with R
  Eigen::VectorXd g_frame = Eigen::VectorXd::Zero(3 * n);
  for (std::size_t u = 0; u < displacements.size(); ++u) {
    const double dE = (e_plus[u] - e_minus[u]) / (2.0 * step);
    const UniqueDisplacement& ud = displacements[u];
    std::vector<char> hit(n, 0);
    for (const DisplacementImage& im : ud.images) {
      if (hit[static_cast<std::size_t>(im.target_atom)]) {
        continue;
      }
      hit[static_cast<std::size_t>(im.target_atom)] = 1;
      g_frame.segment<3>(3 * im.target_atom) += dE * (im.rotation * ud.direction);
    }
  }

  Eigen::VectorXd g = Eigen::VectorXd::Zero(3 * n);
  const Mat3 Rinv = analysis.orientation.transpose();
  for (int i = 0; i < n; ++i) {
    g.segment<3>(3 * i) = Rinv * g_frame.segment<3>(3 * i);
  }
  return g;
}

Eigen::MatrixXd
DisplacementSet::hessian_from_gradients(std::span<const Eigen::VectorXd> g_plus,
                                        std::span<const Eigen::VectorXd> g_minus) const {
  const int n = analysis.original.size();
  const int n3 = 3 * n;
  if (static_cast<int>(g_plus.size()) != static_cast<int>(displacements.size()) ||
      static_cast<int>(g_minus.size()) != static_cast<int>(displacements.size())) {
    throw Error("gradient array size mismatch");
  }

  const Mat3 R = analysis.orientation;
  const Mat3 Rinv = R.transpose();
  std::vector<Eigen::VectorXd> e_cols;
  std::vector<Eigen::VectorXd> c_cols;

  for (std::size_t u = 0; u < displacements.size(); ++u) {
    if (g_plus[u].size() != n3 || g_minus[u].size() != n3) {
      throw Error("gradient vector length mismatch");
    }
    Eigen::VectorXd col_orig = (g_plus[u] - g_minus[u]) / (2.0 * step);
    Eigen::VectorXd col = Eigen::VectorXd::Zero(n3);
    for (int i = 0; i < n; ++i) {
      col.segment<3>(3 * i) = R * col_orig.segment<3>(3 * i);
    }

    Eigen::VectorXd e = Eigen::VectorXd::Zero(n3);
    e.segment<3>(3 * displacements[u].atom) = displacements[u].direction;

    for (const Operation& op : analysis.operations) {
      if (static_cast<int>(op.image.size()) != n) {
        continue;
      }
      Eigen::VectorXd e_img = rotate_pattern(analysis, op, e);
      const double nrm = e_img.norm();
      if (nrm < 1e-12) {
        continue;
      }
      e_img /= nrm;
      if (pattern_seen(e_img, e_cols)) {
        continue;
      }
      e_cols.push_back(e_img);
      c_cols.emplace_back(rotate_pattern(analysis, op, col) / nrm);
    }
  }

  Eigen::MatrixXd H_f = Eigen::MatrixXd::Zero(n3, n3);
  if (!e_cols.empty()) {
    const int m = static_cast<int>(e_cols.size());
    Eigen::MatrixXd E(n3, m);
    Eigen::MatrixXd C(n3, m);
    for (int k = 0; k < m; ++k) {
      E.col(k) = e_cols[static_cast<std::size_t>(k)];
      C.col(k) = c_cols[static_cast<std::size_t>(k)];
    }
    // H E = C  =>  H = C pinv(E)
    H_f = C * E.completeOrthogonalDecomposition().pseudoInverse();
    H_f = 0.5 * (H_f + H_f.transpose());
  }

  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(n3, n3);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      H.block<3, 3>(3 * i, 3 * j) = Rinv * H_f.block<3, 3>(3 * i, 3 * j) * R;
    }
  }
  return H;
}

Eigen::VectorXd symmetrize_gradient(const Analysis& analysis, const Eigen::VectorXd& gradient) {
  const int n = analysis.original.size();
  if (gradient.size() != 3 * n) {
    throw Error("gradient size mismatch");
  }
  if (analysis.operations.empty()) {
    return gradient;
  }

  const Mat3 R = analysis.orientation;
  const Mat3 Rinv = R.transpose();
  Eigen::VectorXd g_f = Eigen::VectorXd::Zero(3 * n);
  for (int i = 0; i < n; ++i) {
    g_f.segment<3>(3 * i) = R * gradient.segment<3>(3 * i);
  }

  Eigen::VectorXd avg = Eigen::VectorXd::Zero(3 * n);
  int nop = 0;
  for (const Operation& op : analysis.operations) {
    if (static_cast<int>(op.image.size()) != n) {
      continue;
    }
    ++nop;
    for (int i = 0; i < n; ++i) {
      const int j = op.image[static_cast<std::size_t>(i)];
      avg.segment<3>(3 * j) += op.matrix * g_f.segment<3>(3 * i);
    }
  }
  if (nop > 0) {
    avg /= static_cast<double>(nop);
  }
  Eigen::VectorXd out = Eigen::VectorXd::Zero(3 * n);
  for (int i = 0; i < n; ++i) {
    out.segment<3>(3 * i) = Rinv * avg.segment<3>(3 * i);
  }
  return out;
}

Eigen::MatrixXd symmetrize_hessian(const Analysis& analysis, const Eigen::MatrixXd& hessian) {
  const int n = analysis.original.size();
  const int n3 = 3 * n;
  if (hessian.rows() != n3 || hessian.cols() != n3) {
    throw Error("hessian size mismatch");
  }
  if (analysis.operations.empty()) {
    return hessian;
  }

  const Mat3 R = analysis.orientation;
  const Mat3 Rinv = R.transpose();
  Eigen::MatrixXd H_f = Eigen::MatrixXd::Zero(n3, n3);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      H_f.block<3, 3>(3 * i, 3 * j) = R * hessian.block<3, 3>(3 * i, 3 * j) * Rinv;
    }
  }

  Eigen::MatrixXd avg = Eigen::MatrixXd::Zero(n3, n3);
  int nop = 0;
  for (const Operation& op : analysis.operations) {
    if (static_cast<int>(op.image.size()) != n) {
      continue;
    }
    ++nop;
    for (int i = 0; i < n; ++i) {
      const int ii = op.image[static_cast<std::size_t>(i)];
      for (int j = 0; j < n; ++j) {
        const int jj = op.image[static_cast<std::size_t>(j)];
        avg.block<3, 3>(3 * ii, 3 * jj) +=
            op.matrix * H_f.block<3, 3>(3 * i, 3 * j) * op.matrix.transpose();
      }
    }
  }
  if (nop > 0) {
    avg /= static_cast<double>(nop);
  }
  Eigen::MatrixXd out = Eigen::MatrixXd::Zero(n3, n3);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      out.block<3, 3>(3 * i, 3 * j) = Rinv * avg.block<3, 3>(3 * i, 3 * j) * R;
    }
  }
  return out;
}

Eigen::MatrixXd project_external(const Analysis& analysis, const Eigen::MatrixXd& hessian) {
  const int n = analysis.original.size();
  const int n3 = 3 * n;
  if (hessian.rows() != n3 || hessian.cols() != n3) {
    throw Error("hessian size mismatch");
  }

  std::vector<Eigen::VectorXd> modes;
  for (int c = 0; c < 3; ++c) {
    Eigen::VectorXd t = Eigen::VectorXd::Zero(n3);
    for (int i = 0; i < n; ++i) {
      t(3 * i + c) = 1.0;
    }
    if (t.norm() > 1e-12) {
      modes.push_back(t.normalized());
    }
  }

  const Vec3 com = analysis.center_of_mass;
  auto rot_mode = [&](const Vec3& axis) {
    Eigen::VectorXd r = Eigen::VectorXd::Zero(n3);
    for (int i = 0; i < n; ++i) {
      Vec3 pos = analysis.original.atoms[static_cast<std::size_t>(i)].position - com;
      r.segment<3>(3 * i) = axis.cross(pos);
    }
    return r;
  };
  if (analysis.linear) {
    for (const Vec3& ax : {Vec3(Vec3::UnitX()), Vec3(Vec3::UnitY())}) {
      Eigen::VectorXd r = rot_mode(ax);
      if (r.norm() > 1e-8) {
        modes.push_back(r.normalized());
      }
    }
  } else {
    for (int c = 0; c < 3; ++c) {
      Eigen::VectorXd r = rot_mode(Vec3::Unit(c));
      if (r.norm() > 1e-8) {
        modes.push_back(r.normalized());
      }
    }
  }

  std::vector<Eigen::VectorXd> basis;
  for (Eigen::VectorXd v : modes) {
    for (const Eigen::VectorXd& b : basis) {
      v -= b.dot(v) * b;
    }
    if (v.norm() > 1e-8) {
      basis.push_back(v.normalized());
    }
  }

  Eigen::MatrixXd H = hessian;
  for (const Eigen::VectorXd& b : basis) {
    H = H - (H * b) * b.transpose() - b * (b.transpose() * H) +
        b * (b.transpose() * H * b) * b.transpose();
  }
  return H;
}

} // namespace molsymm
