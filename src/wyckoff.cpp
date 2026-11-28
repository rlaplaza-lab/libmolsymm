#include "molsymm/wyckoff.hpp"

#include "internal.hpp"
#include "molsymm/error.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <sstream>
#include <utility>

namespace molsymm::detail {
namespace {

constexpr double kMatTol = 1e-8;
constexpr int kMaxOrder = 240;

Operation from_matrix(const Mat3& M) {
  Operation op;
  op.matrix = M;
  return op;
}

void retag_operations(std::vector<Operation>& ops) {
  std::vector<Operation> tmp;
  tmp.reserve(ops.size());
  for (const Operation& op : ops) {
    tmp.push_back(from_matrix(op.matrix));
  }
  close_group_abstract(tmp, kMatTol, kMaxOrder);
  ops = std::move(tmp);
}

void conjugate(std::vector<Operation>& ops, const Mat3& R) {
  for (Operation& op : ops) {
    op.matrix = R * op.matrix * R.transpose();
  }
}

void tag_planes(std::vector<Operation>& ops, const std::string& name, int n) {
  const bool even_axial = n >= 2 && (n % 2 == 0) && !name.empty() &&
                          (name[0] == 'C' || name[0] == 'D') &&
                          (name.back() == 'v' || name.back() == 'h' || name.back() == 'd');
  for (Operation& op : ops) {
    if (op.kind != OperationKind::Reflection || op.axis.norm() < 1e-12) {
      continue;
    }
    const Vec3 nrm = op.axis.normalized();
    const double dz = std::abs(nrm.dot(Vec3::UnitZ()));
    if (dz > 0.9) {
      op.plane = PlaneKind::Horizontal;
    } else if (dz < 0.2) {
      op.plane =
          (!even_axial || std::abs(nrm.dot(Vec3::UnitX())) >= std::abs(nrm.dot(Vec3::UnitY())))
              ? PlaneKind::Vertical
              : PlaneKind::Dihedral;
    } else {
      op.plane = PlaneKind::Dihedral;
    }
  }
}

Mat3 basis_from_z_and_secondary(Vec3 z, Vec3 secondary) {
  z.normalize();
  Vec3 y = z.cross(secondary);
  if (y.norm() < 1e-8) {
    secondary = (std::abs(z.x()) < 0.9) ? Vec3::UnitX() : Vec3::UnitY();
    y = z.cross(secondary);
  }
  y.normalize();
  Vec3 x = y.cross(z);
  Mat3 R;
  R.row(0) = x.transpose();
  R.row(1) = y.transpose();
  R.row(2) = z.transpose();
  if (R.determinant() < 0.0) {
    R.row(1) *= -1.0;
  }
  return R;
}

bool is_c2(const Operation& op) {
  return op.kind == OperationKind::ProperRotation &&
         (op.order == 2 || (op.order == 4 && op.power == 2));
}

void cotton_orient(std::vector<Operation>& ops, const std::string& name) {
  if (!name.empty() && (name[0] == 'C' || name[0] == 'D' || name[0] == 'S') && name != "Cs" &&
      name != "Ci" && name.find("inf") == std::string::npos) {
    int best_n = 1;
    Vec3 primary = Vec3::UnitZ();
    for (const Operation& op : ops) {
      if ((op.kind == OperationKind::ProperRotation ||
           op.kind == OperationKind::ImproperRotation) &&
          op.order > best_n) {
        best_n = op.order;
        primary = op.axis.normalized();
      }
    }
    Vec3 secondary = Vec3::UnitX();
    for (const Operation& op : ops) {
      if (op.kind == OperationKind::ProperRotation && op.order == 2 &&
          std::abs(op.axis.normalized().dot(primary)) < 0.2) {
        secondary = op.axis.normalized();
        break;
      }
    }
    if (std::abs(primary.z()) < 0.99 || std::abs(secondary.x()) < 0.99) {
      conjugate(ops, basis_from_z_and_secondary(primary, secondary));
      retag_operations(ops);
    }
    double ang = 0.0;
    bool found = false;
    for (const Operation& op : ops) {
      if (op.kind == OperationKind::Reflection &&
          std::abs(op.axis.normalized().dot(Vec3::UnitZ())) < 0.2) {
        const Vec3 n = op.axis.normalized();
        ang = std::atan2(n.x(), -n.y());
        found = true;
        break;
      }
    }
    if (!found) {
      for (const Operation& op : ops) {
        if (op.kind == OperationKind::ProperRotation && op.order == 2 &&
            std::abs(op.axis.normalized().dot(Vec3::UnitZ())) < 0.2) {
          const Vec3 a = op.axis.normalized();
          ang = std::atan2(a.y(), a.x());
          found = true;
          break;
        }
      }
    }
    if (found && std::abs(ang) > 1e-10) {
      conjugate(ops, rotation_matrix(Vec3::UnitZ(), -ang));
      retag_operations(ops);
    }
    return;
  }

  if (name == "Cs") {
    for (const Operation& op : ops) {
      if (op.kind == OperationKind::Reflection) {
        conjugate(ops, basis_from_z_and_secondary(op.axis.normalized(), Vec3::UnitX()));
        retag_operations(ops);
        break;
      }
    }
    return;
  }

  if (name == "T" || name == "Td" || name == "Th" || name == "O" || name == "Oh" || name == "I" ||
      name == "Ih") {
    const bool tetra = (name == "T" || name == "Td" || name == "Th");
    int target_n = 2;
    if (name == "O" || name == "Oh") {
      target_n = 4;
    } else if (name == "I" || name == "Ih") {
      target_n = 5;
    }
    Vec3 primary = Vec3::UnitZ();
    double best_abs = -1.0;
    for (const Operation& op : ops) {
      if (op.kind != OperationKind::ProperRotation) {
        continue;
      }
      const bool match =
          tetra ? is_c2(op) : (op.order == target_n && (op.power == 1 || op.power == target_n - 1));
      if (!match) {
        continue;
      }
      Vec3 a = op.axis.normalized();
      if (a.z() < 0.0) {
        a = -a;
      }
      if (std::abs(a.z()) > best_abs) {
        best_abs = std::abs(a.z());
        primary = a;
      }
    }
    Vec3 secondary = Vec3::UnitX();
    for (const Operation& op : ops) {
      if (is_c2(op) && std::abs(op.axis.normalized().dot(primary)) < 0.2) {
        secondary = op.axis.normalized();
        break;
      }
    }
    if (best_abs >= 0.0) {
      conjugate(ops, basis_from_z_and_secondary(primary, secondary));
      retag_operations(ops);
    }
    for (const Operation& op : ops) {
      if (!is_c2(op) || std::abs(op.axis.normalized().dot(Vec3::UnitZ())) > 0.2) {
        continue;
      }
      const double ang = std::atan2(op.axis.normalized().y(), op.axis.normalized().x());
      if (std::abs(ang) > 1e-10) {
        conjugate(ops, rotation_matrix(Vec3::UnitZ(), -ang));
        retag_operations(ops);
      }
      break;
    }
  }
}

std::vector<Operation> generators_for(const std::string& name) {
  const Vec3 z = Vec3::UnitZ();
  const Vec3 x = Vec3::UnitX();
  const Mat3 inv = -Mat3::Identity();
  const Mat3 sh = reflection_matrix(z);
  const Mat3 sv = reflection_matrix(Vec3::UnitY());
  const int n = parse_n(name);
  auto rot = [](const Vec3& ax, double ang) { return rotation_matrix(ax, ang); };

  if (name == "C1") {
    return {};
  }
  if (name == "Cs") {
    return {from_matrix(sh)};
  }
  if (name == "Ci") {
    return {from_matrix(inv)};
  }
  if (name == "T") {
    return {from_matrix(rot(Vec3(1, 1, 1), 2.0 * kPi / 3.0)), from_matrix(rot(z, kPi))};
  }
  if (name == "Th") {
    return {from_matrix(rot(Vec3(1, 1, 1), 2.0 * kPi / 3.0)),
            from_matrix(rot(z, kPi)),
            from_matrix(inv)};
  }
  if (name == "Td") {
    return {from_matrix(rot(Vec3(1, 1, 1), 2.0 * kPi / 3.0)), from_matrix(sh * rot(z, 0.5 * kPi))};
  }
  if (name == "O") {
    return {from_matrix(rot(Vec3(1, 1, 1), 2.0 * kPi / 3.0)), from_matrix(rot(z, 0.5 * kPi))};
  }
  if (name == "Oh") {
    return {from_matrix(rot(Vec3(1, 1, 1), 2.0 * kPi / 3.0)),
            from_matrix(rot(z, 0.5 * kPi)),
            from_matrix(inv)};
  }
  if (name == "I" || name == "Ih") {
    const double phi = 0.5 * (1.0 + std::sqrt(5.0));
    std::vector<Operation> gens{from_matrix(rot(Vec3(1, 1, 1), 2.0 * kPi / 3.0)),
                                from_matrix(rot(Vec3(0, 1, phi), 2.0 * kPi / 5.0))};
    if (name == "Ih") {
      gens.push_back(from_matrix(inv));
    }
    return gens;
  }
  if (!name.empty() && name[0] == 'S' && n >= 2) {
    return {from_matrix(sh * rot(z, 2.0 * kPi / n))};
  }
  if (!name.empty() && name[0] == 'C' && n >= 2) {
    const Mat3 cn = rot(z, 2.0 * kPi / n);
    if (name.back() == 'v') {
      return {from_matrix(cn), from_matrix(sv)};
    }
    if (name.back() == 'h') {
      return {from_matrix(cn), from_matrix(sh)};
    }
    return {from_matrix(cn)};
  }
  if (!name.empty() && name[0] == 'D' && n >= 2) {
    const Mat3 cn = rot(z, 2.0 * kPi / n);
    const Mat3 c2x = rot(x, kPi);
    if (name.back() == 'h') {
      return {from_matrix(cn), from_matrix(c2x), from_matrix(sh)};
    }
    if (name.back() == 'd') {
      return {from_matrix(sh * rot(z, kPi / n)), from_matrix(c2x)};
    }
    return {from_matrix(cn), from_matrix(c2x)};
  }
  throw Error("unsupported or infinite point group: " + name);
}

} // namespace

std::vector<Operation> realize_point_group(const std::string& name) {
  if (name == "Cinfv" || name == "Dinfh" || name == "Kh") {
    throw Error("asymmetric unit is undefined for infinite point groups");
  }
  const int n = parse_n(name);
  const int expect = expected_group_order(name, n);
  if (expect < 0) {
    throw Error("unknown point group: " + name);
  }
  std::vector<Operation> ops = generators_for(name);
  close_group_abstract(ops, kMatTol, kMaxOrder);
  cotton_orient(ops, name);
  if (static_cast<int>(ops.size()) != expect) {
    throw Error("failed to realize point group " + name + " (order " + std::to_string(ops.size()) +
                " != " + std::to_string(expect) + ")");
  }
  tag_planes(ops, name, n);
  return ops;
}

std::string name_of_operations(const std::vector<Operation>& ops) {
  bool has_inversion = false;
  bool has_sigma_h = false;
  bool has_sigma_v = false;
  int highest_cn = 1;
  int n_c3 = 0;
  int n_c4 = 0;
  int n_c5 = 0;
  std::vector<Vec3> c2_axes;
  Vec3 primary = Vec3::UnitZ();

  for (const Operation& op : ops) {
    if (op.kind == OperationKind::Inversion) {
      has_inversion = true;
    }
    if (op.kind == OperationKind::ProperRotation && op.order > highest_cn) {
      highest_cn = op.order;
      primary = op.axis.normalized();
    }
  }
  for (const Operation& op : ops) {
    if (op.kind == OperationKind::ProperRotation && op.power == 1) {
      if (op.order == 3) {
        ++n_c3;
      }
      if (op.order == 4) {
        ++n_c4;
      }
      if (op.order == 5) {
        ++n_c5;
      }
      if (op.order == 2 || (op.order >= 2 && op.power == op.order / 2)) {
        Vec3 a = op.axis.normalized();
        canonicalize_direction(a);
        bool seen = false;
        for (const Vec3& e : c2_axes) {
          if (parallel_axis(e, a, 1e-2)) {
            seen = true;
            break;
          }
        }
        if (!seen) {
          c2_axes.push_back(a);
        }
      }
    }
    if (op.kind == OperationKind::Reflection) {
      if (std::abs(op.axis.normalized().dot(primary)) > 0.9) {
        has_sigma_h = true;
      } else {
        has_sigma_v = true;
      }
    }
  }
  int n_perp_c2 = 0;
  for (const Vec3& ax : c2_axes) {
    if (std::abs(ax.dot(primary)) < 0.2) {
      ++n_perp_c2;
    }
  }
  return classify_point_group(ops,
                              false,
                              has_inversion,
                              has_sigma_h,
                              has_sigma_v,
                              highest_cn,
                              n_perp_c2,
                              n_c3,
                              n_c4,
                              n_c5,
                              n_c3 >= 4);
}

} // namespace molsymm::detail

namespace molsymm {
namespace {

constexpr double kPosTol = 1e-3;

std::vector<Operation> stabilizer(const std::vector<Operation>& ops, const Vec3& p, double tol) {
  std::vector<Operation> stab;
  for (const Operation& op : ops) {
    if ((op.matrix * p - p).norm() <= tol) {
      stab.push_back(op);
    }
  }
  return stab;
}

Mat3 site_projector_of(const std::vector<Operation>& stab) {
  Mat3 P = Mat3::Zero();
  for (const Operation& op : stab) {
    P += op.matrix;
  }
  return P / static_cast<double>(stab.size());
}

int stratum_dof(const Mat3& projector) {
  Eigen::SelfAdjointEigenSolver<Mat3> es(0.5 * (projector + projector.transpose()));
  int dof = 0;
  for (int k = 0; k < 3; ++k) {
    if (es.eigenvalues()(k) > 0.5) {
      ++dof;
    }
  }
  return dof;
}

Vec3 eigenspace_axis(const Mat3& P, bool plus) {
  Eigen::SelfAdjointEigenSolver<Mat3> es(0.5 * (P + P.transpose()));
  for (int k = 0; k < 3; ++k) {
    if (plus ? es.eigenvalues()(k) > 0.5 : es.eigenvalues()(k) <= 0.5) {
      Vec3 v = es.eigenvectors().col(k).normalized();
      detail::canonicalize_direction(v);
      return v;
    }
  }
  return Vec3::UnitZ();
}

bool axes_equivalent(const Vec3& a, const Vec3& b, const std::vector<Operation>& ops, double tol) {
  const Vec3 ua = a.normalized();
  const Vec3 ub = b.normalized();
  return std::any_of(ops.begin(), ops.end(), [&](const Operation& op) {
    const Vec3 img = op.matrix * ua;
    return (img - ub).norm() <= tol || (img + ub).norm() <= tol;
  });
}

std::vector<Vec3> expand_orbit(const std::vector<Operation>& ops, const Vec3& p, double tol) {
  std::vector<Vec3> orbit;
  for (const Operation& op : ops) {
    const Vec3 q = op.matrix * p;
    if (std::none_of(
            orbit.begin(), orbit.end(), [&](const Vec3& r) { return (r - q).norm() <= tol; })) {
      orbit.push_back(q);
    }
  }
  return orbit;
}

std::vector<Vec3>
chamber_of(const std::vector<Operation>& ops, const Mat3& projector, const Vec3& test, double tol) {
  const Vec3 t = projector * test;
  if (t.norm() < tol) {
    return {};
  }
  std::vector<Vec3> normals;
  for (const Operation& op : ops) {
    const Vec3 img = projector * (op.matrix * t);
    Vec3 n = projector * (img - t);
    if (n.norm() < tol) {
      continue;
    }
    n.normalize();
    if (n.dot(t) < 0.0) {
      n = -n;
    }
    if (n.dot(t) <= tol) {
      continue;
    }
    if (std::none_of(
            normals.begin(), normals.end(), [&](const Vec3& e) { return (e - n).norm() < 1e-6; })) {
      normals.push_back(n);
    }
  }
  return normals;
}

struct SiteSeed {
  Vec3 test;
  Mat3 projector;
  int multiplicity = 1;
  int dof = 3;
  std::string site_symmetry;
  Vec3 sort_key = Vec3::Zero();
};

} // namespace

std::vector<WyckoffSite> wyckoff_positions(const std::string& point_group) {
  const std::vector<Operation> ops = detail::realize_point_group(point_group);
  const int order = static_cast<int>(ops.size());

  std::vector<Vec3> axis_dirs;
  for (const Operation& op : ops) {
    if (op.kind != OperationKind::ProperRotation || op.order < 2 || op.axis.norm() < 1e-12) {
      continue;
    }
    Vec3 a = op.axis.normalized();
    detail::canonicalize_direction(a);
    if (std::none_of(axis_dirs.begin(), axis_dirs.end(), [&](const Vec3& e) {
          return axes_equivalent(e, a, ops, 1e-6);
        })) {
      axis_dirs.push_back(a);
    }
  }

  std::vector<Vec3> plane_normals;
  for (const Operation& op : ops) {
    if (op.kind != OperationKind::Reflection || op.axis.norm() < 1e-12) {
      continue;
    }
    Vec3 n = op.axis.normalized();
    detail::canonicalize_direction(n);
    if (std::none_of(plane_normals.begin(), plane_normals.end(), [&](const Vec3& e) {
          return axes_equivalent(e, n, ops, 1e-6);
        })) {
      plane_normals.push_back(n);
    }
  }

  std::vector<SiteSeed> seeds;

  auto add_seed = [&](const Vec3& raw, bool require_nonzero) {
    if (require_nonzero && raw.norm() < kPosTol) {
      return;
    }
    const auto stab = stabilizer(ops, raw, kPosTol);
    if (stab.empty()) {
      return;
    }
    const Mat3 P = site_projector_of(stab);
    const Vec3 snapped = P * raw;
    const int dof = stratum_dof(P);
    if (require_nonzero && snapped.norm() < kPosTol) {
      return;
    }
    const int mult = order / static_cast<int>(stab.size());
    if (mult * static_cast<int>(stab.size()) != order) {
      return;
    }
    for (const SiteSeed& s : seeds) {
      if (s.multiplicity != mult || s.dof != dof) {
        continue;
      }
      if (dof == 0 || dof == 3) {
        return;
      }
      if (dof == 1 && axes_equivalent(eigenspace_axis(P, true), s.sort_key, ops, 1e-5)) {
        return;
      }
      if (dof == 2 && axes_equivalent(eigenspace_axis(P, false), s.sort_key, ops, 1e-5)) {
        return;
      }
    }

    SiteSeed seed;
    seed.test = (snapped.norm() >= kPosTol) ? snapped : raw;
    seed.projector = P;
    seed.multiplicity = mult;
    seed.dof = dof;
    seed.site_symmetry = detail::name_of_operations(stab);
    if (dof == 0) {
      seed.sort_key = Vec3::Zero();
    } else if (dof == 1) {
      seed.sort_key = eigenspace_axis(P, true);
      if (seed.test.norm() < kPosTol) {
        seed.test = seed.sort_key;
      }
    } else if (dof == 2) {
      seed.sort_key = eigenspace_axis(P, false);
    } else {
      seed.sort_key = Vec3(1, 2, 3);
    }
    seeds.push_back(seed);
  };

  add_seed(Vec3::Zero(), false);
  const int origin_stab = (!seeds.empty() && seeds[0].dof == 0) ? order / seeds[0].multiplicity : 0;

  for (const Vec3& a : axis_dirs) {
    add_seed(a, true);
  }

  if (!seeds.empty() && seeds[0].dof == 0) {
    int max_other = 1;
    for (std::size_t i = 1; i < seeds.size(); ++i) {
      max_other = std::max(max_other, order / seeds[i].multiplicity);
    }
    if (origin_stab <= max_other) {
      seeds.erase(seeds.begin());
    }
  }

  for (const Vec3& n : plane_normals) {
    Vec3 u = n.cross(Vec3::UnitZ());
    if (u.norm() < 1e-6) {
      u = n.cross(Vec3::UnitX());
    }
    u.normalize();
    const Vec3 v = n.cross(u).normalized();
    Vec3 p = u + 0.37 * v;
    p -= p.dot(n) * n;
    for (int attempt = 0; attempt < 8; ++attempt) {
      bool on_axis = false;
      for (const Vec3& a : axis_dirs) {
        if (a.cross(p).norm() < 1e-3 * p.norm()) {
          on_axis = true;
          break;
        }
      }
      if (!on_axis && p.norm() > kPosTol) {
        break;
      }
      p = u + (0.37 + 0.11 * attempt) * v + 0.05 * attempt * u;
      p -= p.dot(n) * n;
    }
    add_seed(p.normalized(), true);
  }

  {
    Vec3 g(0.91, 0.37, 0.53);
    for (int attempt = 0; attempt < 16; ++attempt) {
      bool bad = false;
      for (const Vec3& a : axis_dirs) {
        if (a.cross(g).norm() < 0.05 * g.norm()) {
          bad = true;
        }
      }
      for (const Vec3& n : plane_normals) {
        if (std::abs(n.dot(g.normalized())) < 0.05) {
          bad = true;
        }
      }
      if (!bad) {
        break;
      }
      g = Vec3(0.91 + 0.07 * attempt, 0.37 - 0.05 * attempt, 0.53 + 0.03 * attempt);
    }
    add_seed(g, true);
  }

  std::sort(seeds.begin(), seeds.end(), [](const SiteSeed& a, const SiteSeed& b) {
    if (a.multiplicity != b.multiplicity) {
      return a.multiplicity < b.multiplicity;
    }
    if (a.dof != b.dof) {
      return a.dof < b.dof;
    }
    for (int k = 0; k < 3; ++k) {
      if (std::abs(a.sort_key(k) - b.sort_key(k)) > 1e-8) {
        return a.sort_key(k) < b.sort_key(k);
      }
    }
    return false;
  });

  std::vector<WyckoffSite> sites;
  sites.reserve(seeds.size());
  std::map<int, int> letter_count;
  for (const SiteSeed& seed : seeds) {
    WyckoffSite site;
    const int letter_i = letter_count[seed.multiplicity]++;
    site.label = std::to_string(seed.multiplicity) + static_cast<char>('a' + letter_i);
    site.site_symmetry = seed.site_symmetry;
    site.multiplicity = seed.multiplicity;
    site.degrees_of_freedom = seed.dof;
    site.projector = seed.projector;
    site.chamber = chamber_of(ops, seed.projector, seed.test, kPosTol);
    sites.push_back(std::move(site));
  }
  return sites;
}

namespace {

const WyckoffSite* find_site(const std::vector<WyckoffSite>& sites, const std::string& label) {
  for (const WyckoffSite& s : sites) {
    if (s.label == label) {
      return &s;
    }
  }
  return nullptr;
}

bool is_origin_site(const WyckoffSite& s) {
  return s.degrees_of_freedom == 0;
}

void partitions_into(int remaining, // NOLINT(misc-no-recursion)
                     const std::vector<int>& mults,
                     int min_idx,
                     std::vector<int>& current,
                     std::vector<std::vector<int>>& out) {
  if (remaining == 0) {
    out.push_back(current);
    return;
  }
  for (int i = min_idx; i < static_cast<int>(mults.size()); ++i) {
    const int m = mults[static_cast<std::size_t>(i)];
    if (m > remaining) {
      continue;
    }
    current.push_back(m);
    partitions_into(remaining - m, mults, i, current, out);
    current.pop_back();
  }
}

struct OrbitChoice {
  std::string label;
  int multiplicity = 1;
  int dof = 3;
  int atomic_number = 0;
};

void enumerate_element_assignments(const std::vector<WyckoffSite>& sites,
                                   int z,
                                   int count,
                                   bool origin_available,
                                   std::vector<std::vector<OrbitChoice>>& out) {
  std::vector<int> mults;
  for (const WyckoffSite& s : sites) {
    if (is_origin_site(s) && !origin_available) {
      continue;
    }
    if (std::find(mults.begin(), mults.end(), s.multiplicity) == mults.end()) {
      mults.push_back(s.multiplicity);
    }
  }
  std::sort(mults.begin(), mults.end());
  std::vector<std::vector<int>> parts;
  std::vector<int> part_buf;
  partitions_into(count, mults, 0, part_buf, parts);

  for (const std::vector<int>& part : parts) {
    std::map<int, int> need;
    for (int m : part) {
      ++need[m];
    }
    std::map<int, std::vector<const WyckoffSite*>> options;
    for (const WyckoffSite& s : sites) {
      if (is_origin_site(s) && !origin_available) {
        continue;
      }
      if (need.contains(s.multiplicity)) {
        options[s.multiplicity].push_back(&s);
      }
    }
    bool ok = true;
    for (const auto& [m, k] : need) {
      const auto& opts = options[m];
      if (opts.empty()) {
        ok = false;
        break;
      }
      const bool only_origin = std::all_of(
          opts.begin(), opts.end(), [](const WyckoffSite* s) { return is_origin_site(*s); });
      if (only_origin && k > 1) {
        ok = false;
        break;
      }
    }
    if (!ok) {
      continue;
    }

    std::vector<int> mult_keys;
    mult_keys.reserve(need.size());
    for (const auto& [m, _] : need) {
      mult_keys.push_back(m);
    }
    std::vector<std::vector<std::vector<const WyckoffSite*>>> per_mult_combos;
    for (int m : mult_keys) {
      const int k = need[m];
      const auto& opts = options[m];
      std::vector<std::vector<const WyckoffSite*>> combos;
      std::vector<const WyckoffSite*> chosen;
      const int nopts = static_cast<int>(opts.size());
      std::function<void(int, int, bool)> rec = [&](int start, int left, bool origin_used) {
        if (left == 0) {
          combos.push_back(chosen);
          return;
        }
        for (int i = start; i < nopts; ++i) {
          const WyckoffSite* s = opts[static_cast<std::size_t>(i)];
          const bool is_orig = is_origin_site(*s);
          if (is_orig && origin_used) {
            continue;
          }
          chosen.push_back(s);
          rec(is_orig ? i + 1 : i, left - 1, origin_used || is_orig);
          chosen.pop_back();
        }
      };
      rec(0, k, false);
      per_mult_combos.push_back(std::move(combos));
    }

    std::vector<std::size_t> cursor(per_mult_combos.size(), 0);
    while (true) {
      std::vector<OrbitChoice> assignment;
      bool origin_used = false;
      bool valid = true;
      for (std::size_t mi = 0; mi < per_mult_combos.size(); ++mi) {
        for (const WyckoffSite* s : per_mult_combos[mi][cursor[mi]]) {
          if (is_origin_site(*s)) {
            if (origin_used || !origin_available) {
              valid = false;
              break;
            }
            origin_used = true;
          }
          assignment.push_back({s->label, s->multiplicity, s->degrees_of_freedom, z});
        }
        if (!valid) {
          break;
        }
      }
      if (valid) {
        std::sort(
            assignment.begin(), assignment.end(), [](const OrbitChoice& a, const OrbitChoice& b) {
              if (a.multiplicity != b.multiplicity) {
                return a.multiplicity < b.multiplicity;
              }
              return a.label < b.label;
            });
        out.push_back(std::move(assignment));
      }
      int pos = static_cast<int>(cursor.size()) - 1;
      while (pos >= 0) {
        ++cursor[static_cast<std::size_t>(pos)];
        if (cursor[static_cast<std::size_t>(pos)] <
            per_mult_combos[static_cast<std::size_t>(pos)].size()) {
          break;
        }
        cursor[static_cast<std::size_t>(pos)] = 0;
        --pos;
      }
      if (pos < 0) {
        break;
      }
    }
  }
}

WyckoffPlan plan_from_orbits(const std::string& point_group,
                             int order,
                             const std::vector<OrbitChoice>& orbits) {
  WyckoffPlan plan;
  plan.point_group = point_group;
  plan.order = order;
  plan.orbits.reserve(orbits.size());
  for (const OrbitChoice& oc : orbits) {
    plan.orbits.push_back({oc.atomic_number, oc.label, oc.multiplicity, oc.dof});
  }
  return plan;
}

int plan_dof(const WyckoffPlan& p) {
  int d = 0;
  for (const auto& o : p.orbits) {
    d += o.degrees_of_freedom;
  }
  return d;
}

} // namespace

std::vector<WyckoffPlan> wyckoff_plans(const std::string& point_group,
                                       const std::vector<int>& atomic_numbers) {
  const auto sites = wyckoff_positions(point_group);
  const int order = detail::expected_group_order(point_group, detail::parse_n(point_group));
  if (atomic_numbers.empty()) {
    return {};
  }

  std::map<int, int> counts;
  for (int z : atomic_numbers) {
    ++counts[z];
  }
  std::vector<int> zs;
  zs.reserve(counts.size());
  for (const auto& [z, _] : counts) {
    zs.push_back(z);
  }

  struct Partial {
    std::vector<OrbitChoice> orbits;
    bool origin_used = false;
  };
  std::vector<Partial> acc{{{}, false}};

  for (int z : zs) {
    std::vector<Partial> next;
    for (const Partial& part : acc) {
      std::vector<std::vector<OrbitChoice>> assigns;
      enumerate_element_assignments(sites, z, counts[z], !part.origin_used, assigns);
      for (auto& asg : assigns) {
        Partial p = part;
        for (const OrbitChoice& oc : asg) {
          const WyckoffSite* s = find_site(sites, oc.label);
          if (s != nullptr && is_origin_site(*s)) {
            p.origin_used = true;
          }
        }
        p.orbits.insert(p.orbits.end(), asg.begin(), asg.end());
        next.push_back(std::move(p));
      }
    }
    acc = std::move(next);
    if (acc.empty()) {
      return {};
    }
  }

  std::vector<WyckoffPlan> plans;
  std::vector<std::string> seen;
  for (Partial& p : acc) {
    std::sort(p.orbits.begin(), p.orbits.end(), [](const OrbitChoice& a, const OrbitChoice& b) {
      if (a.atomic_number != b.atomic_number) {
        return a.atomic_number < b.atomic_number;
      }
      return a.label < b.label;
    });
    std::ostringstream key;
    for (const OrbitChoice& oc : p.orbits) {
      key << oc.atomic_number << ':' << oc.label << ';';
    }
    const std::string k = key.str();
    if (std::find(seen.begin(), seen.end(), k) != seen.end()) {
      continue;
    }
    seen.push_back(k);
    plans.push_back(plan_from_orbits(point_group, order, p.orbits));
  }

  std::sort(plans.begin(), plans.end(), [](const WyckoffPlan& a, const WyckoffPlan& b) {
    if (a.orbits.size() != b.orbits.size()) {
      return a.orbits.size() < b.orbits.size();
    }
    return plan_dof(a) < plan_dof(b);
  });
  return plans;
}

std::vector<WyckoffPlan> wyckoff_plans(const std::string& point_group, int n_sites) {
  if (n_sites < 0) {
    throw Error("n_sites must be non-negative");
  }
  if (n_sites == 0) {
    return {};
  }
  return wyckoff_plans(point_group, std::vector<int>(static_cast<std::size_t>(n_sites), 0));
}

WyckoffPlan wyckoff_plan(const std::string& point_group,
                         const std::vector<std::string>& site_labels,
                         const std::vector<int>& atomic_numbers) {
  if (!atomic_numbers.empty() && atomic_numbers.size() != site_labels.size()) {
    throw Error("atomic_numbers size must match site_labels");
  }
  const auto sites = wyckoff_positions(point_group);
  const int order = detail::expected_group_order(point_group, detail::parse_n(point_group));
  bool origin_used = false;
  std::vector<OrbitChoice> orbits;
  for (std::size_t i = 0; i < site_labels.size(); ++i) {
    const WyckoffSite* s = find_site(sites, site_labels[i]);
    if (s == nullptr) {
      throw Error("unknown Wyckoff site label: " + site_labels[i]);
    }
    if (is_origin_site(*s)) {
      if (origin_used) {
        throw Error("origin Wyckoff site used more than once");
      }
      origin_used = true;
    }
    orbits.push_back({s->label,
                      s->multiplicity,
                      s->degrees_of_freedom,
                      atomic_numbers.empty() ? 0 : atomic_numbers[i]});
  }
  return plan_from_orbits(point_group, order, orbits);
}

Analysis assemble_from_au(const std::string& point_group,
                          const std::vector<int>& representative_numbers,
                          const std::vector<Vec3>& au_positions,
                          const std::vector<std::string>& site_labels) {
  if (representative_numbers.size() != au_positions.size() ||
      au_positions.size() != site_labels.size()) {
    throw Error("assemble_from_au: numbers, positions, and site labels must match");
  }
  return wyckoff_plan(point_group, site_labels, representative_numbers).assemble(au_positions);
}

Analysis WyckoffPlan::assemble(const std::vector<Vec3>& au_positions) const {
  if (au_positions.size() != orbits.size()) {
    throw Error("assemble: expected one position per orbit");
  }
  const auto sites = wyckoff_positions(point_group);
  std::vector<Operation> ops = detail::realize_point_group(point_group);

  std::vector<int> zs;
  std::vector<Vec3> positions;

  for (std::size_t oi = 0; oi < orbits.size(); ++oi) {
    const WyckoffOrbit& orb = orbits[oi];
    const WyckoffSite* site = find_site(sites, orb.label);
    if (site == nullptr) {
      throw Error("assemble: unknown site " + orb.label);
    }
    const Vec3& raw = au_positions[oi];
    Vec3 p;
    if (site->degrees_of_freedom == 0) {
      if (raw.norm() > kPosTol) {
        throw Error("assemble: origin site requires a near-zero position");
      }
      p = Vec3::Zero();
    } else {
      if (raw.norm() < kPosTol) {
        throw Error("assemble: non-origin site requires a non-zero position");
      }
      const std::vector<Vec3> probe = expand_orbit(ops, raw, kPosTol);
      double best_res = 1e300;
      std::size_t best_k = 0;
      for (std::size_t k = 0; k < probe.size(); ++k) {
        const double res = (probe[k] - site->projector * probe[k]).norm();
        if (res < best_res) {
          best_res = res;
          best_k = k;
        }
      }
      if (best_res > kPosTol) {
        throw Error("assemble: position is off the " + orb.label + " stratum");
      }
      p = site->projector * probe[best_k];
    }

    std::vector<Vec3> orbit = expand_orbit(ops, p, kPosTol);
    if (static_cast<int>(orbit.size()) != orb.multiplicity) {
      throw Error("assemble: orbit size " + std::to_string(orbit.size()) + " != multiplicity " +
                  std::to_string(orb.multiplicity) + " for " + orb.label);
    }
    std::size_t best = 0;
    double best_d = (orbit[0] - raw).norm();
    for (std::size_t k = 1; k < orbit.size(); ++k) {
      const double d = (orbit[k] - raw).norm();
      if (d < best_d) {
        best_d = d;
        best = k;
      }
    }
    if (best != 0) {
      std::swap(orbit[0], orbit[best]);
    }

    for (const Vec3& q : orbit) {
      for (const Vec3& prev : positions) {
        if ((q - prev).norm() <= kPosTol) {
          throw Error("assemble: coincident atoms from distinct orbits");
        }
      }
    }
    for (const Vec3& q : orbit) {
      zs.push_back(orb.atomic_number);
      positions.push_back(q);
    }
  }

  Molecule mol;
  mol.atoms.resize(zs.size());
  for (std::size_t i = 0; i < zs.size(); ++i) {
    mol.atoms[i].atomic_number = zs[i];
    mol.atoms[i].position = positions[i];
  }
  mol.resolve_masses();

  const Vec3 com = center_of_mass(mol);
  const bool com_invariant = std::all_of(ops.begin(), ops.end(), [&](const Operation& op) {
    return (op.matrix * com - com).norm() <= kPosTol;
  });
  if (com_invariant) {
    for (Atom& a : mol.atoms) {
      a.position -= com;
    }
  }

  const auto by_z = detail::build_element_index(mol.atoms);
  for (Operation& op : ops) {
    op.image = detail::permutation(mol.atoms, by_z, op.matrix, kPosTol);
    if (op.image.empty()) {
      throw Error("assemble: group operation does not permute the assembled molecule");
    }
  }

  Analysis result;
  result.point_group = point_group;
  result.order = order;
  result.symmetry_number = detail::symmetry_number_for(point_group, detail::parse_n(point_group));
  result.molecule = mol;
  result.original = std::move(mol);
  result.operations = std::move(ops);
  result.matrix_tol = detail::matrix_tol_from_angle(1e-2);
  detail::assign_conjugacy_classes(result.operations, result.matrix_tol);
  result.equivalence_sets =
      detail::orbits_from_operations(result.operations, result.molecule.size());
  result.character_table =
      detail::build_character_table(point_group, detail::parse_n(point_group), result.operations);
  result.subgroups = detail::constructive_subgroups(point_group, detail::parse_n(point_group));
  return result;
}

} // namespace molsymm
