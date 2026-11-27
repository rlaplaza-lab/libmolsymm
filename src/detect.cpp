#include "molsymm/analysis.hpp"

#include "internal.hpp"
#include "molsymm/error.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <sstream>
#include <unordered_map>
#include <utility>
#include <vector>

namespace molsymm {
namespace {

void canonicalize_axis(Vec3& u) {
  if (u.z() < -1e-12 || (std::abs(u.z()) < 1e-12 && u.y() < -1e-12) ||
      (std::abs(u.z()) < 1e-12 && std::abs(u.y()) < 1e-12 && u.x() < 0)) {
    u = -u;
  }
}

struct AxisBucketKey {
  std::int32_t ix = 0;
  std::int32_t iy = 0;
  std::int32_t iz = 0;
  bool operator==(const AxisBucketKey& o) const { return ix == o.ix && iy == o.iy && iz == o.iz; }
};

struct AxisBucketHash {
  std::size_t operator()(const AxisBucketKey& k) const {
    return (static_cast<std::size_t>(static_cast<std::uint32_t>(k.ix)) * 73856093U) ^
           (static_cast<std::size_t>(static_cast<std::uint32_t>(k.iy)) * 19349663U) ^
           (static_cast<std::size_t>(static_cast<std::uint32_t>(k.iz)) * 83492791U);
  }
};

AxisBucketKey axis_bucket(const Vec3& u, double inv_bin) {
  return {static_cast<std::int32_t>(std::floor(u.x() * inv_bin + 0.5)),
          static_cast<std::int32_t>(std::floor(u.y() * inv_bin + 0.5)),
          static_cast<std::int32_t>(std::floor(u.z() * inv_bin + 0.5))};
}

struct AxisSet {
  std::vector<Vec3> axes;
  std::unordered_map<AxisBucketKey, std::vector<int>, AxisBucketHash> buckets;
  double inv_bin;

  explicit AxisSet(double angle_tol) : inv_bin(1.0 / std::max(angle_tol, 1e-4)) {}

  void add(const Vec3& a, double angle_tol) {
    if (a.norm() < 1e-10) {
      return;
    }
    Vec3 u = a.normalized();
    canonicalize_axis(u);
    const AxisBucketKey key = axis_bucket(u, inv_bin);
    for (int dx = -1; dx <= 1; ++dx) {
      for (int dy = -1; dy <= 1; ++dy) {
        for (int dz = -1; dz <= 1; ++dz) {
          auto it = buckets.find({key.ix + dx, key.iy + dy, key.iz + dz});
          if (it == buckets.end()) {
            continue;
          }
          for (int idx : it->second) {
            if (detail::parallel_axis(axes[static_cast<std::size_t>(idx)], u, angle_tol)) {
              return;
            }
          }
        }
      }
    }
    buckets[key].push_back(static_cast<int>(axes.size()));
    axes.push_back(u);
  }
};

void add_unique_axis(std::vector<Vec3>& axes, const Vec3& a, double angle_tol) {
  if (a.norm() < 1e-10) {
    return;
  }
  Vec3 u = a.normalized();
  canonicalize_axis(u);
  for (const Vec3& existing : axes) {
    if (detail::parallel_axis(existing, u, angle_tol)) {
      return;
    }
  }
  axes.push_back(u);
}

Operation make_cn(const Vec3& axis,
                  int n,
                  int power,
                  const std::vector<Atom>& atoms,
                  const detail::ElementIndex& by_z,
                  double tol) {
  Operation op;
  op.kind = OperationKind::ProperRotation;
  op.order = n;
  op.power = power;
  op.axis = axis.normalized();
  op.matrix = detail::rotation_matrix(op.axis, 2.0 * detail::kPi * power / n);
  op.image = detail::permutation(atoms, by_z, op.matrix, tol);
  return op;
}

Operation make_sigma(const Vec3& normal,
                     PlaneKind plane,
                     const std::vector<Atom>& atoms,
                     const detail::ElementIndex& by_z,
                     double tol) {
  Operation op;
  op.kind = OperationKind::Reflection;
  op.order = 1;
  op.axis = normal.normalized();
  op.plane = plane;
  op.matrix = detail::reflection_matrix(op.axis);
  op.image = detail::permutation(atoms, by_z, op.matrix, tol);
  return op;
}

Operation make_sn(const Vec3& axis,
                  int n,
                  int power,
                  const std::vector<Atom>& atoms,
                  const detail::ElementIndex& by_z,
                  double tol) {
  Operation op;
  op.kind = OperationKind::ImproperRotation;
  op.order = n;
  op.power = power;
  op.axis = axis.normalized();
  op.matrix = detail::improper_rotation_matrix(op.axis, 2.0 * detail::kPi * power / n);
  op.image = detail::permutation(atoms, by_z, op.matrix, tol);
  return op;
}

bool atoms_collinear(const Molecule& mol, double tol) {
  if (mol.size() <= 2) {
    return true;
  }
  int i0 = 0;
  double rmax = 0.0;
  for (int i = 0; i < mol.size(); ++i) {
    const double r = mol.atoms[static_cast<std::size_t>(i)].position.norm();
    if (r > rmax) {
      rmax = r;
      i0 = i;
    }
  }
  if (rmax < tol) {
    return true;
  }
  const Vec3 axis = mol.atoms[static_cast<std::size_t>(i0)].position.normalized();
  for (int i = 0; i < mol.size(); ++i) {
    const Vec3& p = mol.atoms[static_cast<std::size_t>(i)].position;
    if (p.norm() < tol) {
      continue;
    }
    const double cross = axis.cross(p).norm();
    if (cross > tol) {
      return false;
    }
  }
  return true;
}

struct AxisOrder {
  int n = 1;
  std::vector<int> image;
};

AxisOrder highest_order_on_axis(const std::vector<Atom>& atoms,
                                const detail::ElementIndex& by_z,
                                const Vec3& axis,
                                int max_n,
                                double tol) {
  // Scan upward: skip n if a proper divisor failed; after first hit g, only multiples of g.
  std::vector<char> failed(static_cast<std::size_t>(max_n) + 1, 0);
  int best_n = 1;
  std::vector<int> best_image;
  for (int n = 2; n <= max_n; ++n) {
    bool ruled_out = false;
    for (int d = 2; d < n; ++d) {
      if (n % d == 0 && failed[static_cast<std::size_t>(d)]) {
        ruled_out = true;
        break;
      }
    }
    if (ruled_out) {
      failed[static_cast<std::size_t>(n)] = 1;
      continue;
    }
    // Higher orders on this axis must be multiples of the current hit; do not mark
    // non-multiples failed (that would poison composites like C6 after C2 via C3).
    if (best_n > 1 && n % best_n != 0) {
      continue;
    }
    Mat3 R = detail::rotation_matrix(axis, 2.0 * detail::kPi / n);
    std::vector<int> image = detail::permutation(atoms, by_z, R, tol);
    if (!image.empty()) {
      best_n = n;
      best_image = std::move(image);
    } else {
      failed[static_cast<std::size_t>(n)] = 1;
    }
  }
  if (best_n < 2) {
    return {};
  }
  return {best_n, std::move(best_image)};
}

std::vector<int> compose_perm(const std::vector<int>& p, int times) {
  std::vector<int> out(p.size());
  for (std::size_t i = 0; i < p.size(); ++i) {
    out[i] = static_cast<int>(i);
  }
  for (int t = 0; t < times; ++t) {
    std::vector<int> next(out.size());
    for (std::size_t i = 0; i < out.size(); ++i) {
      next[i] = p[static_cast<std::size_t>(out[i])];
    }
    out.swap(next);
  }
  return out;
}

} // namespace

namespace detail {

std::string classify_point_group(const std::vector<Operation>& ops,
                                 bool linear,
                                 bool has_inversion,
                                 bool has_sigma_h,
                                 bool has_sigma_v,
                                 int highest_cn,
                                 int n_perp_c2,
                                 int n_c3,
                                 int n_c4,
                                 int n_c5,
                                 bool spherical_top) {
  if (linear) {
    return has_inversion || has_sigma_h ? "Dinfh" : "Cinfv";
  }

  auto name_n = [](char family, int n, const char* suffix) {
    std::ostringstream os;
    os << family << n << suffix;
    return os.str();
  };

  if (spherical_top && n_c3 >= 4) {
    if (n_c5 >= 6) {
      return has_inversion ? "Ih" : "I";
    }
    if (n_c4 >= 3) {
      return has_inversion ? "Oh" : "O";
    }
    bool has_s4 = false;
    bool has_sigma = false;
    for (const Operation& op : ops) {
      if (op.kind == OperationKind::ImproperRotation && op.order == 4) {
        has_s4 = true;
      }
      if (op.kind == OperationKind::Reflection) {
        has_sigma = true;
      }
    }
    if (has_s4 || (has_sigma && !has_inversion)) {
      return "Td";
    }
    if (has_inversion) {
      return "Th";
    }
    return "T";
  }

  if (highest_cn >= 2 && n_perp_c2 >= highest_cn) {
    if (has_sigma_h) {
      return name_n('D', highest_cn, "h");
    }
    // Dnd has σd and typically S_{2n}
    bool has_sigma_d = false;
    bool has_s2n = false;
    for (const Operation& op : ops) {
      if (op.kind == OperationKind::Reflection) {
        has_sigma_d = true;
      } else if (op.kind == OperationKind::ImproperRotation && op.order == 2 * highest_cn) {
        has_s2n = true;
      }
      if (has_sigma_d && has_s2n) {
        break;
      }
    }
    if (has_sigma_d || has_s2n) {
      return name_n('D', highest_cn, "d");
    }
    return name_n('D', highest_cn, "");
  }

  if (highest_cn >= 2) {
    if (has_sigma_h) {
      return name_n('C', highest_cn, "h");
    }
    if (has_sigma_v) {
      return name_n('C', highest_cn, "v");
    }
    for (const Operation& op : ops) {
      if (op.kind == OperationKind::ImproperRotation && op.order == 2 * highest_cn) {
        return name_n('S', 2 * highest_cn, "");
      }
    }
    return name_n('C', highest_cn, "");
  }

  if (has_inversion) {
    return "Ci";
  }
  for (const Operation& op : ops) {
    if (op.kind == OperationKind::Reflection) {
      return "Cs";
    }
  }
  return "C1";
}

int symmetry_number_for(const std::string& name, int n) {
  if (name == "C1" || name == "Cs" || name == "Ci" || name == "Cinfv" || name == "Kh") {
    return 1;
  }
  if (name == "Dinfh") {
    return 2;
  }
  if (name == "T" || name == "Td" || name == "Th") {
    return 12;
  }
  if (name == "O" || name == "Oh") {
    return 24;
  }
  if (name == "I" || name == "Ih") {
    return 60;
  }
  if (!name.empty() && name[0] == 'D') {
    return 2 * n;
  }
  if (!name.empty() && name[0] == 'S') {
    return std::max(n / 2, 1);
  }
  if (!name.empty() && name[0] == 'C') {
    return std::max(n, 1);
  }
  return 1;
}

static Vec3 sign_principal_axis(Vec3 v) {
  if (v.norm() < 1e-15) {
    return Vec3::UnitZ();
  }
  v.normalize();
  int best = 0;
  for (int i = 1; i < 3; ++i) {
    if (std::abs(v(i)) > std::abs(v(best)) + 1e-15) {
      best = i;
    }
  }
  if (v(best) < 0.0) {
    v = -v;
  }
  return v;
}

static void apply_orientation(Molecule& mol, std::vector<Operation>& ops, const Mat3& R) {
  for (Atom& a : mol.atoms) {
    a.position = R * a.position;
  }
  for (Operation& op : ops) {
    op.matrix = R * op.matrix * R.transpose();
    op.axis = R * op.axis;
    if (op.axis.norm() > 1e-12) {
      op.axis.normalize();
    }
  }
}

static Mat3 basis_from_z_and_secondary(Vec3 z, Vec3 secondary) {
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

static Vec3 pick_secondary_in_plane(const std::vector<Operation>& ops,
                                    const Vec3& primary,
                                    const Vec3& preferred_x) {
  Vec3 best = Vec3::Zero();
  double best_score = -1.0;
  auto consider = [&](Vec3 cand) {
    cand = cand - cand.dot(primary) * primary;
    if (cand.norm() < 1e-8) {
      return;
    }
    cand.normalize();
    if (cand.dot(preferred_x) < 0.0) {
      cand = -cand;
    }
    const double score = cand.dot(preferred_x);
    if (score > best_score) {
      best_score = score;
      best = cand;
    }
  };
  for (const Operation& op : ops) {
    if (op.kind == OperationKind::ProperRotation && op.order == 2 &&
        std::abs(op.axis.normalized().dot(primary)) < 0.2) {
      consider(op.axis.normalized());
    }
  }
  if (best_score >= 0.0) {
    return best;
  }
  for (const Operation& op : ops) {
    if (op.kind == OperationKind::Reflection && std::abs(op.axis.normalized().dot(primary)) < 0.2) {
      consider(op.axis.cross(primary));
    }
  }
  if (best_score >= 0.0) {
    return best;
  }
  Vec3 fallback = preferred_x - preferred_x.dot(primary) * primary;
  if (fallback.norm() < 1e-8) {
    fallback = (std::abs(primary.x()) < 0.9) ? Vec3::UnitX() : Vec3::UnitY();
    fallback = fallback - fallback.dot(primary) * primary;
  }
  return fallback.normalized();
}

static Vec3 pick_spherical_primary(const std::vector<Operation>& ops,
                                   const Vec3& input_z,
                                   const Vec3& input_x) {
  int best_n = 1;
  for (const Operation& op : ops) {
    if (op.kind == OperationKind::ProperRotation && op.order > best_n) {
      best_n = op.order;
    }
  }
  Vec3 best = Vec3::UnitZ();
  double best_abs = -1.0;
  double best_signed = 0.0;
  double best_abs_x = -1.0;
  for (const Operation& op : ops) {
    if (op.kind != OperationKind::ProperRotation || op.order != best_n) {
      continue;
    }
    Vec3 a = op.axis.normalized();
    double oz = a.dot(input_z);
    if (oz < 0.0) {
      a = -a;
      oz = -oz;
    }
    const double abs_z = std::abs(a.dot(input_z));
    const double abs_x = std::abs(a.dot(input_x));
    const bool better = abs_z > best_abs + 1e-12 ||
                        (std::abs(abs_z - best_abs) <= 1e-12 &&
                         (oz > best_signed + 1e-12 ||
                          (std::abs(oz - best_signed) <= 1e-12 && abs_x > best_abs_x + 1e-12)));
    if (better) {
      best_abs = abs_z;
      best_signed = oz;
      best_abs_x = abs_x;
      best = a;
    }
  }
  return best;
}

void align_to_frame(Molecule& mol,
                    std::vector<Operation>& ops,
                    Mat3& orientation,
                    const Tolerances& tol,
                    const Vec3& sevals,
                    bool spherical,
                    bool linear,
                    const std::string& point_group,
                    const Vec3& input_x_in_frame,
                    const Vec3& input_z_in_frame) {
  const double imean =
      (std::abs(sevals(0)) + std::abs(sevals(1)) + std::abs(sevals(2))) / 3.0 + 1e-12;
  const bool deg_xy = std::abs(sevals(0) - sevals(1)) / imean < tol.inertia;
  const bool deg_yz = std::abs(sevals(1) - sevals(2)) / imean < tol.inertia;

  Vec3 preferred_x = input_x_in_frame;
  if (preferred_x.norm() < 1e-8) {
    preferred_x = Vec3::UnitX();
  } else {
    preferred_x.normalize();
  }

  Vec3 z_axis = Vec3::UnitZ();
  Vec3 secondary = Vec3::UnitX();

  if (spherical) {
    z_axis = pick_spherical_primary(ops, input_z_in_frame, preferred_x);
    secondary = pick_secondary_in_plane(ops, z_axis, preferred_x);
  } else if (linear || (deg_yz && !deg_xy)) {
    z_axis = Vec3::UnitX();
    secondary = pick_secondary_in_plane(ops, z_axis, preferred_x);
    if (secondary.norm() < 1e-8) {
      secondary = Vec3::UnitY();
    }
  } else if (deg_xy && !deg_yz) {
    z_axis = Vec3::UnitZ();
    secondary = pick_secondary_in_plane(ops, z_axis, preferred_x);
  } else if (point_group == "D2" || point_group == "D2h" || point_group == "C1" ||
             point_group == "Ci") {
    z_axis = Vec3::UnitZ();
    secondary = Vec3::UnitX();
  } else if (point_group == "Cs") {
    for (const Operation& op : ops) {
      if (op.kind == OperationKind::Reflection) {
        z_axis = op.axis.normalized();
        if (z_axis.dot(input_z_in_frame) < 0.0) {
          z_axis = -z_axis;
        }
        break;
      }
    }
    secondary = preferred_x;
  } else {
    int best_n = 1;
    Vec3 primary = Vec3::UnitZ();
    for (const Operation& op : ops) {
      if (op.kind == OperationKind::ProperRotation && op.order > best_n) {
        best_n = op.order;
        primary = op.axis.normalized();
      }
    }
    if (best_n >= 2) {
      const std::array<Vec3, 3> axes{{Vec3::UnitX(), Vec3::UnitY(), Vec3::UnitZ()}};
      int best_axis = 0;
      double best_abs = -1.0;
      for (int i = 0; i < 3; ++i) {
        const double a = std::abs(primary.dot(axes[static_cast<std::size_t>(i)]));
        if (a > best_abs) {
          best_abs = a;
          best_axis = i;
        }
      }
      z_axis = axes[static_cast<std::size_t>(best_axis)];
      if (primary.dot(z_axis) < 0.0) {
        z_axis = -z_axis;
      }
      secondary = (best_axis == 0) ? Vec3::UnitY() : Vec3::UnitX();
    }
  }

  Mat3 R = basis_from_z_and_secondary(z_axis, secondary);
  orientation = R;
  apply_orientation(mol, ops, R);
}

} // namespace detail

Analysis analyze(const Molecule& molecule, const Tolerances& tol) {
  if (molecule.size() == 0) {
    throw Error("molecule has no atoms");
  }

  Analysis result;
  result.original = molecule;
  result.original.resolve_masses();

  if (molecule.size() == 1) {
    result.point_group = "Kh";
    result.order = 0;
    result.symmetry_number = 1;
    result.infinite = true;
    result.molecule = result.original;
    result.center_of_mass = result.molecule.atoms[0].position;
    result.molecule.atoms[0].position = Vec3::Zero();
    result.equivalence_sets = {EquivalenceSet{{0}}};
    result.subgroups = {"C1", "Cs", "Ci", "Kh"};
    result.character_table = detail::build_character_table("Kh", 0, {});
    return result;
  }

  Molecule mol = molecule;
  mol.resolve_masses();
  const Vec3 com = center_of_mass(mol);
  result.center_of_mass = com;
  for (Atom& a : mol.atoms) {
    a.position -= com;
  }

  Mat3 I = inertia_tensor(mol, Vec3::Zero());
  Eigen::SelfAdjointEigenSolver<Mat3> inertia_solver(I);
  Vec3 evals = inertia_solver.eigenvalues();
  Mat3 evecs = inertia_solver.eigenvectors();
  std::array<int, 3> order{{0, 1, 2}};
  std::sort(order.begin(), order.end(), [&](int a, int b) { return evals(a) < evals(b); });
  Mat3 sorted = Mat3::Zero();
  Vec3 sevals = Vec3::Zero();
  for (int i = 0; i < 3; ++i) {
    sevals(i) = evals(order[static_cast<std::size_t>(i)]);
    sorted.col(i) = detail::sign_principal_axis(evecs.col(order[static_cast<std::size_t>(i)]));
  }
  if (sorted.determinant() < 0) {
    sorted.col(1) *= -1;
  }

  Mat3 to_frame = sorted.transpose();
  const Vec3 input_x_in_frame = to_frame * Vec3::UnitX();
  const Vec3 input_z_in_frame = to_frame * Vec3::UnitZ();
  for (Atom& a : mol.atoms) {
    a.position = to_frame * a.position;
  }

  const double imean =
      (std::abs(sevals(0)) + std::abs(sevals(1)) + std::abs(sevals(2))) / 3.0 + 1e-12;
  const bool spherical = std::abs(sevals(0) - sevals(1)) / imean < tol.inertia &&
                         std::abs(sevals(1) - sevals(2)) / imean < tol.inertia;
  result.linear = atoms_collinear(mol, tol.primary);

  result.equivalence_sets = build_equivalence_sets(mol, tol.primary);

  const detail::ElementIndex by_z = detail::build_element_index(mol.atoms);

  std::vector<Operation> ops;
  {
    Operation e;
    e.kind = OperationKind::Identity;
    e.matrix = Mat3::Identity();
    e.image.resize(static_cast<std::size_t>(mol.size()));
    for (int i = 0; i < mol.size(); ++i) {
      e.image[static_cast<std::size_t>(i)] = i;
    }
    ops.push_back(e);
  }

  {
    const Mat3 inv_mat = detail::inversion_matrix();
    std::vector<int> image = detail::permutation(mol.atoms, by_z, inv_mat, tol.primary);
    if (!image.empty()) {
      Operation inv;
      inv.kind = OperationKind::Inversion;
      inv.matrix = inv_mat;
      inv.image = std::move(image);
      ops.push_back(inv);
    }
  }

  AxisSet axis_set(tol.angle);
  axis_set.add(Vec3::UnitX(), tol.angle);
  axis_set.add(Vec3::UnitY(), tol.angle);
  axis_set.add(Vec3::UnitZ(), tol.angle);

  for (const EquivalenceSet& eq_set : result.equivalence_sets) {
    for (std::size_t i = 0; i < eq_set.atoms.size(); ++i) {
      const Vec3& pi = mol.atoms[static_cast<std::size_t>(eq_set.atoms[i])].position;
      if (pi.norm() > tol.primary) {
        axis_set.add(pi, tol.angle);
      }
      for (std::size_t j = i + 1; j < eq_set.atoms.size(); ++j) {
        const Vec3& pj = mol.atoms[static_cast<std::size_t>(eq_set.atoms[j])].position;
        axis_set.add(pi + pj, tol.angle);
        axis_set.add(pi.cross(pj), tol.angle);
        axis_set.add(pj - pi, tol.angle);
      }
    }
  }
  const std::vector<Vec3>& axes = axis_set.axes;

  if (result.linear) {
    Operation c2 = make_cn(Vec3::UnitX(), 2, 1, mol.atoms, by_z, tol.primary);
    if (!c2.image.empty()) {
      ops.push_back(c2);
    }
    Operation sig = make_sigma(Vec3::UnitY(), PlaneKind::Vertical, mol.atoms, by_z, tol.primary);
    if (!sig.image.empty()) {
      ops.push_back(sig);
    }
    Operation sigh = make_sigma(Vec3::UnitX(), PlaneKind::Horizontal, mol.atoms, by_z, tol.primary);
    if (!sigh.image.empty()) {
      ops.push_back(sigh);
    }
  } else {
    for (const Vec3& axis : axes) {
      AxisOrder ao = highest_order_on_axis(mol.atoms, by_z, axis, tol.max_axis_order, tol.primary);
      if (ao.n >= 2) {
        const std::vector<int> gen_image = ao.image;
        for (int k = 1; k < ao.n; ++k) {
          Operation op;
          op.kind = OperationKind::ProperRotation;
          op.order = ao.n;
          op.power = k;
          op.axis = axis.normalized();
          op.matrix = detail::rotation_matrix(op.axis, 2.0 * detail::kPi * k / ao.n);
          if (k == 1) {
            op.image = gen_image;
          } else {
            op.image = compose_perm(gen_image, k);
          }
          ops.push_back(std::move(op));
        }
        for (int sn = 4; sn <= 2 * ao.n; sn += 2) {
          Operation op = make_sn(axis, sn, 1, mol.atoms, by_z, tol.primary);
          if (!op.image.empty()) {
            ops.push_back(op);
          }
        }
      }
      Operation sig = make_sigma(axis, PlaneKind::None, mol.atoms, by_z, tol.primary);
      if (!sig.image.empty()) {
        ops.push_back(sig);
      }
    }

    constexpr double kPairSigmaAngle = 1e-6;
    for (const EquivalenceSet& eq_set : result.equivalence_sets) {
      for (std::size_t i = 0; i < eq_set.atoms.size(); ++i) {
        for (std::size_t j = i + 1; j < eq_set.atoms.size(); ++j) {
          const Vec3& pi = mol.atoms[static_cast<std::size_t>(eq_set.atoms[i])].position;
          const Vec3& pj = mol.atoms[static_cast<std::size_t>(eq_set.atoms[j])].position;
          Vec3 nrm = (pi - pj);
          if (nrm.norm() < tol.primary) {
            continue;
          }
          bool already = false;
          for (const Vec3& ax : axes) {
            if (detail::parallel_axis(ax, nrm, kPairSigmaAngle)) {
              already = true;
              break;
            }
          }
          if (already) {
            continue;
          }
          Operation sig = make_sigma(nrm, PlaneKind::None, mol.atoms, by_z, tol.primary);
          if (!sig.image.empty()) {
            ops.push_back(sig);
          }
        }
      }
    }
  }

  const double mat_tol = detail::matrix_tol_from_angle(tol.angle);
  detail::dedupe_operations(ops, mat_tol);

  if (!result.linear) {
    detail::close_group(ops, mol.atoms, tol.primary, mat_tol, tol.max_group_order);
  }

  bool has_inversion = false;
  bool has_sigma_h = false;
  bool has_sigma_v = false;
  int highest_cn = 1;
  Vec3 primary_axis = Vec3::UnitZ();
  std::vector<Vec3> c2_axes;
  std::vector<Vec3> c3a;
  std::vector<Vec3> c4a;
  std::vector<Vec3> c5a;

  for (Operation& op : ops) {
    if (op.kind == OperationKind::Inversion) {
      has_inversion = true;
    }
    if (op.kind == OperationKind::ProperRotation) {
      if (op.order > highest_cn) {
        highest_cn = op.order;
        primary_axis = op.axis.normalized();
      }
      if (op.order == 2) {
        add_unique_axis(c2_axes, op.axis, tol.angle);
      } else if (op.order == 3) {
        add_unique_axis(c3a, op.axis, tol.angle);
      } else if (op.order == 4) {
        add_unique_axis(c4a, op.axis, tol.angle);
      } else if (op.order == 5) {
        add_unique_axis(c5a, op.axis, tol.angle);
      }
    }
  }
  const int n_c3 = static_cast<int>(c3a.size());
  const int n_c4 = static_cast<int>(c4a.size());
  const int n_c5 = static_cast<int>(c5a.size());

  int n_perp_c2 = 0;
  for (const Vec3& ax : c2_axes) {
    if (std::abs(ax.dot(primary_axis)) < 0.2) {
      ++n_perp_c2;
    }
  }

  for (Operation& op : ops) {
    if (op.kind != OperationKind::Reflection) {
      continue;
    }
    const double d = std::abs(op.axis.normalized().dot(primary_axis));
    if (d > 0.9) {
      op.plane = PlaneKind::Horizontal;
      has_sigma_h = true;
    } else if (d < 0.2) {
      op.plane = PlaneKind::Vertical;
      has_sigma_v = true;
    } else {
      op.plane = PlaneKind::Dihedral;
      has_sigma_v = true;
    }
  }

  result.point_group = detail::classify_point_group(ops,
                                                    result.linear,
                                                    has_inversion,
                                                    has_sigma_h,
                                                    has_sigma_v,
                                                    highest_cn,
                                                    n_perp_c2,
                                                    n_c3,
                                                    n_c4,
                                                    n_c5,
                                                    spherical);

  int n_group = detail::parse_n(result.point_group);
  if (n_group == 0) {
    n_group = highest_cn;
  }

  Mat3 extra = Mat3::Identity();
  detail::align_to_frame(mol,
                         ops,
                         extra,
                         tol,
                         sevals,
                         spherical,
                         result.linear,
                         result.point_group,
                         input_x_in_frame,
                         input_z_in_frame);

  {
    const bool even_axial = n_group >= 2 && (n_group % 2 == 0) && !result.point_group.empty() &&
                            (result.point_group[0] == 'C' || result.point_group[0] == 'D') &&
                            (result.point_group.back() == 'v' || result.point_group.back() == 'h');
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

  result.orientation = extra * to_frame;
  result.molecule = std::move(mol);
  result.operations = std::move(ops);
  result.infinite =
      result.point_group == "Cinfv" || result.point_group == "Dinfh" || result.point_group == "Kh";
  result.order = result.infinite ? 0 : static_cast<int>(result.operations.size());
  result.symmetry_number = detail::symmetry_number_for(result.point_group, n_group);

  result.matrix_tol = mat_tol;
  if (!result.infinite) {
    detail::assign_conjugacy_classes(result.operations, mat_tol);
  }

  if (!result.infinite && !result.operations.empty()) {
    result.equivalence_sets =
        detail::orbits_from_operations(result.operations, result.molecule.size());
  }

  result.character_table =
      detail::build_character_table(result.point_group, n_group, result.operations);
  result.subgroups = detail::constructive_subgroups(result.point_group, n_group);

  result.rmsd = result.symmetrize();

  return result;
}

} // namespace molsymm
