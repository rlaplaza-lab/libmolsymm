#include "molsymm/analysis.hpp"

#include "internal.hpp"
#include "molsymm/error.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <sstream>
#include <utility>

namespace molsymm {
namespace {

void add_unique_axis(std::vector<Vec3>& axes, const Vec3& a, double angle_tol) {
  if (a.norm() < 1e-10) {
    return;
  }
  Vec3 u = a.normalized();
  if (u.z() < -1e-12 || (std::abs(u.z()) < 1e-12 && u.y() < -1e-12) ||
      (std::abs(u.z()) < 1e-12 && std::abs(u.y()) < 1e-12 && u.x() < 0)) {
    u = -u;
  }
  for (const Vec3& existing : axes) {
    if (detail::parallel_axis(existing, u, angle_tol)) {
      return;
    }
  }
  axes.push_back(u);
}

Operation make_cn(const Vec3& axis, int n, int power, const std::vector<Atom>& atoms, double tol) {
  Operation op;
  op.kind = OperationKind::ProperRotation;
  op.order = n;
  op.power = power;
  op.axis = axis.normalized();
  op.matrix = detail::rotation_matrix(op.axis, 2.0 * detail::kPi * power / n);
  op.image = detail::permutation(atoms, op.matrix, tol);
  return op;
}

Operation
make_sigma(const Vec3& normal, PlaneKind plane, const std::vector<Atom>& atoms, double tol) {
  Operation op;
  op.kind = OperationKind::Reflection;
  op.order = 1;
  op.axis = normal.normalized();
  op.plane = plane;
  op.matrix = detail::reflection_matrix(op.axis);
  op.image = detail::permutation(atoms, op.matrix, tol);
  return op;
}

Operation make_sn(const Vec3& axis, int n, int power, const std::vector<Atom>& atoms, double tol) {
  Operation op;
  op.kind = OperationKind::ImproperRotation;
  op.order = n;
  op.power = power;
  op.axis = axis.normalized();
  op.matrix = detail::improper_rotation_matrix(op.axis, 2.0 * detail::kPi * power / n);
  op.image = detail::permutation(atoms, op.matrix, tol);
  return op;
}

bool atoms_collinear(const Molecule& mol, double tol) {
  if (mol.size() <= 2) {
    return true;
  }
  // Find two atoms farthest from origin / each other for axis
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

int highest_order_on_axis(const std::vector<Atom>& atoms, const Vec3& axis, int max_n, double tol) {
  int best = 1;
  for (int n = max_n; n >= 2; --n) {
    Mat3 R = detail::rotation_matrix(axis, 2.0 * detail::kPi / n);
    if (detail::is_symmetry(atoms, R, tol)) {
      return n;
    }
  }
  return best;
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

  // Polyhedral
  if (spherical_top && n_c3 >= 4) {
    if (n_c5 >= 6) {
      return has_inversion ? "Ih" : "I";
    }
    if (n_c4 >= 3) {
      return has_inversion ? "Oh" : "O";
    }
    // T family
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
    for (const Operation& op : ops) {
      if (op.kind == OperationKind::Reflection) {
        has_sigma_d = true;
        break;
      }
    }
    bool has_s2n = false;
    for (const Operation& op : ops) {
      if (op.kind == OperationKind::ImproperRotation && op.order == 2 * highest_cn) {
        has_s2n = true;
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
    // Sn?
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
  if (!name.empty() && (name[0] == 'C' || name[0] == 'S')) {
    return std::max(n, 1);
  }
  return 1;
}

void align_to_frame(Molecule& mol,
                    std::vector<Operation>& ops,
                    Mat3& orientation,
                    const Tolerances& tol) {
  // Find highest-order Cn axis
  Vec3 primary = Vec3::UnitZ();
  int best_n = 1;
  for (const Operation& op : ops) {
    if (op.kind == OperationKind::ProperRotation && op.order > best_n) {
      best_n = op.order;
      primary = op.axis.normalized();
    }
  }
  if (best_n == 1) {
    // Use principal axes already applied; nothing more
    orientation = Mat3::Identity();
    return;
  }

  // Secondary direction: C2' or atom in a plane
  Vec3 secondary = Vec3::UnitX();
  bool found = false;
  for (const Operation& op : ops) {
    if (op.kind == OperationKind::ProperRotation && op.order == 2) {
      if (std::abs(op.axis.normalized().dot(primary)) < 0.2) {
        secondary = op.axis.normalized();
        found = true;
        break;
      }
    }
  }
  if (!found) {
    for (const Operation& op : ops) {
      if (op.kind == OperationKind::Reflection) {
        // in-plane direction perpendicular to primary and to normal×primary
        Vec3 inplane = op.axis.cross(primary);
        if (inplane.norm() > 1e-8) {
          secondary = inplane.normalized();
          found = true;
          break;
        }
      }
    }
  }
  if (!found) {
    for (const Atom& a : mol.atoms) {
      Vec3 proj = a.position - a.position.dot(primary) * primary;
      if (proj.norm() > tol.primary) {
        secondary = proj.normalized();
        break;
      }
    }
  }

  Vec3 z = primary.normalized();
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
  orientation = R;

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
    result.molecule = molecule;
    result.molecule.resolve_masses();
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
  // Sort eigenvalues ascending, keep eigenvectors
  std::array<int, 3> order{{0, 1, 2}};
  std::sort(order.begin(), order.end(), [&](int a, int b) { return evals(a) < evals(b); });
  Mat3 sorted = Mat3::Zero();
  Vec3 sevals = Vec3::Zero();
  for (int i = 0; i < 3; ++i) {
    sevals(i) = evals(order[static_cast<std::size_t>(i)]);
    sorted.col(i) = evecs.col(order[static_cast<std::size_t>(i)]);
  }
  if (sorted.determinant() < 0) {
    sorted.col(0) *= -1;
  }
  result.principal_axes = sorted;

  // Rotate into principal-axis frame (columns of sorted are lab axes)
  Mat3 to_frame = sorted.transpose();
  for (Atom& a : mol.atoms) {
    a.position = to_frame * a.position;
  }

  const double imean =
      (std::abs(sevals(0)) + std::abs(sevals(1)) + std::abs(sevals(2))) / 3.0 + 1e-12;
  const bool spherical = std::abs(sevals(0) - sevals(1)) / imean < tol.inertia &&
                         std::abs(sevals(1) - sevals(2)) / imean < tol.inertia;
  const bool linear_inertia =
      sevals(0) / imean < tol.inertia && std::abs(sevals(1) - sevals(2)) / imean < tol.inertia;
  const bool collinear = atoms_collinear(mol, tol.primary * 10.0);
  result.linear = linear_inertia || collinear;

  result.equivalence_sets = build_equivalence_sets(mol, tol.primary);

  std::vector<Operation> ops;
  // Identity
  {
    Operation e;
    e.kind = OperationKind::Identity;
    e.matrix = Mat3::Identity();
    e.image = detail::permutation(mol.atoms, e.matrix, tol.primary);
    ops.push_back(e);
  }

  // Inversion
  if (detail::is_symmetry(mol.atoms, detail::inversion_matrix(), tol.primary)) {
    Operation inv;
    inv.kind = OperationKind::Inversion;
    inv.matrix = detail::inversion_matrix();
    inv.image = detail::permutation(mol.atoms, inv.matrix, tol.primary);
    ops.push_back(inv);
  }

  std::vector<Vec3> axes;
  add_unique_axis(axes, Vec3::UnitX(), tol.angle);
  add_unique_axis(axes, Vec3::UnitY(), tol.angle);
  add_unique_axis(axes, Vec3::UnitZ(), tol.angle);

  for (const EquivalenceSet& eq_set : result.equivalence_sets) {
    for (std::size_t i = 0; i < eq_set.atoms.size(); ++i) {
      const Vec3& pi = mol.atoms[static_cast<std::size_t>(eq_set.atoms[i])].position;
      if (pi.norm() > tol.primary) {
        add_unique_axis(axes, pi, tol.angle);
      }
      for (std::size_t j = i + 1; j < eq_set.atoms.size(); ++j) {
        const Vec3& pj = mol.atoms[static_cast<std::size_t>(eq_set.atoms[j])].position;
        add_unique_axis(axes, pi + pj, tol.angle);
        add_unique_axis(axes, pi.cross(pj), tol.angle);
        add_unique_axis(axes, pj - pi, tol.angle);
      }
    }
  }

  if (result.linear) {
    // Keep finite generators: C2 perp, σv, and i/σh
    // Ensure molecule along z
    if (collinear) {
      int i0 = 0;
      double rmax = 0;
      for (int i = 0; i < mol.size(); ++i) {
        double r = mol.atoms[static_cast<std::size_t>(i)].position.norm();
        if (r > rmax) {
          rmax = r;
          i0 = i;
        }
      }
      if (rmax > tol.primary) {
        Vec3 axis = mol.atoms[static_cast<std::size_t>(i0)].position.normalized();
        Vec3 y = axis.cross(Vec3::UnitX());
        if (y.norm() < 1e-8) {
          y = axis.cross(Vec3::UnitY());
        }
        y.normalize();
        Vec3 x = y.cross(axis);
        Mat3 R;
        R.row(0) = x.transpose();
        R.row(1) = y.transpose();
        R.row(2) = axis.transpose();
        for (Atom& a : mol.atoms) {
          a.position = R * a.position;
        }
        to_frame = R * to_frame;
      }
    }
    Operation c2 = make_cn(Vec3::UnitX(), 2, 1, mol.atoms, tol.primary);
    if (!c2.image.empty()) {
      ops.push_back(c2);
    }
    Operation sig = make_sigma(Vec3::UnitY(), PlaneKind::Vertical, mol.atoms, tol.primary);
    if (!sig.image.empty()) {
      ops.push_back(sig);
    }
    Operation sigh = make_sigma(Vec3::UnitZ(), PlaneKind::Horizontal, mol.atoms, tol.primary);
    if (!sigh.image.empty()) {
      ops.push_back(sigh);
    }
  } else {
    for (const Vec3& axis : axes) {
      const int n = highest_order_on_axis(mol.atoms, axis, tol.max_axis_order, tol.primary);
      if (n >= 2) {
        for (int k = 1; k < n; ++k) {
          Operation op = make_cn(axis, n, k, mol.atoms, tol.primary);
          if (!op.image.empty()) {
            ops.push_back(op);
          }
        }
        // Sn along same axis
        for (int sn = 2; sn <= 2 * n; sn += 2) {
          Operation op = make_sn(axis, sn, 1, mol.atoms, tol.primary);
          if (!op.image.empty()) {
            ops.push_back(op);
          }
        }
      }
      // Reflection with this normal
      Operation sig = make_sigma(axis, PlaneKind::None, mol.atoms, tol.primary);
      if (!sig.image.empty()) {
        ops.push_back(sig);
      }
    }

    // Pair midplanes
    for (const EquivalenceSet& eq_set : result.equivalence_sets) {
      for (std::size_t i = 0; i < eq_set.atoms.size(); ++i) {
        for (std::size_t j = i + 1; j < eq_set.atoms.size(); ++j) {
          const Vec3& pi = mol.atoms[static_cast<std::size_t>(eq_set.atoms[i])].position;
          const Vec3& pj = mol.atoms[static_cast<std::size_t>(eq_set.atoms[j])].position;
          Vec3 nrm = (pi - pj);
          if (nrm.norm() < tol.primary) {
            continue;
          }
          Operation sig = make_sigma(nrm, PlaneKind::None, mol.atoms, tol.primary);
          if (!sig.image.empty()) {
            ops.push_back(sig);
          }
        }
      }
    }
  }

  // Deduplicate before closing
  {
    std::vector<Operation> uniq;
    for (const Operation& op : ops) {
      bool found = false;
      for (const Operation& u : uniq) {
        if (detail::matrices_equal(u.matrix, op.matrix, tol.primary * 10.0)) {
          found = true;
          break;
        }
      }
      if (!found) {
        uniq.push_back(op);
      }
    }
    ops.swap(uniq);
  }

  if (!result.linear) {
    detail::close_group(ops, mol.atoms, tol.primary, tol.max_group_order);
  }

  // Recount features
  bool has_inversion = false;
  bool has_sigma_h = false;
  bool has_sigma_v = false;
  int highest_cn = 1;
  Vec3 primary_axis = Vec3::UnitZ();
  int n_c3 = 0;
  int n_c4 = 0;
  int n_c5 = 0;
  std::vector<Vec3> c2_axes;

  for (Operation& op : ops) {
    if (op.kind == OperationKind::Inversion) {
      has_inversion = true;
    }
    if (op.kind == OperationKind::ProperRotation && op.order > highest_cn) {
      highest_cn = op.order;
      primary_axis = op.axis.normalized();
    }
    if (op.kind == OperationKind::ProperRotation && op.order == 2) {
      add_unique_axis(c2_axes, op.axis, tol.angle);
    }
  }
  // Count unique high-order axes
  std::vector<Vec3> c3a;
  std::vector<Vec3> c4a;
  std::vector<Vec3> c5a;
  for (const Operation& op : ops) {
    if (op.kind != OperationKind::ProperRotation) {
      continue;
    }
    if (op.order == 3) {
      add_unique_axis(c3a, op.axis, tol.angle);
    }
    if (op.order == 4) {
      add_unique_axis(c4a, op.axis, tol.angle);
    }
    if (op.order == 5) {
      add_unique_axis(c5a, op.axis, tol.angle);
    }
  }
  n_c3 = static_cast<int>(c3a.size());
  n_c4 = static_cast<int>(c4a.size());
  n_c5 = static_cast<int>(c5a.size());

  int n_perp_c2 = 0;
  for (const Vec3& ax : c2_axes) {
    if (std::abs(ax.dot(primary_axis)) < 0.2) {
      ++n_perp_c2;
    }
  }

  // Label reflection planes relative to primary axis
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

  // Parse n from name
  int n_group = highest_cn;
  if (result.point_group.size() >= 2 &&
      (result.point_group[0] == 'C' || result.point_group[0] == 'D' ||
       result.point_group[0] == 'S')) {
    std::size_t i = 1;
    if (std::isdigit(static_cast<unsigned char>(result.point_group[i]))) {
      n_group = 0;
      while (i < result.point_group.size() &&
             std::isdigit(static_cast<unsigned char>(result.point_group[i]))) {
        n_group = n_group * 10 + (result.point_group[i] - '0');
        ++i;
      }
    }
  }

  Mat3 extra = Mat3::Identity();
  if (!result.linear && result.point_group != "C1" && result.point_group != "Cs" &&
      result.point_group != "Ci") {
    detail::align_to_frame(mol, ops, extra, tol);
  }
  result.orientation = extra * to_frame;
  result.molecule = std::move(mol);
  result.operations = std::move(ops);
  result.infinite =
      result.point_group == "Cinfv" || result.point_group == "Dinfh" || result.point_group == "Kh";
  result.order = result.infinite ? 0 : static_cast<int>(result.operations.size());
  result.symmetry_number = detail::symmetry_number_for(result.point_group, n_group);

  if (!result.infinite) {
    detail::assign_conjugacy_classes(result.operations, tol.primary);
  }

  result.character_table =
      detail::build_character_table(result.point_group, n_group, result.operations);
  result.subgroups = detail::constructive_subgroups(result.point_group, n_group);

  result.rmsd = result.symmetrize();

  return result;
}

} // namespace molsymm
