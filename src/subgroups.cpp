#include "internal.hpp"

#include "molsymm/analysis.hpp"
#include "molsymm/error.hpp"

#include <algorithm>
#include <cmath>
#include <ranges>
#include <utility>

namespace molsymm::detail {

std::vector<std::string> constructive_subgroups(const std::string& name, int n) {
  std::vector<std::string> s{"C1"};
  auto add = [&](const std::string& g) {
    if (std::find(s.begin(), s.end(), g) == s.end()) {
      s.push_back(g);
    }
  };
  auto cn = [](int k) { return "C" + std::to_string(k); };
  auto cnv = [](int k) { return "C" + std::to_string(k) + "v"; };
  auto cnh = [](int k) { return "C" + std::to_string(k) + "h"; };
  auto dn = [](int k) { return "D" + std::to_string(k); };
  auto dnh = [](int k) { return "D" + std::to_string(k) + "h"; };
  auto dnd = [](int k) { return "D" + std::to_string(k) + "d"; };

  if (name == "Cs" || name == "Ci") {
    add(name);
    return s;
  }
  if (name == "Cinfv") {
    add("Cs");
    add("Cinfv");
    return s;
  }
  if (name == "Dinfh") {
    add("Cs");
    add("Ci");
    add("Cinfv");
    add("Dinfh");
    return s;
  }
  if (name == "Kh") {
    add("Cs");
    add("Ci");
    add("Kh");
    return s;
  }
  if (name == "Td") {
    for (auto g : {"C2", "C3", "Cs", "S4", "D2", "C2v", "C3v", "D2d", "T", "Td"}) {
      add(g);
    }
    return s;
  }
  if (name == "T") {
    for (auto g : {"C2", "C3", "D2", "T"}) {
      add(g);
    }
    return s;
  }
  if (name == "Th") {
    for (auto g : {"C2", "C3", "Ci", "Cs", "C2v", "C2h", "S6", "D2", "D2h", "T", "Th"}) {
      add(g);
    }
    return s;
  }
  if (name == "O") {
    for (auto g : {"C2", "C3", "C4", "D2", "D3", "D4", "T", "O"}) {
      add(g);
    }
    return s;
  }
  if (name == "Oh") {
    for (auto g :
         {"C2", "C3",  "C4",  "Ci", "Cs",  "C2v", "C2h", "C3v", "C4v", "C4h", "S4", "S6",
          "D2", "D2h", "D2d", "D3", "D3d", "D4",  "D4h", "Td",  "Th",  "T",   "O",  "Oh"}) {
      add(g);
    }
    return s;
  }
  if (name == "I" || name == "Ih") {
    for (auto g : {"C2", "C3", "C5", "D2", "D3", "D5", "T", "I"}) {
      add(g);
    }
    if (name == "Ih") {
      for (auto g :
           {"Ci", "Cs", "C2v", "C2h", "C3v", "C5v", "S6", "S10", "D2h", "D3d", "D5d", "Th", "Ih"}) {
        add(g);
      }
    }
    return s;
  }

  // Axial families. Divisor loops start at 2 so C1v / D1 / D1h are never listed.
  if (!name.empty() && name[0] == 'C') {
    if (n >= 2) {
      for (int k = 2; k <= n; ++k) {
        if (n % k == 0) {
          add(cn(k));
        }
      }
    }
    if (name.back() == 'v') {
      add("Cs");
      if (n >= 2) {
        add(cnv(n));
      }
      for (int k = 2; k < n; ++k) {
        if (n % k == 0) {
          add(cnv(k));
        }
      }
    } else if (name.back() == 'h') {
      add("Cs");
      if (n % 2 == 0) {
        add("Ci");
      }
      if (n >= 2) {
        add(cnh(n));
      }
      add(cn(n));
      for (int k = 2; k < n; ++k) {
        if (n % k == 0) {
          add(cnh(k));
          add(cn(k));
        }
      }
    } else if (n >= 2) {
      add(cn(n));
    }
    return s;
  }
  if (!name.empty() && name[0] == 'S' && n >= 2) {
    if (n / 2 >= 2) {
      add(cn(n / 2));
    }
    add("S" + std::to_string(n));
    if (n % 4 == 2) {
      add("Ci");
    }
    return s;
  }
  if (!name.empty() && name[0] == 'D') {
    for (int k = 2; k <= n; ++k) {
      if (n % k == 0) {
        add(cn(k));
        add(dn(k));
      }
    }
    add("C2");
    if (name.back() == 'h') {
      add("Cs");
      add(cnh(n));
      add(cnv(n));
      add(dn(n));
      add(dnh(n));
      add("C2v");
      if (n % 2 == 0) {
        add("Ci");
        add("C2h");
        add("D2h");
      }
      for (int k = 2; k < n; ++k) {
        if (n % k == 0) {
          add(dnh(k));
          add(cnv(k));
        }
      }
    } else if (name.back() == 'd') {
      add("Cs");
      add(cnv(n));
      add(dn(n));
      add(dnd(n));
      add("S" + std::to_string(2 * n));
      add("C2v");
      if (n % 2 == 1) {
        add("Ci");
        add("C2h");
      } else {
        add("D2d");
      }
    } else {
      add(dn(n));
    }
    return s;
  }
  add(name);
  return s;
}

static int expected_group_order(const std::string& name, int n) {
  if (name == "C1") {
    return 1;
  }
  if (name == "Cs" || name == "Ci") {
    return 2;
  }
  if (name == "T") {
    return 12;
  }
  if (name == "Td" || name == "Th" || name == "O") {
    return 24;
  }
  if (name == "Oh") {
    return 48;
  }
  if (name == "I") {
    return 60;
  }
  if (name == "Ih") {
    return 120;
  }
  if (name.empty()) {
    return -1;
  }
  if (name[0] == 'S' && n > 0) {
    return n;
  }
  if (name[0] == 'C' && n > 0) {
    if (name.back() == 'v' || name.back() == 'h') {
      return 2 * n;
    }
    return n;
  }
  if (name[0] == 'D' && n > 0) {
    if (name.back() == 'h' || name.back() == 'd') {
      return 4 * n;
    }
    return 2 * n;
  }
  return -1;
}

} // namespace molsymm::detail

namespace molsymm {
namespace {

constexpr double kSubgroupMatTol = 5e-3;
constexpr double kAxisDotTol = 0.2;

bool axes_parallel(const Vec3& a, const Vec3& b) {
  if (a.norm() < 1e-12 || b.norm() < 1e-12) {
    return false;
  }
  return std::abs(a.normalized().dot(b.normalized())) > 0.9;
}

bool axes_perpendicular(const Vec3& a, const Vec3& b) {
  if (a.norm() < 1e-12 || b.norm() < 1e-12) {
    return false;
  }
  return std::abs(a.normalized().dot(b.normalized())) < kAxisDotTol;
}

bool is_180_rotation(const Operation& op) {
  if (op.kind != OperationKind::ProperRotation) {
    return false;
  }
  if (op.order == 2) {
    return true;
  }
  // Face-axis C4^2 in O/Oh is stored with order 4, power 2, trace −1.
  return op.order >= 2 && ((op.order % 2 == 0 && op.power == op.order / 2) ||
                           std::abs(op.matrix.trace() + 1.0) < 1e-4);
}

bool is_cn_generator(const Operation& op, int cn) {
  if (op.kind != OperationKind::ProperRotation || cn < 2 || op.axis.norm() < 1e-12) {
    return false;
  }
  if (op.order == cn && (op.power == 1 || op.power == cn - 1)) {
    return true;
  }
  if (op.order > cn && op.order % cn == 0 && op.power > 0) {
    const int step = op.order / cn;
    return op.power == step || op.power == op.order - step;
  }
  const double ang = std::acos(std::clamp((op.matrix.trace() - 1.0) * 0.5, -1.0, 1.0));
  return std::abs(ang - 2.0 * detail::kPi / cn) < 1e-3;
}

bool is_sn_generator(const Operation& op, int sn) {
  return op.kind == OperationKind::ImproperRotation && op.axis.norm() > 1e-12 && op.order == sn &&
         (op.power == 1 || op.power == sn - 1);
}

const Operation* find_in_parent(const std::vector<Operation>& parent, const Mat3& M) {
  for (const Operation& op : parent) {
    if (detail::matrices_equal(op.matrix, M, kSubgroupMatTol)) {
      return &op;
    }
  }
  return nullptr;
}

std::vector<Operation> close_inside_parent(const std::vector<Operation>& parent,
                                           const std::vector<Operation>& gens) {
  std::vector<Operation> kept;
  kept.reserve(parent.size());
  for (const Operation& g : gens) {
    const Operation* hit = find_in_parent(parent, g.matrix);
    if (hit != nullptr) {
      kept.push_back(*hit);
    }
  }
  for (const Operation& op : parent) {
    if (op.kind == OperationKind::Identity) {
      kept.push_back(op);
      break;
    }
  }
  detail::dedupe_operations(kept, kSubgroupMatTol);

  bool grew = true;
  while (grew) {
    grew = false;
    const int n = static_cast<int>(kept.size());
    for (int i = 0; i < n && !grew; ++i) {
      for (int j = 0; j < n; ++j) {
        const Mat3 M =
            kept[static_cast<std::size_t>(i)].matrix * kept[static_cast<std::size_t>(j)].matrix;
        if (find_in_parent(kept, M) != nullptr) {
          continue;
        }
        const Operation* hit = find_in_parent(parent, M);
        if (hit == nullptr) {
          continue;
        }
        kept.push_back(*hit);
        grew = true;
        break;
      }
    }
  }
  return kept;
}

struct KindCounts {
  int identity = 0;
  int proper = 0;
  int improper = 0;
  int reflection = 0;
  int inversion = 0;
};

KindCounts count_kinds(const std::vector<Operation>& ops) {
  KindCounts c;
  for (const Operation& op : ops) {
    switch (op.kind) {
    case OperationKind::Identity:
      ++c.identity;
      break;
    case OperationKind::ProperRotation:
      ++c.proper;
      break;
    case OperationKind::ImproperRotation:
      ++c.improper;
      break;
    case OperationKind::Reflection:
      ++c.reflection;
      break;
    case OperationKind::Inversion:
      ++c.inversion;
      break;
    }
  }
  return c;
}

bool fingerprint_matches(const std::string& name, int n, const std::vector<Operation>& ops) {
  const KindCounts c = count_kinds(ops);
  if (c.identity != 1) {
    return false;
  }
  if (name == "Cs") {
    return c.reflection == 1 && c.proper == 0 && c.improper == 0 && c.inversion == 0;
  }
  if (name == "Ci") {
    return c.inversion == 1 && c.proper == 0 && c.improper == 0 && c.reflection == 0;
  }
  if (name == "T") {
    return c.proper == 11 && c.improper == 0 && c.reflection == 0 && c.inversion == 0;
  }
  if (name == "Td") {
    return c.proper == 11 && c.improper == 6 && c.reflection == 6 && c.inversion == 0;
  }
  if (name == "Th") {
    return c.proper == 11 && c.improper == 8 && c.reflection == 3 && c.inversion == 1;
  }
  if (name == "O") {
    return c.proper == 23 && c.improper == 0 && c.reflection == 0 && c.inversion == 0;
  }
  if (name == "Oh") {
    return c.proper == 23 && c.improper == 14 && c.reflection == 9 && c.inversion == 1;
  }
  if (name == "I") {
    return c.proper == 59 && c.improper == 0 && c.reflection == 0 && c.inversion == 0;
  }
  if (name == "Ih") {
    return c.proper == 59 && c.improper == 44 && c.reflection == 15 && c.inversion == 1;
  }
  if (name.empty() || n < 2) {
    return false;
  }
  if (name[0] == 'S') {
    return c.proper + c.improper + c.inversion == n - 1 && c.reflection == 0;
  }
  if (name[0] == 'C') {
    if (name.back() == 'v') {
      return c.proper == n - 1 && c.reflection == n && c.improper == 0 && c.inversion == 0;
    }
    if (name.back() == 'h') {
      if (n % 2 == 0) {
        return c.proper == n - 1 && c.reflection == 1 && c.inversion == 1 && c.improper == n - 2;
      }
      return c.proper == n - 1 && c.reflection == 1 && c.inversion == 0 && c.improper == n - 1;
    }
    return c.proper == n - 1 && c.reflection == 0 && c.improper == 0 && c.inversion == 0;
  }
  if (name[0] == 'D') {
    // Dn contributes 2n − 1 proper rotations (Cn powers + n perpendicular C2).
    if (name.back() == 'h') {
      if (n % 2 == 0) {
        // Even Dnh: S_n / S_n^{n-1} only among Sn (i and σh are separate) → n − 2 impropers.
        return c.proper == 2 * n - 1 && c.reflection == n + 1 && c.inversion == 1 &&
               c.improper == n - 2;
      }
      return c.proper == 2 * n - 1 && c.reflection == n + 1 && c.inversion == 0 &&
             c.improper == n - 1;
    }
    if (name.back() == 'd') {
      if (n % 2 == 1) {
        return c.proper == 2 * n - 1 && c.reflection == n && c.inversion == 1 &&
               c.improper == n - 1;
      }
      return c.proper == 2 * n - 1 && c.reflection == n && c.inversion == 0 && c.improper == n;
    }
    return c.proper == 2 * n - 1 && c.reflection == 0 && c.improper == 0 && c.inversion == 0;
  }
  return false;
}

void apply_frame(Molecule& mol, std::vector<Operation>& ops, Mat3& orientation, const Mat3& Q) {
  for (Atom& a : mol.atoms) {
    a.position = Q * a.position;
  }
  for (Operation& op : ops) {
    op.matrix = Q * op.matrix * Q.transpose();
    op.axis = Q * op.axis;
    if (op.axis.norm() > 1e-12) {
      op.axis.normalize();
    }
  }
  orientation = Q * orientation;
}

Mat3 basis_z_from_axis(Vec3 z, Vec3 secondary) {
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

bool needs_principal_z(const std::string& name) {
  if (name.empty()) {
    return false;
  }
  if (name == "T" || name == "Td" || name == "Th" || name == "O" || name == "Oh" || name == "I" ||
      name == "Ih" || name == "Cs" || name == "Ci" || name == "C1") {
    return false;
  }
  return name[0] == 'C' || name[0] == 'D' || name[0] == 'S';
}

Vec3 principal_axis_of(const std::string& name, int n, const std::vector<Operation>& ops) {
  if (name[0] == 'S') {
    for (const Operation& op : ops) {
      if (is_sn_generator(op, n)) {
        return op.axis.normalized();
      }
    }
  }
  int best_order = 0;
  Vec3 best = Vec3::UnitZ();
  for (const Operation& op : ops) {
    if (op.kind == OperationKind::ProperRotation && op.axis.norm() > 1e-12 &&
        op.order > best_order) {
      best_order = op.order;
      best = op.axis.normalized();
    }
  }
  if (best_order >= 2) {
    return best;
  }
  for (const Operation& op : ops) {
    if (is_180_rotation(op) && op.axis.norm() > 1e-12) {
      return op.axis.normalized();
    }
  }
  return Vec3::UnitZ();
}

void reorient_principal_to_z(Analysis& a, const std::string& name, int n) {
  if (!needs_principal_z(name)) {
    return;
  }
  const Vec3 z = principal_axis_of(name, n, a.operations);
  if (std::abs(z.z()) > 0.9) {
    return;
  }
  Vec3 secondary = Vec3::UnitX();
  for (const Operation& op : a.operations) {
    if (op.kind == OperationKind::ProperRotation && is_180_rotation(op) &&
        axes_perpendicular(op.axis, z)) {
      secondary = op.axis.normalized();
      break;
    }
    if (op.kind == OperationKind::Reflection && axes_perpendicular(op.axis, z)) {
      secondary = op.axis.cross(z).normalized();
      if (secondary.norm() > 1e-8) {
        break;
      }
    }
  }
  const Mat3 Q = basis_z_from_axis(z, secondary);
  apply_frame(a.molecule, a.operations, a.orientation, Q);

  // Retag vertical / horizontal / dihedral planes in the new frame.
  const int n_group = std::max(n, 2);
  const bool even_axial = (n_group % 2 == 0) && (name.back() == 'v' || name.back() == 'h' ||
                                                 name.back() == 'd' || name[0] == 'D');
  for (Operation& op : a.operations) {
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

bool try_realize(const std::string& name,
                 int n,
                 const std::vector<Operation>& parent,
                 const std::vector<Operation>& gens,
                 int expect,
                 std::vector<Operation>& out) {
  std::vector<Operation> kept = close_inside_parent(parent, gens);
  if (static_cast<int>(kept.size()) != expect) {
    return false;
  }
  if (!fingerprint_matches(name, n, kept)) {
    return false;
  }
  out = std::move(kept);
  return true;
}

bool realize_subgroup(const std::string& name,
                      int n,
                      const std::vector<Operation>& parent,
                      std::vector<Operation>& out) {
  const int expect = detail::expected_group_order(name, n);
  if (expect < 0) {
    return false;
  }

  auto attempt = [&](const std::vector<Operation>& gens) {
    return try_realize(name, n, parent, gens, expect, out);
  };

  if (name == "Cs") {
    return std::ranges::any_of(parent, [&](const Operation& op) {
      return op.kind == OperationKind::Reflection && attempt({op});
    });
  }
  if (name == "Ci") {
    return std::ranges::any_of(parent, [&](const Operation& op) {
      return op.kind == OperationKind::Inversion && attempt({op});
    });
  }

  if (name == "T") {
    for (const Operation& c3 : parent) {
      if (c3.kind != OperationKind::ProperRotation || c3.order != 3) {
        continue;
      }
      for (const Operation& c2 : parent) {
        if (!is_180_rotation(c2)) {
          continue;
        }
        // Prefer face-axis C2 (‖ a C4 in O/Oh); skip edge-axis C2 that generate O.
        const bool face =
            detail::axis_parallel_to_c4(c2.axis, parent) || (c2.order == 4 && is_180_rotation(c2));
        const bool td_c2 = c2.order == 2 && !detail::axis_parallel_to_c4(c2.axis, parent);
        if (!face && !td_c2) {
          continue;
        }
        if (face && attempt({c3, c2})) {
          return true;
        }
      }
      for (const Operation& c2 : parent) {
        if (c2.order == 2 && is_180_rotation(c2) && attempt({c3, c2})) {
          return true;
        }
      }
    }
    return false;
  }

  if (name == "Th") {
    const Operation* inv = nullptr;
    for (const Operation& op : parent) {
      if (op.kind == OperationKind::Inversion) {
        inv = &op;
        break;
      }
    }
    if (inv == nullptr) {
      return false;
    }
    for (const Operation& c3 : parent) {
      if (c3.kind != OperationKind::ProperRotation || c3.order != 3) {
        continue;
      }
      // Prefer face-axis C2 (stored as C4^2 in O/Oh); edge C2 does not generate Th.
      for (const Operation& c2 : parent) {
        if (!is_180_rotation(c2)) {
          continue;
        }
        const bool face =
            detail::axis_parallel_to_c4(c2.axis, parent) || (c2.order == 4 && is_180_rotation(c2));
        if (face && attempt({c3, *inv, c2})) {
          return true;
        }
      }
      for (const Operation& c2 : parent) {
        if (is_180_rotation(c2) && attempt({c3, *inv, c2})) {
          return true;
        }
      }
    }
    return false;
  }

  if (name == "Td" || name == "O" || name == "Oh" || name == "I" || name == "Ih") {
    std::vector<Operation> gens;
    auto add_first = [&](auto pred) {
      for (const Operation& op : parent) {
        if (pred(op)) {
          gens.push_back(op);
          return true;
        }
      }
      return false;
    };
    if (name == "Td") {
      add_first([](const Operation& op) {
        return op.kind == OperationKind::ProperRotation && op.order == 3;
      });
      if (!add_first([](const Operation& op) {
            return op.kind == OperationKind::ImproperRotation && op.order == 4;
          })) {
        add_first([](const Operation& op) { return op.kind == OperationKind::Reflection; });
      }
    } else if (name == "O") {
      add_first([](const Operation& op) {
        return op.kind == OperationKind::ProperRotation && op.order == 4 && op.power == 1;
      });
      add_first([](const Operation& op) {
        return op.kind == OperationKind::ProperRotation && op.order == 3;
      });
    } else if (name == "Oh") {
      add_first([](const Operation& op) {
        return op.kind == OperationKind::ProperRotation && op.order == 4 && op.power == 1;
      });
      add_first([](const Operation& op) {
        return op.kind == OperationKind::ProperRotation && op.order == 3;
      });
      add_first([](const Operation& op) { return op.kind == OperationKind::Inversion; });
    } else if (name == "I") {
      add_first([](const Operation& op) {
        return op.kind == OperationKind::ProperRotation && op.order == 5 && op.power == 1;
      });
      add_first([](const Operation& op) {
        return op.kind == OperationKind::ProperRotation && op.order == 3;
      });
    } else { // Ih
      add_first([](const Operation& op) {
        return op.kind == OperationKind::ProperRotation && op.order == 5 && op.power == 1;
      });
      add_first([](const Operation& op) {
        return op.kind == OperationKind::ProperRotation && op.order == 3;
      });
      add_first([](const Operation& op) { return op.kind == OperationKind::Inversion; });
    }
    return attempt(gens);
  }

  if (!name.empty() && name[0] == 'S' && n >= 2) {
    for (int pass = 0; pass < 2; ++pass) {
      const bool require_z = (pass == 0);
      for (const Operation& op : parent) {
        if (!is_sn_generator(op, n)) {
          continue;
        }
        if (require_z && !detail::axis_on_z(op)) {
          continue;
        }
        if (attempt({op})) {
          return true;
        }
      }
    }
    return false;
  }

  if (!name.empty() && name[0] == 'C' && name.back() != 's' && name.back() != 'i') {
    // Prefer the parent principal axis (already on z) so Dnh→C2v keeps the ring axis.
    for (int pass = 0; pass < 2; ++pass) {
      const bool require_z = (pass == 0);
      for (const Operation& cn : parent) {
        if (!is_cn_generator(cn, n)) {
          continue;
        }
        if (require_z && !detail::axis_on_z(cn)) {
          continue;
        }
        if (name.back() != 'v' && name.back() != 'h') {
          if (attempt({cn})) {
            return true;
          }
          continue;
        }
        for (const Operation& sig : parent) {
          if (sig.kind != OperationKind::Reflection || sig.axis.norm() < 1e-12) {
            continue;
          }
          if (name.back() == 'v' && axes_perpendicular(sig.axis, cn.axis) && attempt({cn, sig})) {
            return true;
          }
          if (name.back() == 'h' && axes_parallel(sig.axis, cn.axis) && attempt({cn, sig})) {
            return true;
          }
        }
      }
    }
    return false;
  }

  if (!name.empty() && name[0] == 'D') {
    for (int pass = 0; pass < 2; ++pass) {
      const bool require_z = (pass == 0);
      for (const Operation& cn : parent) {
        if (n >= 2 && !is_cn_generator(cn, n)) {
          continue;
        }
        if (n < 2) {
          continue;
        }
        if (require_z && !detail::axis_on_z(cn)) {
          continue;
        }
        for (const Operation& c2 : parent) {
          if (!is_180_rotation(c2) || !axes_perpendicular(c2.axis, cn.axis)) {
            continue;
          }
          if (name.back() != 'h' && name.back() != 'd') {
            if (attempt({cn, c2})) {
              return true;
            }
            continue;
          }
          for (const Operation& sig : parent) {
            if (sig.kind != OperationKind::Reflection || sig.axis.norm() < 1e-12) {
              continue;
            }
            if (name.back() == 'h') {
              if (axes_parallel(sig.axis, cn.axis) && attempt({cn, c2, sig})) {
                return true;
              }
            } else { // d
              if (axes_perpendicular(sig.axis, cn.axis) && attempt({cn, c2, sig})) {
                return true;
              }
            }
          }
        }
      }
    }
    // D2 has no unique Cn>2; try any pair of perpendicular C2 axes.
    if (n == 2 && (name == "D2" || name == "D2h" || name == "D2d")) {
      for (std::size_t i = 0; i < parent.size(); ++i) {
        if (!is_180_rotation(parent[i])) {
          continue;
        }
        for (std::size_t j = i + 1; j < parent.size(); ++j) {
          if (!is_180_rotation(parent[j]) || !axes_perpendicular(parent[i].axis, parent[j].axis)) {
            continue;
          }
          if (name == "D2") {
            if (attempt({parent[i], parent[j]})) {
              return true;
            }
            continue;
          }
          for (const Operation& sig : parent) {
            if (sig.kind != OperationKind::Reflection || sig.axis.norm() < 1e-12) {
              continue;
            }
            if (name == "D2h" && axes_parallel(sig.axis, parent[i].axis) &&
                attempt({parent[i], parent[j], sig})) {
              return true;
            }
            if (name == "D2d" && axes_perpendicular(sig.axis, parent[i].axis) &&
                attempt({parent[i], parent[j], sig})) {
              return true;
            }
          }
        }
      }
    }
    return false;
  }

  return false;
}

} // namespace

void Analysis::select_subgroup(const std::string& name) {
  auto it = std::find(subgroups.begin(), subgroups.end(), name);
  if (it == subgroups.end()) {
    throw Error("subgroup not listed: " + name);
  }
  if (name == point_group) {
    return;
  }
  if (name == "C1") {
    operations.erase(
        std::remove_if(operations.begin(),
                       operations.end(),
                       [](const Operation& op) { return op.kind != OperationKind::Identity; }),
        operations.end());
    point_group = "C1";
    order = 1;
    infinite = false;
    symmetry_number = 1;
    detail::assign_conjugacy_classes(operations, 1e-6);
    character_table = detail::build_character_table("C1", 1, operations);
    subgroups = detail::constructive_subgroups("C1", 1);
    return;
  }
  if (name == "Cinfv") {
    // Continuous axial group: keep the C∞v projector, drop the D∞h inversion average.
    operations.erase(std::remove_if(operations.begin(),
                                    operations.end(),
                                    [](const Operation& op) {
                                      return op.kind == OperationKind::Inversion ||
                                             (op.kind == OperationKind::Reflection &&
                                              op.plane == PlaneKind::Horizontal);
                                    }),
                     operations.end());
    point_group = "Cinfv";
    infinite = true;
    order = 0;
    symmetry_number = 1;
    character_table = detail::build_character_table("Cinfv", 0, operations);
    subgroups = detail::constructive_subgroups("Cinfv", 0);
    return;
  }

  const int n = detail::parse_n(name);
  std::vector<Operation> kept;
  if (!realize_subgroup(name, n, operations, kept)) {
    throw Error("cannot realize subgroup " + name);
  }

  operations = std::move(kept);
  point_group = name;
  infinite = false;
  order = static_cast<int>(operations.size());
  reorient_principal_to_z(*this, name, n);
  detail::assign_conjugacy_classes(operations, 1e-6);
  character_table = detail::build_character_table(point_group, n, operations);
  subgroups = detail::constructive_subgroups(point_group, n);
  symmetry_number = detail::symmetry_number_for(point_group, n);
}

} // namespace molsymm
