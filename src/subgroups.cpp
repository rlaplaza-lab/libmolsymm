#include "internal.hpp"

#include "molsymm/analysis.hpp"
#include "molsymm/error.hpp"

#include <algorithm>
#include <cmath>
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
    for (auto g : {"C2", "C3", "S4", "D2", "C2v", "C3v", "D2d", "T", "Td"}) {
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
    for (auto g : {"C2", "C3", "Ci", "Cs", "C2h", "C3h", "D2", "D2h", "T", "Th"}) {
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
    for (auto g : {"C2",  "C3",  "C4", "Ci",  "Cs", "C2v", "C2h", "C3v", "C4v", "C4h", "D2",
                   "D2h", "D2d", "D3", "D3d", "D4", "D4h", "Td",  "Th",  "T",   "O",   "Oh"}) {
      add(g);
    }
    return s;
  }
  if (name == "I" || name == "Ih") {
    for (auto g : {"C2", "C3", "C5", "D2", "D3", "D5", "T", "I"}) {
      add(g);
    }
    if (name == "Ih") {
      for (auto g : {"Ci", "Cs", "C2h", "C3v", "C5v", "D2h", "D3d", "D5d", "Th", "Ih"}) {
        add(g);
      }
    }
    return s;
  }

  // Axial families
  if (!name.empty() && name[0] == 'C') {
    if (n >= 2) {
      for (int k = 1; k <= n; ++k) {
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
      for (int k = 1; k < n; ++k) {
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
    } else if (n >= 2) {
      add(cn(n));
    }
    return s;
  }
  if (!name.empty() && name[0] == 'S' && n >= 2) {
    add(cn(n / 2));
    add("S" + std::to_string(n));
    if (n % 4 == 2) {
      add("C" + std::to_string(n / 2) + "h");
    }
    return s;
  }
  if (!name.empty() && name[0] == 'D') {
    for (int k = 1; k <= n; ++k) {
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
      for (int k = 1; k < n; ++k) {
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

// Slightly loose: parent geometries may be imperfect (e.g. cyclopropane at 1e-2 Å).
constexpr double kSubgroupMatTol = 5e-3;

bool is_cn_generator(const Operation& op, int cn) {
  if (op.kind != OperationKind::ProperRotation || cn < 2 || !detail::axis_on_z(op)) {
    return false;
  }
  // Primitive Cn^1, including Cn stored as a higher-axis power (C6^2 → C3).
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

const Operation* find_in_parent(const std::vector<Operation>& parent, const Mat3& M) {
  for (const Operation& op : parent) {
    if (detail::matrices_equal(op.matrix, M, kSubgroupMatTol)) {
      return &op;
    }
  }
  return nullptr;
}

/// Close a generating set inside the parent by matching products to parent elements.
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

std::vector<Operation>
seed_generators(const std::string& name, int n, const std::vector<Operation>& parent) {
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

  if (name == "Cs") {
    add_first([](const Operation& op) { return op.kind == OperationKind::Reflection; });
    return gens;
  }
  if (name == "Ci") {
    add_first([](const Operation& op) { return op.kind == OperationKind::Inversion; });
    return gens;
  }
  if (name == "T") {
    // C3 + face-axis C2 (parallel to a C4 in O/Oh); else any C2 (Td/Th).
    add_first([](const Operation& op) {
      return op.kind == OperationKind::ProperRotation && op.order == 3;
    });
    if (!add_first([&](const Operation& op) {
          return op.kind == OperationKind::ProperRotation && op.order == 2 &&
                 detail::axis_parallel_to_c4(op.axis, parent);
        })) {
      add_first([](const Operation& op) {
        return op.kind == OperationKind::ProperRotation && op.order == 2;
      });
    }
    return gens;
  }
  if (!name.empty() && name[0] == 'S' && n >= 2) {
    add_first([&](const Operation& op) {
      return op.kind == OperationKind::ImproperRotation && detail::axis_on_z(op) && op.order == n &&
             op.power == 1;
    });
    return gens;
  }
  if (!name.empty() && name[0] == 'C' && name.back() != 's' && name.back() != 'i') {
    add_first([&](const Operation& op) { return is_cn_generator(op, n); });
    if (name.back() == 'v') {
      add_first([](const Operation& op) {
        return op.kind == OperationKind::Reflection &&
               (op.plane == PlaneKind::Vertical || op.plane == PlaneKind::Dihedral);
      });
    }
    if (name.back() == 'h') {
      add_first([](const Operation& op) {
        return op.kind == OperationKind::Reflection && op.plane == PlaneKind::Horizontal;
      });
    }
    return gens;
  }
  if (!name.empty() && name[0] == 'D') {
    if (n >= 2) {
      add_first([&](const Operation& op) { return is_cn_generator(op, n); });
    }
    add_first([&](const Operation& op) {
      return op.kind == OperationKind::ProperRotation && op.order == 2 && !detail::axis_on_z(op);
    });
    if (name.back() == 'h') {
      add_first([](const Operation& op) {
        return op.kind == OperationKind::Reflection && op.plane == PlaneKind::Horizontal;
      });
    }
    if (name.back() == 'd') {
      add_first([](const Operation& op) {
        return op.kind == OperationKind::Reflection &&
               (op.plane == PlaneKind::Dihedral || op.plane == PlaneKind::Vertical);
      });
    }
    return gens;
  }
  return gens;
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
    character_table = detail::build_character_table("C1", 1, operations);
    subgroups = detail::constructive_subgroups("C1", 1);
    return;
  }

  const int n = detail::parse_n(name);
  std::vector<Operation> gens = seed_generators(name, n, operations);
  std::vector<Operation> kept = close_inside_parent(operations, gens);

  const int expect = detail::expected_group_order(name, n);
  if (expect < 0 || static_cast<int>(kept.size()) != expect) {
    throw Error("cannot realize subgroup " + name + " (have " + std::to_string(kept.size()) +
                " operations, need " + std::to_string(expect) + ")");
  }

  operations = std::move(kept);
  point_group = name;
  infinite = false;
  order = static_cast<int>(operations.size());
  detail::assign_conjugacy_classes(operations, 1e-6);
  character_table = detail::build_character_table(point_group, n, operations);
  subgroups = detail::constructive_subgroups(point_group, n);
  symmetry_number = detail::symmetry_number_for(point_group, n);
}

} // namespace molsymm
