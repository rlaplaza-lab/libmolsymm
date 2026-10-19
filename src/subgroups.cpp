#include "internal.hpp"

#include "molsymm/analysis.hpp"
#include "molsymm/error.hpp"

#include <algorithm>

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
                   "D2h", "D2d", "D3", "D3d", "D4", "D4h", "Td",  "Th",  "T",   "O",   "Oh"})
      add(g);
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
      add("Ci");
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
    if (n % 4 == 0) {
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
    if (name.back() == 'h') {
      add("Cs");
      add("Ci");
      add(cnh(n));
      add(cnv(n));
      add(dn(n));
      add(dnh(n));
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
    } else {
      add(dn(n));
    }
    return s;
  }
  add(name);
  return s;
}

} // namespace molsymm::detail

namespace molsymm {

void Analysis::select_subgroup(const std::string& name) {
  auto it = std::find(subgroups.begin(), subgroups.end(), name);
  if (it == subgroups.end()) {
    throw Error("subgroup not listed: " + name);
  }
  // Restrict operations to those generating the subgroup is non-trivial.
  // For practical use: if requesting C1, keep only identity; if same as current, no-op.
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
  // Filter by matching classified subgroup name via re-detect on same geometry with
  // forced operation subset is heavy; keep operations that belong to the subgroup
  // family by regenerating from name for common cases.
  point_group = name;
  int n = 1;
  for (char c : name) {
    if (c >= '0' && c <= '9') {
      n = n * 10 + (c - '0');
    }
  }
  // Re-parse properly
  n = 0;
  for (std::size_t i = 0; i < name.size(); ++i) {
    if (name[i] >= '0' && name[i] <= '9') {
      n = 0;
      while (i < name.size() && name[i] >= '0' && name[i] <= '9') {
        n = n * 10 + (name[i] - '0');
        ++i;
      }
      break;
    }
  }
  // Keep identity + operations compatible with reduced group by order constraints
  std::vector<Operation> kept;
  for (const Operation& op : operations) {
    if (op.kind == OperationKind::Identity) {
      kept.push_back(op);
      continue;
    }
    if (name == "Cs" && op.kind == OperationKind::Reflection) {
      kept.push_back(op);
      break;
    }
    if (name == "Ci" && op.kind == OperationKind::Inversion) {
      kept.push_back(op);
      continue;
    }
    if (!name.empty() && name[0] == 'C' && name.back() != 's' && name.back() != 'i') {
      if (op.kind == OperationKind::ProperRotation && n > 0 && op.order <= n &&
          (n % op.order == 0)) {
        if (std::abs(op.axis.normalized().dot(Vec3::UnitZ())) > 0.9) {
          kept.push_back(op);
        }
      }
      if (name.back() == 'v' && op.kind == OperationKind::Reflection &&
          op.plane == PlaneKind::Vertical)
        kept.push_back(op);
      if (name.back() == 'h' &&
          ((op.kind == OperationKind::Reflection && op.plane == PlaneKind::Horizontal) ||
           op.kind == OperationKind::Inversion))
        kept.push_back(op);
    }
    if (!name.empty() && name[0] == 'D') {
      if (op.kind == OperationKind::ProperRotation) {
        kept.push_back(op);
      }
      if (name.back() == 'h' &&
          (op.kind == OperationKind::Reflection || op.kind == OperationKind::Inversion))
        kept.push_back(op);
      if (name.back() == 'd' &&
          (op.kind == OperationKind::Reflection || op.kind == OperationKind::ImproperRotation))
        kept.push_back(op);
    }
  }
  // Deduplicate
  std::vector<Operation> uniq;
  for (const Operation& op : kept) {
    bool found = false;
    for (const Operation& u : uniq) {
      if (detail::matrices_equal(u.matrix, op.matrix, 1e-6)) {
        found = true;
        break;
      }
    }
    if (!found) {
      uniq.push_back(op);
    }
  }
  operations = uniq;
  order = infinite ? 0 : static_cast<int>(operations.size());
  detail::assign_conjugacy_classes(operations, 1e-6);
  character_table = detail::build_character_table(point_group, n, operations);
  subgroups = detail::constructive_subgroups(point_group, n);
  symmetry_number = detail::symmetry_number_for(point_group, n);
}

} // namespace molsymm
