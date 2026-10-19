#include "internal.hpp"

#include "molsymm/error.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>

namespace molsymm::detail {
namespace {

void add_irrep(CharacterTable& ct,
               const std::string& name,
               int dim,
               const std::vector<double>& chars) {
  ct.irrep_names.push_back(name);
  ct.irrep_dims.push_back(dim);
  ct.characters.push_back(chars);
}

CharacterTable table_c1() {
  CharacterTable ct;
  ct.class_labels = {"E"};
  ct.class_sizes = {1};
  add_irrep(ct, "A", 1, {1});
  return ct;
}

CharacterTable table_cs() {
  CharacterTable ct;
  ct.class_labels = {"E", "σ"};
  ct.class_sizes = {1, 1};
  add_irrep(ct, "A'", 1, {1, 1});
  add_irrep(ct, "A''", 1, {1, -1});
  return ct;
}

CharacterTable table_ci() {
  CharacterTable ct;
  ct.class_labels = {"E", "i"};
  ct.class_sizes = {1, 1};
  add_irrep(ct, "Ag", 1, {1, 1});
  add_irrep(ct, "Au", 1, {1, -1});
  return ct;
}

CharacterTable table_cn(int n) {
  CharacterTable ct;
  ct.class_labels.emplace_back("E");
  ct.class_sizes.push_back(1);
  for (int k = 1; k < n; ++k) {
    std::ostringstream os;
    os << "C_" << n;
    if (k > 1) {
      os << "^" << k;
    }
    ct.class_labels.push_back(os.str());
    ct.class_sizes.push_back(1);
  }
  // A
  std::vector<double> a(n, 1.0);
  add_irrep(ct, "A", 1, a);
  if (n % 2 == 0) {
    std::vector<double> b(n);
    for (int k = 0; k < n; ++k) {
      b[static_cast<std::size_t>(k)] = (k % 2 == 0) ? 1.0 : -1.0;
    }
    add_irrep(ct, "B", 1, b);
  }
  const int n_e = (n - 1) / 2;
  for (int j = 1; j <= n_e; ++j) {
    std::vector<double> e(n);
    e[0] = 2.0;
    for (int k = 1; k < n; ++k) {
      e[static_cast<std::size_t>(k)] = 2.0 * std::cos(2.0 * kPi * j * k / n);
    }
    std::ostringstream os;
    os << "E";
    if (n_e > 1) {
      os << "_" << j;
    }
    add_irrep(ct, os.str(), 2, e);
  }
  return ct;
}

CharacterTable table_cnv(int n) {
  CharacterTable ct;
  // Classes: E, C_n^k (paired), σv (and σd if even)
  ct.class_labels.emplace_back("E");
  ct.class_sizes.push_back(1);
  for (int k = 1; k <= (n - 1) / 2; ++k) {
    std::ostringstream os;
    os << "2C_" << n;
    if (k > 1) {
      os << "^" << k;
    }
    ct.class_labels.push_back(os.str());
    ct.class_sizes.push_back(2);
  }
  if (n % 2 == 0) {
    ct.class_labels.emplace_back("C_2");
    ct.class_sizes.push_back(1);
    ct.class_labels.emplace_back("nσv");
    ct.class_sizes.push_back(n / 2);
    ct.class_labels.emplace_back("nσd");
    ct.class_sizes.push_back(n / 2);
  } else {
    ct.class_labels.emplace_back("nσv");
    ct.class_sizes.push_back(n);
  }

  const int nc = static_cast<int>(ct.class_labels.size());
  auto chars = [&](double e, const std::vector<double>& rot, double sigv, double sigd = 0) {
    std::vector<double> c;
    c.push_back(e);
    for (double v : rot) {
      c.push_back(v);
    }
    c.push_back(sigv);
    if (n % 2 == 0) {
      c.push_back(sigd);
    }
    while (static_cast<int>(c.size()) < nc) {
      c.push_back(0);
    }
    return c;
  };

  auto rot_vals = [&](int j, bool b_irrep) {
    std::vector<double> r;
    for (int k = 1; k <= (n - 1) / 2; ++k) {
      if (b_irrep) {
        r.push_back((k % 2 == 0) ? 1.0 : -1.0);
      } else if (j == 0) {
        r.push_back(1.0);
      } else {
        r.push_back(2.0 * std::cos(2.0 * kPi * j * k / n));
      }
    }
    if (n % 2 == 0) {
      if (b_irrep) {
        r.push_back(-1.0);
      } else if (j == 0) {
        r.push_back(1.0);
      } else {
        r.push_back(2.0 * std::cos(2.0 * kPi * j * (0.5 * n) / n));
      }
    }
    return r;
  };

  add_irrep(ct, "A1", 1, chars(1, rot_vals(0, false), 1, 1));
  add_irrep(ct, "A2", 1, chars(1, rot_vals(0, false), -1, -1));
  if (n % 2 == 0) {
    add_irrep(ct, "B1", 1, chars(1, rot_vals(0, true), 1, -1));
    add_irrep(ct, "B2", 1, chars(1, rot_vals(0, true), -1, 1));
  }
  for (int j = 1; j <= (n - 1) / 2; ++j) {
    std::ostringstream os;
    os << "E";
    if ((n - 1) / 2 > 1) {
      os << "_" << j;
    }
    add_irrep(ct, os.str(), 2, chars(2, rot_vals(j, false), 0, 0));
  }
  return ct;
}

CharacterTable table_td() {
  CharacterTable ct;
  ct.class_labels = {"E", "8C3", "3C2", "6S4", "6σd"};
  ct.class_sizes = {1, 8, 3, 6, 6};
  add_irrep(ct, "A1", 1, {1, 1, 1, 1, 1});
  add_irrep(ct, "A2", 1, {1, 1, 1, -1, -1});
  add_irrep(ct, "E", 2, {2, -1, 2, 0, 0});
  add_irrep(ct, "T1", 3, {3, 0, -1, 1, -1});
  add_irrep(ct, "T2", 3, {3, 0, -1, -1, 1});
  return ct;
}

CharacterTable table_oh() {
  CharacterTable ct;
  ct.class_labels = {"E", "8C3", "6C2", "6C4", "3C2", "i", "6S4", "8S6", "3σh", "6σd"};
  ct.class_sizes = {1, 8, 6, 6, 3, 1, 6, 8, 3, 6};
  add_irrep(ct, "A1g", 1, {1, 1, 1, 1, 1, 1, 1, 1, 1, 1});
  add_irrep(ct, "A2g", 1, {1, 1, -1, -1, 1, 1, -1, 1, 1, -1});
  add_irrep(ct, "Eg", 2, {2, -1, 0, 0, 2, 2, 0, -1, 2, 0});
  add_irrep(ct, "T1g", 3, {3, 0, -1, 1, -1, 3, 1, 0, -1, -1});
  add_irrep(ct, "T2g", 3, {3, 0, 1, -1, -1, 3, -1, 0, -1, 1});
  add_irrep(ct, "A1u", 1, {1, 1, 1, 1, 1, -1, -1, -1, -1, -1});
  add_irrep(ct, "A2u", 1, {1, 1, -1, -1, 1, -1, 1, -1, -1, 1});
  add_irrep(ct, "Eu", 2, {2, -1, 0, 0, 2, -2, 0, 1, -2, 0});
  add_irrep(ct, "T1u", 3, {3, 0, -1, 1, -1, -3, -1, 0, 1, 1});
  add_irrep(ct, "T2u", 3, {3, 0, 1, -1, -1, -3, 1, 0, 1, -1});
  return ct;
}

CharacterTable table_t() {
  CharacterTable ct;
  ct.class_labels = {"E", "4C3", "4C3^2", "3C2"};
  ct.class_sizes = {1, 4, 4, 3};
  // Complex irreps (1E/2E); real parts of ω, ω² so dim² sums to |T|=12.
  add_irrep(ct, "A", 1, {1, 1, 1, 1});
  add_irrep(ct, "1E", 1, {1, -0.5, -0.5, 1});
  add_irrep(ct, "2E", 1, {1, -0.5, -0.5, 1});
  add_irrep(ct, "T", 3, {3, 0, 0, -1});
  return ct;
}

CharacterTable table_th() {
  CharacterTable ct;
  ct.class_labels = {"E", "4C3", "4C3^2", "3C2", "i", "4S6", "4S6^5", "3σh"};
  ct.class_sizes = {1, 4, 4, 3, 1, 4, 4, 3};
  add_irrep(ct, "Ag", 1, {1, 1, 1, 1, 1, 1, 1, 1});
  add_irrep(ct, "1Eg", 1, {1, -0.5, -0.5, 1, 1, -0.5, -0.5, 1});
  add_irrep(ct, "2Eg", 1, {1, -0.5, -0.5, 1, 1, -0.5, -0.5, 1});
  add_irrep(ct, "Tg", 3, {3, 0, 0, -1, 3, 0, 0, -1});
  add_irrep(ct, "Au", 1, {1, 1, 1, 1, -1, -1, -1, -1});
  add_irrep(ct, "1Eu", 1, {1, -0.5, -0.5, 1, -1, 0.5, 0.5, -1});
  add_irrep(ct, "2Eu", 1, {1, -0.5, -0.5, 1, -1, 0.5, 0.5, -1});
  add_irrep(ct, "Tu", 3, {3, 0, 0, -1, -3, 0, 0, 1});
  return ct;
}

CharacterTable double_irreps(const CharacterTable& base,
                             const std::string& even_suffix,
                             const std::string& odd_suffix) {
  CharacterTable ct;
  ct.class_labels = base.class_labels;
  ct.class_sizes = base.class_sizes;
  for (std::size_t i = 0; i < base.irrep_names.size(); ++i) {
    add_irrep(ct, base.irrep_names[i] + even_suffix, base.irrep_dims[i], base.characters[i]);
    add_irrep(ct, base.irrep_names[i] + odd_suffix, base.irrep_dims[i], base.characters[i]);
  }
  return ct;
}

CharacterTable table_o() {
  CharacterTable ct;
  ct.class_labels = {"E", "8C3", "6C2", "6C4", "3C2"};
  ct.class_sizes = {1, 8, 6, 6, 3};
  add_irrep(ct, "A1", 1, {1, 1, 1, 1, 1});
  add_irrep(ct, "A2", 1, {1, 1, -1, -1, 1});
  add_irrep(ct, "E", 2, {2, -1, 0, 0, 2});
  add_irrep(ct, "T1", 3, {3, 0, -1, 1, -1});
  add_irrep(ct, "T2", 3, {3, 0, 1, -1, -1});
  return ct;
}

CharacterTable table_ih() {
  CharacterTable ct;
  ct.class_labels = {"E", "12C5", "12C5^2", "20C3", "15C2", "i", "12S10", "12S10^3", "20S6", "15σ"};
  ct.class_sizes = {1, 12, 12, 20, 15, 1, 12, 12, 20, 15};
  const double phi = 0.5 * (1.0 + std::sqrt(5.0));
  add_irrep(ct, "Ag", 1, {1, 1, 1, 1, 1, 1, 1, 1, 1, 1});
  add_irrep(ct, "T1g", 3, {3, phi, 1 - phi, 0, -1, 3, 1 - phi, phi, 0, -1});
  add_irrep(ct, "T2g", 3, {3, 1 - phi, phi, 0, -1, 3, phi, 1 - phi, 0, -1});
  add_irrep(ct, "Gg", 4, {4, -1, -1, 1, 0, 4, -1, -1, 1, 0});
  add_irrep(ct, "Hg", 5, {5, 0, 0, -1, 1, 5, 0, 0, -1, 1});
  add_irrep(ct, "Au", 1, {1, 1, 1, 1, 1, -1, -1, -1, -1, -1});
  add_irrep(ct, "T1u", 3, {3, phi, 1 - phi, 0, -1, -3, phi - 1, -phi, 0, 1});
  add_irrep(ct, "T2u", 3, {3, 1 - phi, phi, 0, -1, -3, -phi, phi - 1, 0, 1});
  add_irrep(ct, "Gu", 4, {4, -1, -1, 1, 0, -4, 1, 1, -1, 0});
  add_irrep(ct, "Hu", 5, {5, 0, 0, -1, 1, -5, 0, 0, 1, -1});
  return ct;
}

CharacterTable table_i() {
  CharacterTable ct;
  ct.class_labels = {"E", "12C5", "12C5^2", "20C3", "15C2"};
  ct.class_sizes = {1, 12, 12, 20, 15};
  const double phi = 0.5 * (1.0 + std::sqrt(5.0));
  add_irrep(ct, "A", 1, {1, 1, 1, 1, 1});
  add_irrep(ct, "T1", 3, {3, phi, 1 - phi, 0, -1});
  add_irrep(ct, "T2", 3, {3, 1 - phi, phi, 0, -1});
  add_irrep(ct, "G", 4, {4, -1, -1, 1, 0});
  add_irrep(ct, "H", 5, {5, 0, 0, -1, 1});
  return ct;
}

CharacterTable table_cinfv() {
  CharacterTable ct;
  ct.class_labels = {"E", "2C∞(φ)", "∞σv"};
  ct.class_sizes = {1, 2, 0};
  add_irrep(ct, "Σ+", 1, {1, 1, 1});
  add_irrep(ct, "Σ-", 1, {1, 1, -1});
  add_irrep(ct, "Π", 2, {2, 2, 0}); // 2 cos φ stored as harmonic via class label
  add_irrep(ct, "Δ", 2, {2, 2, 0});
  add_irrep(ct, "Φ", 2, {2, 2, 0});
  return ct;
}

CharacterTable table_dinfh() {
  CharacterTable ct;
  ct.class_labels = {"E", "2C∞(φ)", "∞σv", "i", "2S∞(φ)", "∞C2"};
  ct.class_sizes = {1, 2, 0, 1, 2, 0};
  add_irrep(ct, "Σg+", 1, {1, 1, 1, 1, 1, 1});
  add_irrep(ct, "Σg-", 1, {1, 1, -1, 1, 1, -1});
  add_irrep(ct, "Πg", 2, {2, 2, 0, 2, -2, 0});
  add_irrep(ct, "Δg", 2, {2, 2, 0, 2, 2, 0});
  add_irrep(ct, "Σu+", 1, {1, 1, 1, -1, -1, -1});
  add_irrep(ct, "Σu-", 1, {1, 1, -1, -1, -1, 1});
  add_irrep(ct, "Πu", 2, {2, 2, 0, -2, 2, 0});
  add_irrep(ct, "Δu", 2, {2, 2, 0, -2, -2, 0});
  return ct;
}

std::string family(const std::string& name) {
  if (name == "C1" || name == "Cs" || name == "Ci" || name == "Cinfv" || name == "Dinfh" ||
      name == "Kh" || name == "T" || name == "Td" || name == "Th" || name == "O" || name == "Oh" ||
      name == "I" || name == "Ih") {
    return name;
  }
  if (name.size() >= 2 && name[0] == 'S') {
    return "Sn";
  }
  if (name.size() >= 2 && name[0] == 'D') {
    if (name.back() == 'h') {
      return "Dnh";
    }
    if (name.back() == 'd') {
      return "Dnd";
    }
    return "Dn";
  }
  if (name.size() >= 2 && name[0] == 'C') {
    if (name.back() == 'h') {
      return "Cnh";
    }
    if (name.back() == 'v') {
      return "Cnv";
    }
    return "Cn";
  }
  return name;
}

int rotation_power(const Operation& op) {
  if (op.power > 0) {
    return op.power % std::max(op.order, 1);
  }
  if (op.order <= 1 ||
      (op.kind != OperationKind::ProperRotation && op.kind != OperationKind::ImproperRotation)) {
    return 0;
  }
  Mat3 M = op.matrix;
  if (op.kind == OperationKind::ImproperRotation) {
    // S_n = σ · C_n along the same axis; recover the proper part.
    M = reflection_matrix(op.axis) * op.matrix;
  }
  const double tr = M.trace();
  const double c = std::clamp((tr - 1.0) * 0.5, -1.0, 1.0);
  double ang = std::acos(c);
  Vec3 skew(M(2, 1) - M(1, 2), M(0, 2) - M(2, 0), M(1, 0) - M(0, 1));
  if (op.axis.norm() > 1e-12 && skew.dot(op.axis.normalized()) < 0.0) {
    ang = -ang;
  }
  int best = 1;
  double best_err = 1e9;
  for (int k = 1; k < op.order; ++k) {
    double target = 2.0 * kPi * k / op.order;
    if (target > kPi) {
      target -= 2.0 * kPi;
    }
    const double err = std::abs(ang - target);
    if (err < best_err) {
      best_err = err;
      best = k;
    }
  }
  return best;
}

bool axis_parallel_to_c4(const Vec3& axis, const std::vector<Operation>& ops) {
  if (axis.norm() < 1e-12) {
    return false;
  }
  const Vec3 u = axis.normalized();
  for (const Operation& op : ops) {
    if (op.kind == OperationKind::ProperRotation && op.order == 4 && op.axis.norm() > 1e-12) {
      if (std::abs(u.dot(op.axis.normalized())) > 0.9) {
        return true;
      }
    }
  }
  return false;
}

int abstract_class_index(const std::string& fam,
                         const Operation& rep,
                         const std::vector<Operation>& ops) {
  if (fam == "C1") {
    return 0;
  }
  if (fam == "Cs") {
    return (rep.kind == OperationKind::Reflection) ? 1 : 0;
  }
  if (fam == "Ci") {
    return (rep.kind == OperationKind::Inversion) ? 1 : 0;
  }
  if (fam == "Td") {
    if (rep.kind == OperationKind::Identity) {
      return 0;
    }
    if (rep.kind == OperationKind::ProperRotation && rep.order == 3) {
      return 1;
    }
    if (rep.kind == OperationKind::ProperRotation && rep.order == 2) {
      return 2;
    }
    if (rep.kind == OperationKind::ImproperRotation && rep.order == 4) {
      return 3;
    }
    if (rep.kind == OperationKind::Reflection) {
      return 4;
    }
  } else if (fam == "T") {
    if (rep.kind == OperationKind::Identity) {
      return 0;
    }
    if (rep.kind == OperationKind::ProperRotation && rep.order == 3) {
      const int p = rotation_power(rep);
      return (p == 1 || p == 4) ? 1 : 2; // 4C3 vs 4C3^2
    }
    if (rep.kind == OperationKind::ProperRotation && rep.order == 2) {
      return 3;
    }
  } else if (fam == "Th") {
    if (rep.kind == OperationKind::Identity) {
      return 0;
    }
    if (rep.kind == OperationKind::ProperRotation && rep.order == 3) {
      const int p = rotation_power(rep);
      return (p == 1 || p == 4) ? 1 : 2;
    }
    if (rep.kind == OperationKind::ProperRotation && rep.order == 2) {
      return 3;
    }
    if (rep.kind == OperationKind::Inversion) {
      return 4;
    }
    if (rep.kind == OperationKind::ImproperRotation && rep.order == 6) {
      const int p = rotation_power(rep);
      return (p == 1 || p == 5) ? 5 : 6; // 4S6 vs 4S6^5
    }
    if (rep.kind == OperationKind::Reflection) {
      return 7;
    }
  } else if (fam == "O") {
    if (rep.kind == OperationKind::Identity) {
      return 0;
    }
    if (rep.kind == OperationKind::ProperRotation && rep.order == 3) {
      return 1;
    }
    if (rep.kind == OperationKind::ProperRotation && rep.order == 2) {
      return axis_parallel_to_c4(rep.axis, ops) ? 4 : 2;
    }
    if (rep.kind == OperationKind::ProperRotation && rep.order == 4) {
      return 3;
    }
  } else if (fam == "Oh") {
    if (rep.kind == OperationKind::Identity) {
      return 0;
    }
    if (rep.kind == OperationKind::ProperRotation && rep.order == 3) {
      return 1;
    }
    if (rep.kind == OperationKind::ProperRotation && rep.order == 4) {
      return 3;
    }
    if (rep.kind == OperationKind::ProperRotation && rep.order == 2) {
      return axis_parallel_to_c4(rep.axis, ops) ? 4 : 2;
    }
    if (rep.kind == OperationKind::Inversion) {
      return 5;
    }
    if (rep.kind == OperationKind::ImproperRotation && rep.order == 4) {
      return 6;
    }
    if (rep.kind == OperationKind::ImproperRotation && rep.order == 6) {
      return 7;
    }
    if (rep.kind == OperationKind::Reflection) {
      return axis_parallel_to_c4(rep.axis, ops) ? 8 : 9; // σh vs σd
    }
  } else if (fam == "I") {
    if (rep.kind == OperationKind::Identity) {
      return 0;
    }
    if (rep.kind == OperationKind::ProperRotation && rep.order == 5) {
      const int p = rotation_power(rep);
      return (p == 1 || p == 4) ? 1 : 2; // C5 vs C5^2
    }
    if (rep.kind == OperationKind::ProperRotation && rep.order == 3) {
      return 3;
    }
    if (rep.kind == OperationKind::ProperRotation && rep.order == 2) {
      return 4;
    }
  } else if (fam == "Ih") {
    if (rep.kind == OperationKind::Identity) {
      return 0;
    }
    if (rep.kind == OperationKind::ProperRotation && rep.order == 5) {
      const int p = rotation_power(rep);
      return (p == 1 || p == 4) ? 1 : 2;
    }
    if (rep.kind == OperationKind::ProperRotation && rep.order == 3) {
      return 3;
    }
    if (rep.kind == OperationKind::ProperRotation && rep.order == 2) {
      return 4;
    }
    if (rep.kind == OperationKind::Inversion) {
      return 5;
    }
    if (rep.kind == OperationKind::ImproperRotation && rep.order == 10) {
      const int p = rotation_power(rep);
      return (p == 1 || p == 9) ? 6 : 7; // S10 vs S10^3
    }
    if (rep.kind == OperationKind::ImproperRotation && rep.order == 6) {
      return 8;
    }
    if (rep.kind == OperationKind::Reflection) {
      return 9;
    }
  }
  throw Error("unmapped conjugacy class for group " + fam);
}

} // namespace

CharacterTable
build_character_table(const std::string& name, int n, const std::vector<Operation>& ops) {
  CharacterTable ct;
  const std::string fam = family(name);
  if (fam == "Cs") {
    ct = table_cs();
  } else if (fam == "Ci") {
    ct = table_ci();
  } else if (fam == "Td") {
    ct = table_td();
  } else if (fam == "Th") {
    ct = table_th();
  } else if (fam == "T") {
    ct = table_t();
  } else if (fam == "Oh") {
    ct = table_oh();
  } else if (fam == "O") {
    ct = table_o();
  } else if (fam == "Ih") {
    ct = table_ih();
  } else if (fam == "I") {
    ct = table_i();
  } else if (fam == "Cinfv") {
    return table_cinfv();
  } else if (fam == "Dinfh") {
    return table_dinfh();
  } else if (fam == "Kh") {
    ct.class_labels = {"E"};
    ct.class_sizes = {1};
    add_irrep(ct, "S", 1, {1});
    return ct;
  } else if (fam == "Cnv" || fam == "Dn") {
    ct = table_cnv(n > 0 ? n : parse_n(name));
  } else if (fam == "Dnh" || fam == "Dnd") {
    ct = double_irreps(table_cnv(n > 0 ? n : parse_n(name)), "g", "u");
  } else if (fam == "Cn" || fam == "Sn") {
    ct = table_cn(n > 0 ? n : parse_n(name));
  } else if (fam == "Cnh") {
    ct = double_irreps(table_cn(n > 0 ? n : parse_n(name)), "'", "''");
  } else {
    ct = table_c1();
  }

  if (ops.empty() || ct.characters.empty()) {
    return ct;
  }

  int nclass = 0;
  for (const Operation& op : ops) {
    nclass = std::max(nclass, op.class_index + 1);
  }
  if (nclass <= 0) {
    return ct;
  }

  CharacterTable indexed;
  indexed.irrep_names = ct.irrep_names;
  indexed.irrep_dims = ct.irrep_dims;
  indexed.class_labels.resize(nclass);
  indexed.class_sizes.assign(nclass, 0);
  for (const Operation& op : ops) {
    indexed.class_sizes[static_cast<std::size_t>(op.class_index)]++;
  }
  indexed.characters.assign(ct.irrep_names.size(), std::vector<double>(nclass, 0.0));

  const std::string& f = fam;
  const int nn = n > 0 ? n : parse_n(name);

  if (f == "C1" || f == "Cs" || f == "Ci" || f == "Td" || f == "Oh" || f == "T" || f == "Th" ||
      f == "O" || f == "I" || f == "Ih") {
    for (int ic = 0; ic < nclass; ++ic) {
      const Operation* rep = nullptr;
      for (const Operation& op : ops) {
        if (op.class_index == ic) {
          rep = &op;
          break;
        }
      }
      if (!rep) {
        throw Error("empty conjugacy class in " + fam);
      }
      const int abstract = abstract_class_index(f, *rep, ops);
      if (abstract < 0 || abstract >= static_cast<int>(ct.class_labels.size())) {
        throw Error("abstract class out of range for " + fam);
      }
      indexed.class_labels[static_cast<std::size_t>(ic)] =
          ct.class_labels[static_cast<std::size_t>(abstract)];
      for (std::size_t ir = 0; ir < ct.characters.size(); ++ir) {
        indexed.characters[ir][static_cast<std::size_t>(ic)] =
            ct.characters[ir][static_cast<std::size_t>(abstract)];
      }
    }
    return indexed;
  }

  auto strip_parity = [](std::string iname, bool& odd) {
    odd = false;
    if (iname.size() >= 2 && iname.compare(iname.size() - 2, 2, "''") == 0) {
      odd = true;
      iname.resize(iname.size() - 2);
    } else if (!iname.empty() && iname.back() == '\'') {
      iname.pop_back();
    } else if (!iname.empty() && (iname.back() == 'g' || iname.back() == 'u')) {
      odd = (iname.back() == 'u');
      iname.pop_back();
    }
    return iname;
  };

  // Axial: characters from Mulliken formulas per operation
  for (std::size_t ir = 0; ir < indexed.irrep_names.size(); ++ir) {
    bool odd_parity = false;
    const std::string base = strip_parity(indexed.irrep_names[ir], odd_parity);
    for (const Operation& op : ops) {
      double chi = 0.0;
      if (op.kind == OperationKind::Identity) {
        chi = indexed.irrep_dims[ir];
      } else if (f == "Cnv" || f == "Dn" || f == "Dnh" || f == "Dnd") {
        if (base == "A1") {
          chi = 1.0;
        } else if (base == "A2") {
          if (op.kind == OperationKind::Reflection ||
              (op.kind == OperationKind::ProperRotation && op.order == 2 &&
               std::abs(op.axis.z()) < 0.5)) {
            chi = -1.0;
          } else {
            chi = 1.0;
          }
        } else if (base == "B1") {
          if (op.kind == OperationKind::ProperRotation) {
            chi = (op.power * (nn / std::max(op.order, 1)) % 2 == 0) ? 1.0 : -1.0;
          } else if (op.kind == OperationKind::Reflection) {
            chi = (op.plane == PlaneKind::Vertical) ? 1.0 : -1.0;
          } else {
            chi = 1.0;
          }
        } else if (base == "B2") {
          if (op.kind == OperationKind::ProperRotation) {
            chi = (op.power * (nn / std::max(op.order, 1)) % 2 == 0) ? 1.0 : -1.0;
          } else if (op.kind == OperationKind::Reflection) {
            chi = (op.plane == PlaneKind::Vertical) ? -1.0 : 1.0;
          } else {
            chi = 1.0;
          }
        } else if (!base.empty() && base[0] == 'E') {
          int j = 1;
          auto pos = base.find('_');
          if (pos != std::string::npos) {
            j = std::stoi(base.substr(pos + 1));
          }
          if (op.kind == OperationKind::ProperRotation && op.order >= 2) {
            chi = 2.0 * std::cos(2.0 * kPi * j * op.power / op.order);
          } else {
            chi = 0.0;
          }
        }
      } else if (f == "Cn" || f == "Cnh" || f == "Sn") {
        if (base == "A") {
          chi = 1.0;
        } else if (base == "B") {
          chi = (op.kind == OperationKind::ProperRotation && (op.power % 2)) ? -1.0 : 1.0;
        } else if (!base.empty() && base[0] == 'E') {
          int j = 1;
          auto pos = base.find('_');
          if (pos != std::string::npos) {
            j = std::stoi(base.substr(pos + 1));
          }
          if (op.kind == OperationKind::ProperRotation) {
            chi = 2.0 * std::cos(2.0 * kPi * j * op.power / op.order);
          } else {
            chi = 0.0;
          }
        }
      }
      const bool improper = op.kind == OperationKind::Inversion ||
                            op.kind == OperationKind::Reflection ||
                            op.kind == OperationKind::ImproperRotation;
      if (odd_parity && improper) {
        chi = -chi;
      }
      indexed.characters[ir][static_cast<std::size_t>(op.class_index)] = chi;
    }
  }
  for (int ic = 0; ic < nclass; ++ic) {
    for (const Operation& op : ops) {
      if (op.class_index == ic) {
        std::ostringstream os;
        switch (op.kind) {
        case OperationKind::Identity:
          os << "E";
          break;
        case OperationKind::Inversion:
          os << "i";
          break;
        case OperationKind::ProperRotation:
          os << "C" << op.order;
          break;
        case OperationKind::ImproperRotation:
          os << "S" << op.order;
          break;
        case OperationKind::Reflection:
          os << "σ";
          break;
        }
        indexed.class_labels[static_cast<std::size_t>(ic)] = os.str();
        break;
      }
    }
  }
  return indexed;
}

} // namespace molsymm::detail
