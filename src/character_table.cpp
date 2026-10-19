#include "internal.hpp"

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
    std::ostringstream os;
    os << "C_2";
    ct.class_labels.push_back(os.str());
    ct.class_sizes.push_back(1);
  }
  if (n % 2 == 0) {
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
    std::vector<double> rv;
    for (int k = 1; k <= (n - 1) / 2; ++k) {
      rv.push_back(2.0 * std::cos(2.0 * kPi * j * k / n));
    }
    if (n % 2 == 0) {
      rv.push_back(2.0 * std::cos(kPi * j));
    }
    std::ostringstream os;
    os << "E";
    if ((n - 1) / 2 > 1) {
      os << "_" << j;
    }
    add_irrep(ct, os.str(), 2, chars(2, rv, 0, 0));
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
  add_irrep(ct, "A", 1, {1, 1, 1, 1});
  add_irrep(ct, "E", 2, {2, -1, -1, 2}); // real reducible form of complex pair
  add_irrep(ct, "T", 3, {3, 0, 0, -1});
  return ct;
}

CharacterTable table_th() {
  CharacterTable ct;
  ct.class_labels = {"E", "4C3", "4C3^2", "3C2", "i", "4S6", "4S6^5", "3σh"};
  ct.class_sizes = {1, 4, 4, 3, 1, 4, 4, 3};
  add_irrep(ct, "Ag", 1, {1, 1, 1, 1, 1, 1, 1, 1});
  add_irrep(ct, "Eg", 2, {2, -1, -1, 2, 2, -1, -1, 2});
  add_irrep(ct, "Tg", 3, {3, 0, 0, -1, 3, 0, 0, -1});
  add_irrep(ct, "Au", 1, {1, 1, 1, 1, -1, -1, -1, -1});
  add_irrep(ct, "Eu", 2, {2, -1, -1, 2, -2, 1, 1, -2});
  add_irrep(ct, "Tu", 3, {3, 0, 0, -1, -3, 0, 0, 1});
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

int parse_n(const std::string& name) {
  std::size_t i = 0;
  while (i < name.size() && !std::isdigit(static_cast<unsigned char>(name[i]))) {
    ++i;
  }
  if (i >= name.size()) {
    return 0;
  }
  int n = 0;
  while (i < name.size() && std::isdigit(static_cast<unsigned char>(name[i]))) {
    n = n * 10 + (name[i] - '0');
    ++i;
  }
  return n;
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
    ct = table_cinfv();
  } else if (fam == "Dinfh") {
    ct = table_dinfh();
  } else if (fam == "Kh") {
    ct.class_labels = {"E"};
    ct.class_sizes = {1};
    add_irrep(ct, "S", 1, {1});
  } else if (fam == "Cnv" || fam == "Dn" || fam == "Dnh" || fam == "Dnd") {
    // Use Cnv abelianization for axial groups; χ(g) filled per operation below.
    ct = table_cnv(n > 0 ? n : parse_n(name));
  } else if (fam == "Cn" || fam == "Cnh" || fam == "Sn") {
    ct = table_cn(n > 0 ? n : parse_n(name));
  } else {
    // C1 and unknown groups
    ct = table_c1();
  }

  // Build an operation-indexed character lookup: replace table with one column per class_index
  if (!ops.empty() && !ct.characters.empty()) {
    int nclass = 0;
    for (const Operation& op : ops) {
      nclass = std::max(nclass, op.class_index + 1);
    }
    if (nclass > 0) {
      CharacterTable indexed;
      indexed.irrep_names = ct.irrep_names;
      indexed.irrep_dims = ct.irrep_dims;
      indexed.class_labels.resize(nclass);
      indexed.class_sizes.assign(nclass, 0);
      for (const Operation& op : ops) {
        indexed.class_sizes[static_cast<std::size_t>(op.class_index)]++;
      }
      // Character of each class from matrix trace in representation is hard;
      // use geometric formula: χ for vector rep etc. For Mulliken labels we compute
      // from operation kind for standard 1D/2D axial groups.
      indexed.characters.assign(ct.irrep_names.size(), std::vector<double>(nclass, 0.0));

      // Fill by matching each operation to the nearest abstract class via character of
      // the defining representation (3D orthogonal matrix trace & det).
      // For axial groups, assign using order and plane kind.
      const std::string& f = fam;
      const int nn = n > 0 ? n : parse_n(name);

      if (f == "C1" || f == "Cs" || f == "Ci" || f == "Td" || f == "Oh" || f == "T" || f == "Th" ||
          f == "O" || f == "I" || f == "Ih") {
        // Map each op to abstract class by (kind, order) heuristic
        for (int ic = 0; ic < nclass; ++ic) {
          // Find a representative operation
          const Operation* rep = nullptr;
          for (const Operation& op : ops) {
            if (op.class_index == ic) {
              rep = &op;
              break;
            }
          }
          if (!rep) {
            continue;
          }
          int abstract = 0;
          if (f == "Td") {
            if (rep->kind == OperationKind::Identity) {
              abstract = 0;
            } else if (rep->kind == OperationKind::ProperRotation && rep->order == 3) {
              abstract = 1;
            } else if (rep->kind == OperationKind::ProperRotation && rep->order == 2) {
              abstract = 2;
            } else if (rep->kind == OperationKind::ImproperRotation && rep->order == 4) {
              abstract = 3;
            } else if (rep->kind == OperationKind::Reflection) {
              abstract = 4;
            }
          } else if (f == "Oh") {
            if (rep->kind == OperationKind::Identity) {
              abstract = 0;
            } else if (rep->kind == OperationKind::ProperRotation && rep->order == 3) {
              abstract = 1;
            } else if (rep->kind == OperationKind::ProperRotation && rep->order == 2 &&
                       std::abs(rep->axis.normalized().dot(Vec3::UnitZ())) < 0.5 &&
                       std::abs(rep->axis.normalized().dot(Vec3::UnitX())) < 0.5) {
              abstract = 2; // C2 along face diagonal-ish — approximate
            } else if (rep->kind == OperationKind::ProperRotation && rep->order == 4) {
              abstract = 3;
            } else if (rep->kind == OperationKind::ProperRotation && rep->order == 2) {
              abstract = 4;
            } else if (rep->kind == OperationKind::Inversion) {
              abstract = 5;
            } else if (rep->kind == OperationKind::ImproperRotation && rep->order == 4) {
              abstract = 6;
            } else if (rep->kind == OperationKind::ImproperRotation && rep->order == 6) {
              abstract = 7;
            } else if (rep->kind == OperationKind::Reflection &&
                       rep->plane == PlaneKind::Horizontal) {
              abstract = 8;
            } else if (rep->kind == OperationKind::Reflection) {
              abstract = 9;
            }
          } else if (f == "C1") {
            abstract = 0;
          } else if (f == "Cs") {
            abstract = (rep->kind == OperationKind::Reflection) ? 1 : 0;
          } else if (f == "Ci") {
            abstract = (rep->kind == OperationKind::Inversion) ? 1 : 0;
          }
          abstract = std::min(abstract, static_cast<int>(ct.class_labels.size()) - 1);
          indexed.class_labels[static_cast<std::size_t>(ic)] =
              ct.class_labels[static_cast<std::size_t>(abstract)];
          for (std::size_t ir = 0; ir < ct.characters.size(); ++ir) {
            indexed.characters[ir][static_cast<std::size_t>(ic)] =
                ct.characters[ir][static_cast<std::size_t>(abstract)];
          }
        }
        return indexed;
      }

      // Axial: compute characters from formulas per operation
      for (std::size_t ir = 0; ir < indexed.irrep_names.size(); ++ir) {
        const std::string& iname = indexed.irrep_names[ir];
        for (const Operation& op : ops) {
          double chi = 0.0;
          if (op.kind == OperationKind::Identity) {
            chi = indexed.irrep_dims[ir];
          } else if (f == "Cnv" || f == "Dn" || f == "Dnh" || f == "Dnd") {
            // Simplified: A1 all +1; A2 σ=-1; E: 2cos on Cn, 0 on σ
            if (iname == "A1" || iname == "A1g") {
              chi = 1.0;
            } else if (iname == "A2" || iname == "A2g") {
              if (op.kind == OperationKind::Reflection ||
                  (op.kind == OperationKind::ProperRotation && op.order == 2 &&
                   std::abs(op.axis.z()) < 0.5)) {
                chi = -1.0;
              } else {
                chi = 1.0;
              }
            } else if (iname == "B1" || iname == "B1g") {
              if (op.kind == OperationKind::ProperRotation) {
                chi = (op.power * (nn / std::max(op.order, 1)) % 2 == 0) ? 1.0 : -1.0;
              } else if (op.kind == OperationKind::Reflection) {
                chi = (op.plane == PlaneKind::Vertical) ? 1.0 : -1.0;
              } else {
                chi = 1.0;
              }
            } else if (iname == "B2" || iname == "B2g") {
              if (op.kind == OperationKind::ProperRotation) {
                chi = (op.power * (nn / std::max(op.order, 1)) % 2 == 0) ? 1.0 : -1.0;
              } else if (op.kind == OperationKind::Reflection) {
                chi = (op.plane == PlaneKind::Vertical) ? -1.0 : 1.0;
              } else {
                chi = 1.0;
              }
            } else if (!iname.empty() && iname[0] == 'E') {
              int j = 1;
              auto pos = iname.find('_');
              if (pos != std::string::npos) {
                j = std::stoi(iname.substr(pos + 1));
              }
              if (op.kind == OperationKind::ProperRotation && op.order >= 2) {
                chi = 2.0 * std::cos(2.0 * kPi * j * op.power / op.order);
              } else {
                chi = 0.0;
              }
            } else {
              chi = (op.kind == OperationKind::Identity) ? indexed.irrep_dims[ir] : 0.0;
            }
            if (op.kind == OperationKind::Inversion) {
              if (!iname.empty() && iname.back() == 'u') {
                chi = -chi;
              }
            }
          } else if (f == "Cn" || f == "Cnh" || f == "Sn") {
            if (iname == "A") {
              chi = 1.0;
            } else if (iname == "B") {
              chi = (op.kind == OperationKind::ProperRotation && (op.power % 2)) ? -1.0 : 1.0;
            } else if (!iname.empty() && iname[0] == 'E') {
              int j = 1;
              auto pos = iname.find('_');
              if (pos != std::string::npos) {
                j = std::stoi(iname.substr(pos + 1));
              }
              if (op.kind == OperationKind::ProperRotation) {
                chi = 2.0 * std::cos(2.0 * kPi * j * op.power / op.order);
              } else {
                chi = 0.0;
              }
            }
          }
          indexed.characters[ir][static_cast<std::size_t>(op.class_index)] = chi;
        }
      }
      // Labels
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
  }

  return ct;
}

} // namespace molsymm::detail
