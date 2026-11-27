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

CharacterTable table_d2d() {
  CharacterTable ct;
  ct.class_labels = {"E", "2S4", "C2", "2C2'", "2σd"};
  ct.class_sizes = {1, 2, 1, 2, 2};
  add_irrep(ct, "A1", 1, {1, 1, 1, 1, 1});
  add_irrep(ct, "A2", 1, {1, 1, 1, -1, -1});
  add_irrep(ct, "B1", 1, {1, -1, 1, 1, -1});
  add_irrep(ct, "B2", 1, {1, -1, 1, -1, 1});
  add_irrep(ct, "E", 2, {2, 0, -2, 0, 0});
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
  // Real Mulliken E = 1E ⊕ 2E (same real-form compromise as Cn).
  add_irrep(ct, "A", 1, {1, 1, 1, 1});
  add_irrep(ct, "E", 2, {2, -1, -1, 2});
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

int rotation_power(const Operation& op);

// C2'/σv vs C2''/σd for even-n axial groups by azimuth of the axis/normal.
bool is_prime_in_plane(const Vec3& axis, int n) {
  if (axis.norm() < 1e-12 || n < 2) {
    return true;
  }
  const Vec3 u = axis.normalized();
  // Mulliken C2v (molecule in yz): σ(xz) with normal ‖ y is the prime vertical plane,
  // so out-of-plane x transforms as B1.
  if (n == 2) {
    return std::abs(u.y()) >= std::abs(u.x());
  }
  const double phi = std::atan2(u.y(), u.x());
  const double period = kPi / static_cast<double>(n);
  double t = std::fmod(std::abs(phi), 2.0 * period);
  if (t > period) {
    t = 2.0 * period - t;
  }
  return t <= 0.5 * period + 1e-8;
}

int e_index(const std::string& base) {
  auto pos = base.find('_');
  if (pos == std::string::npos) {
    return 1;
  }
  return std::stoi(base.substr(pos + 1));
}

// Character on a proper rotation about z. Remap p into principal_n for B parity when
// C_n^k is stored as a lower-order axis (e.g. C4^2 as C2).
double axial_principal_chi(const std::string& base, int p, int order, int principal_n) {
  const int ord = std::max(order, 1);
  const int k = (principal_n > 0 && principal_n % ord == 0) ? p * (principal_n / ord) : p;
  if (base == "A1" || base == "A" || base == "A2") {
    return 1.0;
  }
  if (base == "B1" || base == "B2" || base == "B") {
    return (k % 2 == 0) ? 1.0 : -1.0;
  }
  if (!base.empty() && base[0] == 'E') {
    return 2.0 * std::cos(2.0 * kPi * e_index(base) * p / ord);
  }
  return 0.0;
}

// C2'/C2'' and σv/σd characters share the same A/B pattern (prime ↔ vertical).
double axial_perp_chi(const std::string& base, bool prime) {
  if (base == "A1" || base == "A") {
    return 1.0;
  }
  if (base == "A2") {
    return -1.0;
  }
  if (base == "B1") {
    return prime ? 1.0 : -1.0;
  }
  if (base == "B2") {
    return prime ? -1.0 : 1.0;
  }
  if (base == "B") {
    return -1.0;
  }
  return 0.0;
}

bool is_e_irrep(const std::string& base) {
  return !base.empty() && base[0] == 'E';
}

double axial_base_dim(const std::string& base) {
  return is_e_irrep(base) ? 2.0 : 1.0;
}

double axial_base_chi(const std::string& base, const Operation& op, int n) {
  if (op.kind == OperationKind::Identity) {
    return axial_base_dim(base);
  }
  if (op.kind == OperationKind::ProperRotation) {
    if (axis_on_z(op) || op.order > 2) {
      return axial_principal_chi(base, rotation_power(op), op.order, n);
    }
    return axial_perp_chi(base, is_prime_in_plane(op.axis, n));
  }
  if (op.kind == OperationKind::Reflection) {
    if (op.plane == PlaneKind::Horizontal || axis_on_z(op)) {
      if (is_e_irrep(base)) {
        return 2.0 * std::cos(kPi * static_cast<double>(e_index(base)));
      }
      return axial_principal_chi(base, 1, 2, n);
    }
    const bool sigma_v =
        (n % 2 == 0) ? is_prime_in_plane(op.axis, n) : (op.plane != PlaneKind::Dihedral);
    return axial_perp_chi(base, sigma_v);
  }
  return 0.0;
}

// Base character on an improper element, before g/u or '/'/'' flip.
double axial_improper_base_chi(
    const std::string& base, const Operation& op, int n, bool via_sigma_h, bool dnd_ci) {
  if (op.kind == OperationKind::Inversion) {
    return axial_base_dim(base);
  }
  if (op.kind == OperationKind::Reflection) {
    if (op.plane == PlaneKind::Horizontal || axis_on_z(op)) {
      if (is_e_irrep(base)) {
        // Odd Dnh/Cnh (σh doubling): E' has +2 on σh. Even (Ci): σh ↔ C2(z).
        if (via_sigma_h) {
          return 2.0;
        }
        return 2.0 * std::cos(kPi * static_cast<double>(e_index(base)));
      }
      return axial_principal_chi(base, 1, 2, n);
    }
    // Ci doubling: σv/σd ↔ C2 along the normal.
    return axial_perp_chi(base, is_prime_in_plane(op.axis, n));
  }
  if (op.kind == OperationKind::ImproperRotation) {
    // Inversion doubling (even Dnh/Cnh, odd Dnd): χ_g(S(θ)) = χ(C(θ+π)) via −S.
    // σh doubling (odd Dnh/Cnh): S = σh·C → χ(S) = χ(C) via rotation_power.
    if (dnd_ci || !via_sigma_h) {
      int ord = 0;
      int pow = 0;
      Vec3 ax = Vec3::UnitZ();
      if (classify_proper_rotation(-op.matrix, ord, pow, ax) && ord >= 2) {
        return axial_principal_chi(base, pow, ord, n);
      }
    }
    return axial_principal_chi(base, rotation_power(op), op.order, n);
  }
  return 0.0;
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
  // Degenerate rows are the φ = 0 snapshot used only for labeling; reduction is geometric.
  add_irrep(ct, "Π", 2, {2, 2, 0});
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
  if (op.order <= 1 ||
      (op.kind != OperationKind::ProperRotation && op.kind != OperationKind::ImproperRotation)) {
    return 0;
  }
  // Proper: trust stored power. Improper: recompute from σ·S (closure often leaves power=1).
  if (op.kind == OperationKind::ProperRotation && op.power > 0) {
    return op.power % std::max(op.order, 1);
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
    // Face-axis C2 may be stored as C4^2 (order 4) when Td is taken from O/Oh.
    if (rep.kind == OperationKind::ProperRotation &&
        (rep.order == 2 || rotation_power(rep) == 2 || std::abs(rep.matrix.trace() + 1.0) < 1e-4)) {
      return 2;
    }
    if (rep.kind == OperationKind::ImproperRotation && rep.order == 4) {
      return 3;
    }
    if (rep.kind == OperationKind::Reflection) {
      return 4;
    }
  } else if (fam == "Dnd") {
    if (rep.kind == OperationKind::Identity) {
      return 0;
    }
    if (rep.kind == OperationKind::ImproperRotation && rep.order == 4) {
      return 1;
    }
    if (rep.kind == OperationKind::ProperRotation) {
      const bool on_z = rep.axis.norm() > 1e-12 && std::abs(rep.axis.normalized().z()) > 0.9;
      // 180° rotations may be stored as C4^2 (order 4) when taken from O/Oh.
      const bool is_180 =
          rep.order == 2 || rotation_power(rep) == 2 || std::abs(rep.matrix.trace() + 1.0) < 1e-4;
      if (is_180) {
        return on_z ? 2 : 3; // C2 vs 2C2'
      }
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
    // Face-axis C2 may be stored as C4^2 (order 4) when T is taken from O/Oh.
    if (rep.kind == OperationKind::ProperRotation &&
        (rep.order == 2 || rotation_power(rep) == 2 || std::abs(rep.matrix.trace() + 1.0) < 1e-4)) {
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
    if (rep.kind == OperationKind::ProperRotation &&
        (rep.order == 2 || rotation_power(rep) == 2 || std::abs(rep.matrix.trace() + 1.0) < 1e-4)) {
      return 3;
    }
    if (rep.kind == OperationKind::Inversion) {
      return 4;
    }
    if (rep.kind == OperationKind::ImproperRotation && rep.order == 6) {
      const int p = rotation_power(rep);
      return (p == 1) ? 5 : 6; // 4S6 vs 4S6^5
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
    if (rep.kind == OperationKind::ProperRotation && rep.order == 4) {
      // C4^2 belongs to the face-axis 3C2 class, not 6C4.
      if (rotation_power(rep) == 2 || std::abs(rep.matrix.trace() + 1.0) < 1e-4) {
        return 4;
      }
      return 3;
    }
    if (rep.kind == OperationKind::ProperRotation && rep.order == 2) {
      return axis_parallel_to_c4(rep.axis, ops) ? 4 : 2;
    }
  } else if (fam == "Oh") {
    if (rep.kind == OperationKind::Identity) {
      return 0;
    }
    if (rep.kind == OperationKind::ProperRotation && rep.order == 3) {
      return 1;
    }
    if (rep.kind == OperationKind::ProperRotation && rep.order == 4) {
      if (rotation_power(rep) == 2 || std::abs(rep.matrix.trace() + 1.0) < 1e-4) {
        return 4; // 3C2
      }
      return 3; // 6C4
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
  const int nn = n > 0 ? n : parse_n(name);
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
    ct = table_cnv(nn);
  } else if (fam == "Dnd") {
    if (nn % 2 == 1) {
      ct = double_irreps(table_cnv(nn), "g", "u"); // inversion → g/u
    } else if (nn == 2) {
      ct = table_d2d();
    } else {
      ct = table_cnv(2 * nn); // even Dnd ≅ C_{2n}v
    }
  } else if (fam == "Dnh") {
    if (nn % 2 == 0) {
      ct = double_irreps(table_cnv(nn), "g", "u");
    } else {
      ct = double_irreps(table_cnv(nn), "'", "''"); // σh only
    }
  } else if (fam == "Cn" || fam == "Sn") {
    ct = table_cn(nn);
  } else if (fam == "Cnh") {
    if (nn % 2 == 0) {
      ct = double_irreps(table_cn(nn), "g", "u");
    } else {
      ct = double_irreps(table_cn(nn), "'", "''");
    }
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

  std::vector<int> rep_idx(static_cast<std::size_t>(nclass), -1);
  CharacterTable indexed;
  indexed.irrep_names = ct.irrep_names;
  indexed.irrep_dims = ct.irrep_dims;
  indexed.class_labels.resize(nclass);
  indexed.class_sizes.assign(nclass, 0);
  for (int i = 0; i < static_cast<int>(ops.size()); ++i) {
    const int ic = ops[static_cast<std::size_t>(i)].class_index;
    indexed.class_sizes[static_cast<std::size_t>(ic)]++;
    if (rep_idx[static_cast<std::size_t>(ic)] < 0) {
      rep_idx[static_cast<std::size_t>(ic)] = i;
    }
  }
  indexed.characters.assign(ct.irrep_names.size(), std::vector<double>(nclass, 0.0));

  const std::string& f = fam;

  if (f == "C1" || f == "Cs" || f == "Ci" || f == "Td" || f == "Oh" || f == "T" || f == "Th" ||
      f == "O" || f == "I" || f == "Ih" || (f == "Dnd" && nn == 2)) {
    for (int ic = 0; ic < nclass; ++ic) {
      const int ri = rep_idx[static_cast<std::size_t>(ic)];
      if (ri < 0) {
        throw Error("empty conjugacy class in " + fam);
      }
      const Operation& rep = ops[static_cast<std::size_t>(ri)];
      const int abstract = abstract_class_index(f, rep, ops);
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
    if (iname.ends_with("''")) {
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

  const bool has_sigma_h_doubling = (f == "Dnh" && nn % 2 == 1) || (f == "Cnh" && nn % 2 == 1);
  const bool has_parity_doubling = has_sigma_h_doubling || (f == "Dnh" && nn % 2 == 0) ||
                                   (f == "Cnh" && nn % 2 == 0) || (f == "Dnd" && nn % 2 == 1);
  const bool is_dihedral_family = (f == "Cnv" || f == "Dn" || f == "Dnh" || f == "Dnd");
  const bool is_cyclic_family = (f == "Cn" || f == "Cnh" || f == "Sn");
  const bool dnd_ci = (f == "Dnd" && nn % 2 == 1);
  const bool even_dnd = (f == "Dnd" && nn > 2 && nn % 2 == 0);
  const int n_char = even_dnd ? (2 * nn) : nn;

  for (std::size_t ir = 0; ir < indexed.irrep_names.size(); ++ir) {
    bool odd_parity = false;
    const std::string base = strip_parity(indexed.irrep_names[ir], odd_parity);
    for (int ic = 0; ic < nclass; ++ic) {
      const int ri = rep_idx[static_cast<std::size_t>(ic)];
      if (ri < 0) {
        continue;
      }
      const Operation& op = ops[static_cast<std::size_t>(ri)];
      double chi = 0.0;

      if (is_cyclic_family) {
        if (op.kind == OperationKind::Identity) {
          chi = indexed.irrep_dims[ir];
        } else if (op.kind == OperationKind::ProperRotation ||
                   (op.kind == OperationKind::ImproperRotation && f == "Sn")) {
          chi = axial_principal_chi(base, rotation_power(op), op.order, nn);
        } else if (op.kind == OperationKind::Inversion && f == "Sn" && nn % 2 == 0) {
          // S_n ≅ C_n: i = S_n^{n/2} has the character of C_n^{n/2}.
          chi = axial_principal_chi(base, nn / 2, nn, nn);
        } else if (op.kind == OperationKind::Inversion || op.kind == OperationKind::Reflection ||
                   op.kind == OperationKind::ImproperRotation) {
          if (!has_parity_doubling) {
            throw Error("unmapped conjugacy class for group " + fam);
          }
          chi = axial_improper_base_chi(base, op, nn, has_sigma_h_doubling, false);
          if (odd_parity) {
            chi = -chi;
          }
        } else {
          throw Error("unmapped conjugacy class for group " + fam);
        }
      } else if (is_dihedral_family) {
        const bool improper = op.kind == OperationKind::Inversion ||
                              op.kind == OperationKind::Reflection ||
                              op.kind == OperationKind::ImproperRotation;
        if (even_dnd) {
          // Even Dnd ≅ C_{2n}v: S_{2n} ↔ C_{2n}, so C_n = C_{2n}^{2} uses scale 2n/order.
          if (op.kind == OperationKind::ImproperRotation) {
            chi = axial_principal_chi(base, rotation_power(op), op.order, n_char);
          } else if (op.kind == OperationKind::ProperRotation && axis_on_z(op)) {
            chi = axial_principal_chi(
                base, rotation_power(op) * (n_char / std::max(op.order, 1)), n_char, n_char);
          } else {
            chi = axial_base_chi(base, op, n_char);
          }
        } else if (improper && has_parity_doubling) {
          chi = axial_improper_base_chi(base, op, n_char, has_sigma_h_doubling, dnd_ci);
          if (odd_parity) {
            chi = -chi;
          }
        } else {
          chi = axial_base_chi(base, op, n_char);
        }
      } else {
        throw Error("unmapped conjugacy class for group " + fam);
      }
      indexed.characters[ir][static_cast<std::size_t>(ic)] = chi;
    }
  }
  for (int ic = 0; ic < nclass; ++ic) {
    const int ri = rep_idx[static_cast<std::size_t>(ic)];
    if (ri < 0) {
      continue;
    }
    const Operation& op = ops[static_cast<std::size_t>(ri)];
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
  }
  return indexed;
}

} // namespace molsymm::detail
