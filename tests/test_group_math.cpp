#include "molsymm/molsymm.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

int g_failed = 0;
int g_passed = 0;

void expect(bool cond, const std::string& msg) {
  if (!cond) {
    std::cerr << "FAIL: " << msg << "\n";
    ++g_failed;
  } else {
    ++g_passed;
  }
}

bool has_subgroup(const molsymm::Analysis& a, const std::string& name) {
  return std::find(a.subgroups.begin(), a.subgroups.end(), name) != a.subgroups.end();
}

int irrep_index(const molsymm::CharacterTable& ct, const std::string& name) {
  for (std::size_t i = 0; i < ct.irrep_names.size(); ++i) {
    if (ct.irrep_names[i] == name) {
      return static_cast<int>(i);
    }
  }
  // Tables use E_1 / E_1g; callers may write E1 / E1g.
  if (name.size() >= 2 && name[0] == 'E' && name[1] != '_') {
    const std::string alt = "E_" + name.substr(1);
    for (std::size_t i = 0; i < ct.irrep_names.size(); ++i) {
      if (ct.irrep_names[i] == alt) {
        return static_cast<int>(i);
      }
    }
  }
  return -1;
}

double chi_on(const molsymm::Analysis& a, const std::string& irrep, const molsymm::Operation& op) {
  const int ir = irrep_index(a.character_table, irrep);
  if (ir < 0 || op.class_index < 0 ||
      op.class_index >=
          static_cast<int>(a.character_table.characters[static_cast<std::size_t>(ir)].size())) {
    return 1e9;
  }
  return a.character_table
      .characters[static_cast<std::size_t>(ir)][static_cast<std::size_t>(op.class_index)];
}

molsymm::Molecule water() {
  std::vector<int> z{8, 1, 1};
  std::vector<double> xyz{0.0, 0.0, 0.1173, 0.0, 0.7572, -0.4692, 0.0, -0.7572, -0.4692};
  return molsymm::Molecule(z, xyz);
}

molsymm::Molecule methane() {
  const double a = 0.627;
  std::vector<int> z{6, 1, 1, 1, 1};
  std::vector<double> xyz{0, 0, 0, a, a, a, a, -a, -a, -a, a, -a, -a, -a, a};
  return molsymm::Molecule(z, xyz);
}

molsymm::Molecule benzene() {
  std::vector<int> z(12, 6);
  for (int i = 6; i < 12; ++i) {
    z[static_cast<std::size_t>(i)] = 1;
  }
  std::vector<double> xyz;
  constexpr double kPi = 3.14159265358979323846;
  for (int i = 0; i < 6; ++i) {
    const double ang = 2.0 * kPi * i / 6.0;
    xyz.push_back(1.39 * std::cos(ang));
    xyz.push_back(1.39 * std::sin(ang));
    xyz.push_back(0.0);
  }
  for (int i = 0; i < 6; ++i) {
    const double ang = 2.0 * kPi * i / 6.0;
    xyz.push_back(2.47 * std::cos(ang));
    xyz.push_back(2.47 * std::sin(ang));
    xyz.push_back(0.0);
  }
  return molsymm::Molecule(z, xyz);
}

molsymm::Molecule sf6() {
  std::vector<int> z{16, 9, 9, 9, 9, 9, 9};
  const double r = 1.56;
  std::vector<double> xyz{0, 0, 0, r, 0, 0, -r, 0, 0, 0, r, 0, 0, -r, 0, 0, 0, r, 0, 0, -r};
  return molsymm::Molecule(z, xyz);
}

molsymm::Molecule co2() {
  std::vector<int> z{6, 8, 8};
  std::vector<double> xyz{0, 0, 0, 0, 0, 1.16, 0, 0, -1.16};
  return molsymm::Molecule(z, xyz);
}

molsymm::Molecule staggered_ethane() {
  std::vector<int> z{1, 6, 1, 1, 6, 1, 1, 1};
  std::vector<double> xyz{1.1851,  -0.0039, 0.9875,  0.7516,  -0.0225, -0.0209, 1.1669, 0.8330,
                          -0.5693, 1.1155,  -0.9329, -0.5145, -0.7516, 0.0225,  0.0209, -1.1669,
                          -0.8334, 0.5687,  -1.1157, 0.9326,  0.5151,  -1.1850, 0.0044, -0.9875};
  return molsymm::Molecule(z, xyz);
}

molsymm::Molecule square_antiprism() {
  constexpr double kPi = 3.14159265358979323846;
  constexpr double r = 1.0;
  constexpr double z = 0.7;
  std::vector<int> atomic(8, 6);
  std::vector<double> xyz;
  for (int i = 0; i < 4; ++i) {
    const double ang = 0.5 * kPi * i;
    xyz.push_back(r * std::cos(ang));
    xyz.push_back(r * std::sin(ang));
    xyz.push_back(-z);
  }
  for (int i = 0; i < 4; ++i) {
    const double ang = 0.5 * kPi * i + 0.25 * kPi;
    xyz.push_back(r * std::cos(ang));
    xyz.push_back(r * std::sin(ang));
    xyz.push_back(z);
  }
  return molsymm::Molecule(atomic, xyz);
}

void expect_orthogonal(const molsymm::Analysis& a, const std::string& label) {
  const auto& ct = a.character_table;
  if (ct.characters.empty() || ct.class_sizes.size() != ct.characters.front().size()) {
    expect(false, label + " missing class sizes");
    return;
  }
  const int nir = static_cast<int>(ct.irrep_names.size());
  const int ncl = static_cast<int>(ct.class_sizes.size());
  for (int i = 0; i < nir; ++i) {
    for (int j = i; j < nir; ++j) {
      double sum = 0.0;
      for (int c = 0; c < ncl; ++c) {
        sum += static_cast<double>(ct.class_sizes[static_cast<std::size_t>(c)]) *
               ct.characters[static_cast<std::size_t>(i)][static_cast<std::size_t>(c)] *
               ct.characters[static_cast<std::size_t>(j)][static_cast<std::size_t>(c)];
      }
      const double want = (i == j) ? static_cast<double>(a.order) : 0.0;
      expect(std::abs(sum - want) < 1e-6,
             label + " orthogonality " + ct.irrep_names[static_cast<std::size_t>(i)] + "/" +
                 ct.irrep_names[static_cast<std::size_t>(j)] + " got " + std::to_string(sum));
    }
  }
}

void test_benzene_e1g_s_columns() {
  auto a = molsymm::analyze(benzene());
  expect(a.point_group == "D6h", "benzene D6h");
  expect(!has_subgroup(a, "D1") && !has_subgroup(a, "D1h") && !has_subgroup(a, "C1v"),
         "benzene has no phantom D1/D1h/C1v");
  bool saw_s6 = false;
  bool saw_s3 = false;
  for (const auto& op : a.operations) {
    if (op.kind != molsymm::OperationKind::ImproperRotation) {
      continue;
    }
    const double tr = op.matrix.trace();
    if (std::abs(tr) < 1e-4) {
      expect(std::abs(chi_on(a, "E1g", op) + 1.0) < 1e-6, "benzene E1g(S6) = −1");
      saw_s6 = true;
    } else if (std::abs(tr + 2.0) < 1e-4) {
      expect(std::abs(chi_on(a, "E1g", op) - 1.0) < 1e-6, "benzene E1g(S3) = +1");
      saw_s3 = true;
    }
  }
  expect(saw_s6 && saw_s3, "benzene has S6 and S3 classes");

  bool saw_c2p = false;
  for (const auto& op : a.operations) {
    if (op.kind != molsymm::OperationKind::ProperRotation || op.axis.norm() < 1e-12) {
      continue;
    }
    const bool half_turn = op.order == 2 || (op.order == 6 && op.power == 3) ||
                           std::abs(op.matrix.trace() + 1.0) < 1e-4;
    if (!half_turn) {
      continue;
    }
    const auto u = op.axis.normalized();
    if (std::abs(u.z()) > 0.2) {
      continue;
    }
    const double phi = std::atan2(u.y(), u.x());
    constexpr double kPi = 3.14159265358979323846;
    const double step = kPi / 3.0;
    double t = std::fmod(std::abs(phi), step);
    if (t > 0.5 * step) {
      t = step - t;
    }
    if (t < 0.15) {
      expect(std::abs(chi_on(a, "B1g", op) - 1.0) < 1e-6, "benzene B1g(C2') = +1");
      saw_c2p = true;
      break;
    }
  }
  expect(saw_c2p, "benzene has in-plane C2' through atoms");
}

void test_d4d_characters() {
  auto a = molsymm::analyze(square_antiprism());
  expect(a.point_group == "D4d", "square antiprism D4d got " + a.point_group);
  expect_orthogonal(a, "D4d");
  bool saw_c4 = false;
  bool saw_c2 = false;
  for (const auto& op : a.operations) {
    if (op.kind != molsymm::OperationKind::ProperRotation || op.axis.norm() < 1e-12) {
      continue;
    }
    if (std::abs(op.axis.normalized().z()) < 0.9) {
      continue;
    }
    if (op.order == 4 && (op.power == 1 || op.power == 3)) {
      expect(std::abs(chi_on(a, "E1", op)) < 1e-6, "D4d E1(C4) = 0");
      saw_c4 = true;
    }
    // C4^1 / C4^3 have tr = 1; half-turns are order 2 or C4^2 (tr = −1).
    const bool principal_c2 = op.order == 2 || (op.order == 4 && op.power == 2) ||
                              std::abs(op.matrix.trace() + 1.0) < 1e-4;
    if (principal_c2) {
      expect(std::abs(chi_on(a, "E1", op) + 2.0) < 1e-6, "D4d E1(C2) = −2");
      saw_c2 = true;
    }
  }
  expect(saw_c4 && saw_c2, "D4d has principal C4 and C2");
}

void test_s6_inversion() {
  auto a = molsymm::analyze(staggered_ethane());
  expect(a.point_group == "D3d", "staggered ethane D3d got " + a.point_group);
  expect(has_subgroup(a, "S6"), "D3d lists S6");
  a.select_subgroup("S6");
  expect(a.point_group == "S6" && a.order == 6, "ethane→S6");
  expect(has_subgroup(a, "Ci") && !has_subgroup(a, "C3h"), "S6 lists Ci not C3h");
  bool saw_i = false;
  for (const auto& op : a.operations) {
    if (op.kind == molsymm::OperationKind::Inversion) {
      expect(std::abs(chi_on(a, "B", op) + 1.0) < 1e-6, "S6 B(i) = −1");
      expect(std::abs(chi_on(a, "E1", op) + 2.0) < 1e-6, "S6 E1(i) = −2");
      saw_i = true;
    }
  }
  expect(saw_i, "S6 contains inversion");
}

void test_water_b1_mulliken() {
  auto a = molsymm::analyze(water());
  expect(a.point_group == "C2v", "water C2v");
  expect(!has_subgroup(a, "C1v"), "water has no C1v");
  bool saw_yz = false;
  bool saw_xz = false;
  for (const auto& op : a.operations) {
    if (op.kind != molsymm::OperationKind::Reflection || op.axis.norm() < 1e-12) {
      continue;
    }
    const auto n = op.axis.normalized();
    if (std::abs(n.x()) > 0.9) {
      expect(std::abs(chi_on(a, "B1", op) + 1.0) < 1e-6, "water B1(σ_yz) = −1");
      expect(std::abs(chi_on(a, "B2", op) - 1.0) < 1e-6, "water B2(σ_yz) = +1");
      saw_yz = true;
    }
    if (std::abs(n.y()) > 0.9) {
      expect(std::abs(chi_on(a, "B1", op) - 1.0) < 1e-6, "water B1(σ_xz) = +1");
      expect(std::abs(chi_on(a, "B2", op) + 1.0) < 1e-6, "water B2(σ_xz) = −1");
      saw_xz = true;
    }
  }
  expect(saw_yz && saw_xz, "water has σ_yz and σ_xz");
}

void test_cubic_subgroups() {
  {
    auto a = molsymm::analyze(methane());
    expect(a.point_group == "Td", "methane Td");
    for (const auto& [name, order] : {std::pair<const char*, int>{"C2", 2},
                                      {"S4", 4},
                                      {"C2v", 4},
                                      {"D2", 4},
                                      {"D2d", 8},
                                      {"T", 12},
                                      {"Td", 24}}) {
      auto sub = a;
      sub.select_subgroup(name);
      expect(sub.point_group == name && sub.order == order,
             std::string("methane→") + name + " order " + std::to_string(order));
      expect(!sub.character_table.irrep_names.empty(), std::string(name) + " has character table");
    }
  }
  {
    auto a = molsymm::analyze(sf6());
    expect(a.point_group == "Oh", "SF6 Oh");
    for (const auto& [name, order] : {std::pair<const char*, int>{"T", 12},
                                      {"O", 24},
                                      {"Oh", 48},
                                      {"D2d", 8},
                                      {"D4h", 16},
                                      {"D3d", 12}}) {
      auto sub = a;
      sub.select_subgroup(name);
      expect(sub.point_group == name && sub.order == order,
             std::string("SF6→") + name + " order " + std::to_string(order));
      expect(!sub.character_table.irrep_names.empty(), std::string(name) + " has character table");
    }
    expect(has_subgroup(a, "Th"), "SF6 lists Th");
    try {
      auto sub = a;
      sub.select_subgroup("Th");
      expect(false, "SF6→Th should throw");
    } catch (const molsymm::Error&) {
      expect(true, "SF6→Th throws");
    }
  }
}

void test_co2_subgroup_projection() {
  molsymm::Molecule mol = co2();
  Eigen::VectorXd forces = Eigen::VectorXd::Zero(9);
  forces(0) = 0.1;
  forces(1) = 0.2;
  forces(5) = 0.3;
  forces(8) = 0.4;

  {
    auto out = molsymm::symmetrize(mol, "C1", &forces, nullptr);
    expect(std::abs(out.forces(0) - 0.1) < 1e-10, "CO2 C1 keeps Fx");
    expect(std::abs(out.forces(1) - 0.2) < 1e-10, "CO2 C1 keeps Fy");
  }
  {
    auto out = molsymm::symmetrize(mol, "Ci", &forces, nullptr);
    expect(out.forces.segment<3>(0).norm() < 1e-10, "CO2 Ci zeros carbon force");
    expect(std::abs(out.forces(5) + out.forces(8)) < 1e-10, "CO2 Ci antisymmetrizes Fz");
  }
  {
    Eigen::VectorXd f2 = Eigen::VectorXd::Zero(9);
    f2(3) = 0.5;
    f2(6) = -0.5;
    auto out = molsymm::symmetrize(mol, "Ci", &f2, nullptr);
    expect(std::abs(out.forces(3) - 0.5) < 1e-10, "CO2 Ci keeps opposite Fx");
    expect(std::abs(out.forces(6) + 0.5) < 1e-10, "CO2 Ci keeps partner −Fx");
  }
  {
    auto a = molsymm::analyze(mol);
    expect(a.point_group == "Dinfh", "CO2 Dinfh");
    a.select_subgroup("Cinfv");
    expect(a.point_group == "Cinfv" && a.infinite, "CO2→Cinfv");
    auto d = molsymm::unique_displacements(a, 1e-3);
    bool saw_sigma = false;
    for (const auto& ud : d.displacements) {
      if (ud.direction.z() > 0.5) {
        expect(!ud.irreps.empty() && ud.irreps[0] == "Σ+", "Cinfv axial irrep Σ+");
        saw_sigma = true;
      }
    }
    expect(saw_sigma, "Cinfv has axial displacement");
  }
  {
    auto a = molsymm::analyze(mol);
    auto d = molsymm::unique_displacements(a, 1e-3);
    bool saw_u = false;
    bool saw_g = false;
    for (const auto& ud : d.displacements) {
      if (ud.direction.z() < 0.5) {
        continue;
      }
      const auto& p = a.molecule.atoms[static_cast<std::size_t>(ud.atom)].position;
      if (p.norm() < 1e-8) {
        expect(!ud.irreps.empty() && ud.irreps[0] == "Σu+", "CO2 carbon z is Σu+");
        saw_u = true;
      } else {
        expect(!ud.irreps.empty() && ud.irreps[0] == "Σg+", "CO2 oxygen z is Σg+");
        saw_g = true;
      }
    }
    expect(saw_u && saw_g, "CO2 has center Σu+ and end Σg+");
  }
}

} // namespace

void test_group_math() {
  test_benzene_e1g_s_columns();
  test_d4d_characters();
  test_s6_inversion();
  test_water_b1_mulliken();
  test_cubic_subgroups();
  test_co2_subgroup_projection();
  if (g_failed) {
    throw std::runtime_error(std::to_string(g_failed) + " group-math failures (" +
                             std::to_string(g_passed) + " checks passed)");
  }
  std::cout << "group-math checks passed (" << g_passed << ")\n";
}
