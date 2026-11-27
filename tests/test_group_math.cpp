#include "molsymm/molsymm.hpp"

#include "fixtures.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using molsymm::test::benzene;
using molsymm::test::co2;
using molsymm::test::kPi;
using molsymm::test::methane;
using molsymm::test::sf6;
using molsymm::test::water;

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

bool has_class_label(const molsymm::Analysis& a, const std::string& label) {
  const auto& labs = a.character_table.class_labels;
  return std::find(labs.begin(), labs.end(), label) != labs.end();
}

void expect_c2v_b_split(const molsymm::Analysis& a, const std::string& tag) {
  bool saw_c2 = false;
  bool saw_sig = false;
  double first_b1_sig = 0.0;
  for (const auto& op : a.operations) {
    if (op.kind == molsymm::OperationKind::ProperRotation) {
      expect(std::abs(chi_on(a, "B1", op) + 1.0) < 1e-6, tag + " B1(C2) = −1");
      expect(std::abs(chi_on(a, "B2", op) + 1.0) < 1e-6, tag + " B2(C2) = −1");
      saw_c2 = true;
    } else if (op.kind == molsymm::OperationKind::Reflection) {
      expect(std::abs(chi_on(a, "B1", op) + chi_on(a, "B2", op)) < 1e-6, tag + " B1(σ) = −B2(σ)");
      if (!saw_sig) {
        first_b1_sig = chi_on(a, "B1", op);
        saw_sig = true;
      } else {
        expect(std::abs(chi_on(a, "B1", op) - first_b1_sig) > 0.5, tag + " mirrors disagree on B1");
      }
    }
  }
  expect(saw_c2 && saw_sig, tag + " has C2 and σ");
}

molsymm::Molecule staggered_ethane() {
  std::vector<int> z{1, 6, 1, 1, 6, 1, 1, 1};
  std::vector<double> xyz{1.1851,  -0.0039, 0.9875,  0.7516,  -0.0225, -0.0209, 1.1669, 0.8330,
                          -0.5693, 1.1155,  -0.9329, -0.5145, -0.7516, 0.0225,  0.0209, -1.1669,
                          -0.8334, 0.5687,  -1.1157, 0.9326,  0.5151,  -1.1850, 0.0044, -0.9875};
  return molsymm::Molecule(z, xyz);
}

molsymm::Molecule square_antiprism() {
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
  expect(!has_subgroup(a, "C2v"), "D3d does not list C2v");
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
    expect(has_subgroup(a, "Cs"), "Td lists Cs");
    for (const auto& [name, order] : {std::pair<const char*, int>{"C2", 2},
                                      {"Cs", 2},
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
    auto t = a;
    t.select_subgroup("T");
    bool saw_c3 = false;
    bool saw_c2 = false;
    for (const auto& op : t.operations) {
      if (op.kind == molsymm::OperationKind::ProperRotation && op.order == 3) {
        expect(std::abs(chi_on(t, "E", op) + 1.0) < 1e-6, "T E(C3) = −1");
        expect(std::abs(chi_on(t, "A", op) - 1.0) < 1e-6, "T A(C3) = +1");
        saw_c3 = true;
      }
      if (op.kind == molsymm::OperationKind::ProperRotation &&
          (op.order == 2 || std::abs(op.matrix.trace() + 1.0) < 1e-4)) {
        expect(std::abs(chi_on(t, "T", op) + 1.0) < 1e-6, "T T(C2) = −1");
        expect(std::abs(chi_on(t, "E", op) - 2.0) < 1e-6, "T E(C2) = +2");
        saw_c2 = true;
      }
    }
    expect(saw_c3 && saw_c2, "T has C3 and C2 classes");
  }
  {
    auto a = molsymm::analyze(sf6());
    expect(a.point_group == "Oh", "SF6 Oh");
    expect(has_subgroup(a, "Th") && has_subgroup(a, "S4") && has_subgroup(a, "S6"),
           "Oh lists Th/S4/S6");
    for (const auto& [name, order] : {std::pair<const char*, int>{"T", 12},
                                      {"Th", 24},
                                      {"Td", 24},
                                      {"O", 24},
                                      {"Oh", 48},
                                      {"Ci", 2},
                                      {"Cs", 2},
                                      {"C2v", 4},
                                      {"C2h", 4},
                                      {"C4v", 8},
                                      {"C4h", 8},
                                      {"S4", 4},
                                      {"S6", 6},
                                      {"D2", 4},
                                      {"D2h", 8},
                                      {"D2d", 8},
                                      {"D4", 8},
                                      {"D4h", 16},
                                      {"D3", 6},
                                      {"D3d", 12}}) {
      auto sub = a;
      sub.select_subgroup(name);
      expect(sub.point_group == name && sub.order == order,
             std::string("SF6→") + name + " order " + std::to_string(order));
      expect(!sub.character_table.irrep_names.empty(), std::string(name) + " has character table");
    }
    for (const char* cubic : {"T", "Th", "O"}) {
      auto sub = a;
      sub.select_subgroup(cubic);
      sub.select_subgroup("C1");
      expect(sub.point_group == "C1" && sub.order == 1, std::string("SF6→") + cubic + "→C1");
      expect(!sub.character_table.irrep_names.empty(),
             std::string(cubic) + "→C1 has character table");
    }
    auto th = a;
    th.select_subgroup("Th");
    expect(has_subgroup(th, "S6") && has_subgroup(th, "C2v") && !has_subgroup(th, "C3h"),
           "Th lists S6 and C2v not C3h");
    {
      auto c2v = th;
      c2v.select_subgroup("C2v");
      expect(c2v.point_group == "C2v" && c2v.order == 4, "Th→C2v order 4");
      expect_c2v_b_split(c2v, "Th→C2v");
    }
    expect(has_class_label(th, "4C3") && has_class_label(th, "4C3^2"),
           "Th labels split 4C3 / 4C3^2");
    expect(has_class_label(th, "4S6") && has_class_label(th, "4S6^5"),
           "Th labels split 4S6 / 4S6^5");
    bool saw_i = false;
    bool saw_s6 = false;
    bool saw_sigma = false;
    for (const auto& op : th.operations) {
      if (op.kind == molsymm::OperationKind::Inversion) {
        expect(std::abs(chi_on(th, "Eg", op) - 2.0) < 1e-6, "Th Eg(i) = +2");
        expect(std::abs(chi_on(th, "Eu", op) + 2.0) < 1e-6, "Th Eu(i) = −2");
        saw_i = true;
      } else if (op.kind == molsymm::OperationKind::ImproperRotation && op.order == 6) {
        // S6 = i·C3²; S6⁵ = i·C3. Real E: χ_u(S6) = −χ(C3²) = +1.
        expect(std::abs(chi_on(th, "Eu", op) - 1.0) < 1e-6, "Th Eu(S6) = +1");
        expect(std::abs(chi_on(th, "Eg", op) + 1.0) < 1e-6, "Th Eg(S6) = −1");
        saw_s6 = true;
      } else if (op.kind == molsymm::OperationKind::Reflection) {
        expect(std::abs(chi_on(th, "Eu", op) + 2.0) < 1e-6, "Th Eu(σh) = −2");
        expect(std::abs(chi_on(th, "Eg", op) - 2.0) < 1e-6, "Th Eg(σh) = +2");
        saw_sigma = true;
      }
    }
    expect(saw_i && saw_s6 && saw_sigma, "Th has i, S6, and σh");

    {
      auto c2 = a;
      c2.select_subgroup("C2");
      for (const auto& op : c2.operations) {
        if (op.kind == molsymm::OperationKind::ProperRotation) {
          expect(std::abs(chi_on(c2, "B", op) + 1.0) < 1e-6, "SF6→C2 B(C2) = −1");
        }
      }
    }
    {
      auto c2v = a;
      c2v.select_subgroup("C2v");
      expect_c2v_b_split(c2v, "SF6→C2v");
    }
    {
      auto d2 = a;
      d2.select_subgroup("D2");
      double ip_a1_b1 = 0.0;
      for (const auto& op : d2.operations) {
        ip_a1_b1 += chi_on(d2, "A1", op) * chi_on(d2, "B1", op);
      }
      expect(std::abs(ip_a1_b1) < 1e-6, "SF6→D2 ⟨A1|B1⟩ = 0");
    }
    {
      auto d4 = a;
      d4.select_subgroup("D4");
      bool saw_face_c2 = false;
      for (const auto& op : d4.operations) {
        if (op.kind != molsymm::OperationKind::ProperRotation) {
          continue;
        }
        const bool half = (op.order == 2) || (op.order == 4 && op.power == 2) ||
                          std::abs(op.matrix.trace() + 1.0) < 1e-4;
        const bool on_z = op.axis.norm() > 1e-12 && std::abs(op.axis.normalized().z()) > 0.9;
        if (half && !on_z) {
          expect(std::abs(chi_on(d4, "A2", op) + 1.0) < 1e-6, "SF6→D4 A2(C2′) = −1");
          saw_face_c2 = true;
        }
      }
      expect(saw_face_c2, "SF6→D4 has off-axis C2′");
    }
    {
      auto t = a;
      t.select_subgroup("T");
      expect(has_class_label(t, "4C3") && has_class_label(t, "4C3^2"),
             "SF6→T labels split 4C3 / 4C3^2");
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
