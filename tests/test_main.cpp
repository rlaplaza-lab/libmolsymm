#include "molsymm/molsymm.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef MOLSYMM_TEST_DATA_DIR
#error "MOLSYMM_TEST_DATA_DIR must be defined"
#endif

namespace {

constexpr double kPi = 3.14159265358979323846;
int g_failed = 0;

void expect(bool cond, const std::string& msg) {
  if (!cond) {
    std::cerr << "FAIL: " << msg << "\n";
    ++g_failed;
  } else {
    std::cout << "ok: " << msg << "\n";
  }
}

int atomic_number_from_symbol(const std::string& sym) {
  static const char* symbols[] = {
      "",   "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al",
      "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co",
      "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb",
      "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs",
      "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm",
      "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg"};
  for (int z = 1; z <= 80; ++z) {
    if (sym == symbols[z]) {
      return z;
    }
  }
  throw std::runtime_error("unknown element symbol: " + sym);
}

molsymm::Molecule load_xyz(const std::string& path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("cannot open " + path);
  }
  int nat = 0;
  std::string line;
  std::getline(in, line);
  nat = std::stoi(line);
  std::getline(in, line); // comment
  std::vector<int> z;
  std::vector<double> xyz;
  z.reserve(static_cast<std::size_t>(nat));
  xyz.reserve(static_cast<std::size_t>(nat) * 3);
  for (int i = 0; i < nat; ++i) {
    if (!std::getline(in, line)) {
      throw std::runtime_error("truncated XYZ in " + path);
    }
    std::istringstream ss(line);
    std::string sym;
    double x = 0.0;
    double y = 0.0;
    double zz = 0.0;
    ss >> sym >> x >> y >> zz;
    if (!ss) {
      throw std::runtime_error("bad XYZ atom line in " + path);
    }
    z.push_back(atomic_number_from_symbol(sym));
    xyz.push_back(x);
    xyz.push_back(y);
    xyz.push_back(zz);
  }
  return molsymm::Molecule(z, xyz);
}

bool has_subgroup(const molsymm::Analysis& a, const std::string& name) {
  return std::find(a.subgroups.begin(), a.subgroups.end(), name) != a.subgroups.end();
}

molsymm::Molecule water() {
  // O-H2, C2v
  std::vector<int> z{8, 1, 1};
  std::vector<double> xyz{
      0.000000,
      0.000000,
      0.117300,
      -0.000000,
      0.757200,
      -0.469200,
      -0.000000,
      -0.757200,
      -0.469200,
  };
  return molsymm::Molecule(z, xyz);
}

molsymm::Molecule ammonia() {
  std::vector<int> z{7, 1, 1, 1};
  std::vector<double> xyz{
      0.000000,
      0.000000,
      0.000000,
      0.000000,
      0.943863,
      -0.382029,
      0.817408,
      -0.471931,
      -0.382029,
      -0.817408,
      -0.471931,
      -0.382029,
  };
  return molsymm::Molecule(z, xyz);
}

molsymm::Molecule methane() {
  const double a = 0.627;
  std::vector<int> z{6, 1, 1, 1, 1};
  std::vector<double> xyz{
      0,
      0,
      0,
      a,
      a,
      a,
      a,
      -a,
      -a,
      -a,
      a,
      -a,
      -a,
      -a,
      a,
  };
  return molsymm::Molecule(z, xyz);
}

molsymm::Molecule benzene() {
  std::vector<int> z(12);
  std::vector<double> xyz(36);
  for (int i = 0; i < 6; ++i) {
    const double ang = 2.0 * kPi * i / 6.0;
    z[i] = 6;
    xyz[3 * i] = 1.39 * std::cos(ang);
    xyz[3 * i + 1] = 1.39 * std::sin(ang);
    xyz[3 * i + 2] = 0.0;
    z[6 + i] = 1;
    xyz[3 * (6 + i)] = 2.47 * std::cos(ang);
    xyz[3 * (6 + i) + 1] = 2.47 * std::sin(ang);
    xyz[3 * (6 + i) + 2] = 0.0;
  }
  return molsymm::Molecule(z, xyz);
}

molsymm::Molecule allene() {
  // D2d
  std::vector<int> z{6, 6, 6, 1, 1, 1, 1};
  std::vector<double> xyz{
      0,    0,     0, 0,    0, 1.31, 0,     0, -1.31, 0.95,  0,
      1.89, -0.95, 0, 1.89, 0, 0.95, -1.89, 0, -0.95, -1.89,
  };
  return molsymm::Molecule(z, xyz);
}

molsymm::Molecule sf6() {
  std::vector<int> z{16, 9, 9, 9, 9, 9, 9};
  const double r = 1.56;
  std::vector<double> xyz{
      0, 0, 0, r, 0, 0, -r, 0, 0, 0, r, 0, 0, -r, 0, 0, 0, r, 0, 0, -r,
  };
  return molsymm::Molecule(z, xyz);
}

molsymm::Molecule co2() {
  std::vector<int> z{6, 8, 8};
  std::vector<double> xyz{0, 0, 0, 0, 0, 1.16, 0, 0, -1.16};
  return molsymm::Molecule(z, xyz);
}

molsymm::Molecule co() {
  std::vector<int> z{6, 8};
  std::vector<double> xyz{0, 0, 0, 0, 0, 1.13};
  return molsymm::Molecule(z, xyz);
}

void test_point_groups() {
  auto check = [](const molsymm::Molecule& mol, const std::string& expect_pg, int sigma) {
    auto a = molsymm::analyze(mol);
    expect(a.point_group == expect_pg, "group " + expect_pg + " got " + a.point_group);
    expect(a.symmetry_number == sigma,
           "sigma " + expect_pg + " expect " + std::to_string(sigma) + " got " +
               std::to_string(a.symmetry_number));
  };

  check(water(), "C2v", 2);
  check(ammonia(), "C3v", 3);
  check(methane(), "Td", 12);
  check(benzene(), "D6h", 12);
  check(allene(), "D2d", 4);
  check(sf6(), "Oh", 24);
  check(co2(), "Dinfh", 2);
  check(co(), "Cinfv", 1);

  {
    molsymm::Molecule m;
    m.atoms.push_back({8, {0, 0, 0}, 0});
    m.resolve_masses();
    auto a = molsymm::analyze(m);
    expect(a.point_group == "Kh", "single atom Kh");
  }

  {
    auto w = water();
    w.atoms[1].position.x() += 0.05;
    molsymm::Tolerances loose;
    loose.primary = 1e-2;
    auto a = molsymm::analyze(w, loose);
    expect(a.point_group == "C2v", "perturbed water recovers C2v got " + a.point_group);
    a.symmetrize();
  }
}

void expect_character_table_consistent(const molsymm::Analysis& a, const std::string& label) {
  const auto& ct = a.character_table;
  expect(!ct.irrep_names.empty(), label + " character table non-empty");
  if (a.infinite || a.order <= 0 || ct.irrep_dims.empty() || ct.characters.empty()) {
    return;
  }
  int dim2 = 0;
  for (int d : ct.irrep_dims) {
    dim2 += d * d;
  }
  expect(dim2 == a.order,
         label + " sum dim^2 expect " + std::to_string(a.order) + " got " + std::to_string(dim2));
  // Totally symmetric irrep is first and +1 on every class.
  expect(!ct.characters.front().empty(), label + " totally symmetric row present");
  bool all_one = true;
  for (double chi : ct.characters.front()) {
    if (std::abs(chi - 1.0) > 1e-8) {
      all_one = false;
      break;
    }
  }
  expect(all_one, label + " totally symmetric row is +1 on every class");
}

void test_subgroups_and_table() {
  auto a = molsymm::analyze(benzene());
  expect(has_subgroup(a, "D6") && has_subgroup(a, "C6v") && has_subgroup(a, "C2v"),
         "D6h subgroups include D6, C6v, C2v");
  expect_character_table_consistent(a, "benzene");

  auto d6 = a;
  d6.select_subgroup("D6");
  expect(d6.point_group == "D6" && d6.order == 12, "benzene select_subgroup D6");
  expect_character_table_consistent(d6, "benzene/D6");

  auto c6v = a;
  c6v.select_subgroup("C6v");
  expect(c6v.point_group == "C6v" && c6v.order == 12, "benzene select_subgroup C6v");
  expect_character_table_consistent(c6v, "benzene/C6v");

  auto m = molsymm::analyze(methane());
  m.select_subgroup("T");
  expect(m.point_group == "T" && m.order == 12, "methane select_subgroup T");
  expect_character_table_consistent(m, "methane/T");
}

void test_displacements_water() {
  auto a = molsymm::analyze(water());
  expect_character_table_consistent(a, "water");
  auto d = molsymm::unique_displacements(a, 1e-3);
  expect(static_cast<int>(d.displacements.size()) == 6,
         "water has 6 unique displacements, got " + std::to_string(d.displacements.size()));
  for (std::size_t u = 0; u < d.displacements.size(); ++u) {
    expect(!d.displacements[u].irreps.empty(),
           "water displacement " + std::to_string(u) + " has irreps");
  }

  // Planted Hessian: identity on 3N (in original frame)
  const int n = a.original.size();
  const int n3 = 3 * n;
  Eigen::MatrixXd H = Eigen::MatrixXd::Identity(n3, n3);
  // Make it totally symmetric by projecting
  H = molsymm::symmetrize_hessian(a, H);

  auto energy = [&](const molsymm::Molecule& m) {
    Eigen::VectorXd x(n3);
    for (int i = 0; i < n; ++i) {
      molsymm::Vec3 dr = m.atoms[static_cast<std::size_t>(i)].position -
                         a.original.atoms[static_cast<std::size_t>(i)].position;
      x.segment<3>(3 * i) = dr;
    }
    return 0.5 * x.dot(H * x);
  };
  auto grad = [&](const molsymm::Molecule& m) {
    Eigen::VectorXd x(n3);
    for (int i = 0; i < n; ++i) {
      molsymm::Vec3 dr = m.atoms[static_cast<std::size_t>(i)].position -
                         a.original.atoms[static_cast<std::size_t>(i)].position;
      x.segment<3>(3 * i) = dr;
    }
    return Eigen::VectorXd(H * x);
  };

  std::vector<double> ep(d.displacements.size());
  std::vector<double> em(d.displacements.size());
  for (std::size_t u = 0; u < d.displacements.size(); ++u) {
    ep[u] = energy(d.plus[u]);
    em[u] = energy(d.minus[u]);
  }
  Eigen::VectorXd g_rec = d.gradient_from_energies(ep, em);
  expect(g_rec.norm() < 1e-6, "water equilibrium reconstructed gradient ~0");

  std::vector<Eigen::VectorXd> gp(d.displacements.size());
  std::vector<Eigen::VectorXd> gm(d.displacements.size());
  for (std::size_t u = 0; u < d.displacements.size(); ++u) {
    gp[u] = grad(d.plus[u]);
    gm[u] = grad(d.minus[u]);
  }
  Eigen::MatrixXd H_rec = d.hessian_from_gradients(gp, gm);
  H_rec = molsymm::symmetrize_hessian(a, H_rec);
  const double err = (H_rec - H).norm() / (H.norm() + 1e-15);
  expect(err < 1e-6, "water hessian rebuild relative err " + std::to_string(err));
}

void test_displacements_methane() {
  auto a = molsymm::analyze(methane());
  auto d = molsymm::unique_displacements(a, 1e-3);
  expect(static_cast<int>(d.displacements.size()) == 3,
         "methane has 3 unique displacements, got " + std::to_string(d.displacements.size()));
}

struct XyzCase {
  const char* file;
  const char* expect_pg;
  int expect_order;
  int expect_sigma;
  int expect_ndisp;
  double primary;
  double max_rmsd;
  std::vector<const char*> subgroups;
};

const std::vector<XyzCase>& xyz_corpus() {
  // 22paracyclophane is D2h chemically, but the published coordinates only recover
  // D2 (missing σh) even at primary=5e-3; assert what analyze returns.
  // cyclopropane needs primary=1e-2 to recover D3h (ideal chemistry) from the
  // IQmol fragment coordinates.
  static const std::vector<XyzCase> cases = {
      {"water.xyz", "C2v", 4, 2, 6, 1e-3, 1e-3, {"Cs", "C2"}},
      {"anion-pi.xyz", "C6v", 12, 6, 8, 1e-3, 1e-3, {"C6", "C3v", "C2v", "Cs"}},
      {"22paracyclophane.xyz", "D2", 4, 4, 24, 5e-3, 1e-3, {"D2", "C2"}},
      {"E_2butene.xyz", "C2h", 4, 2, 15, 1e-3, 1e-3, {"C2", "Ci"}},
      {"pt_complex.xyz", "C2v", 4, 2, 15, 1e-3, 1e-3, {"Cs", "C2"}},
      {"c60.xyz", "Ih", 120, 60, 3, 1e-3, 1e-3, {"I", "C5", "C5v", "D5d", "Th"}},
      {"cubane.xyz", "Oh", 48, 24, 2, 1e-3, 1e-3, {"O", "Td", "C4v", "D4h"}},
      {"adamantane.xyz", "Td", 24, 12, 9, 1e-3, 1e-3, {"T", "C3v", "C2v"}},
      {"coronene.xyz", "D6h", 24, 12, 10, 1e-3, 1e-3, {"D6", "C6v", "C2v"}},
      {"corannulene.xyz", "C5v", 10, 5, 10, 1e-3, 1e-3, {"C5", "Cs"}},
      {"ferrocene.xyz", "D5d", 20, 10, 5, 1e-3, 1e-3, {"D5", "C5v", "S10"}},
      {"cyclopropane.xyz", "D3h", 12, 6, 3, 1e-2, 5e-3, {"D3", "C3v", "C3h"}},
      {"octahedral.xyz", "Oh", 48, 24, 2, 1e-3, 1e-3, {"O", "C4v", "D4h"}},
  };
  return cases;
}

void check_xyz_case(const XyzCase& c) {
  const std::string path = std::string(MOLSYMM_TEST_DATA_DIR) + "/" + c.file;
  molsymm::Molecule mol = load_xyz(path);
  molsymm::Tolerances tol;
  tol.primary = c.primary;
  auto a = molsymm::analyze(mol, tol);
  expect(a.point_group == c.expect_pg,
         std::string(c.file) + " group " + c.expect_pg + " got " + a.point_group);
  expect(a.order == c.expect_order,
         std::string(c.file) + " order expect " + std::to_string(c.expect_order) + " got " +
             std::to_string(a.order));
  expect(a.symmetry_number == c.expect_sigma,
         std::string(c.file) + " sigma expect " + std::to_string(c.expect_sigma) + " got " +
             std::to_string(a.symmetry_number));
  expect_character_table_consistent(a, c.file);
  expect(a.rmsd < c.max_rmsd,
         std::string(c.file) + " rmsd " + std::to_string(a.rmsd) + " max " +
             std::to_string(c.max_rmsd));
  for (const char* sg : c.subgroups) {
    expect(has_subgroup(a, sg), std::string(c.file) + " has subgroup " + sg);
  }
  auto d = molsymm::unique_displacements(a, 1e-3);
  const int n3 = 3 * a.original.size();
  const int nd = static_cast<int>(d.displacements.size());
  expect(nd >= 1 && nd < n3,
         std::string(c.file) + " unique disp in [1,3N) got " + std::to_string(nd) + " of " +
             std::to_string(n3));
  expect(nd == c.expect_ndisp,
         std::string(c.file) + " unique disp expect " + std::to_string(c.expect_ndisp) + " got " +
             std::to_string(nd));
}

void test_xyz_corpus() {
  for (const XyzCase& c : xyz_corpus()) {
    try {
      check_xyz_case(c);
    } catch (const std::exception& ex) {
      expect(false, std::string(c.file) + " exception: " + ex.what());
    }
  }
}

molsymm::Molecule rattle_molecule(molsymm::Molecule mol, double amplitude, unsigned seed) {
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> dist(-amplitude, amplitude);
  for (auto& atom : mol.atoms) {
    atom.position.x() += dist(rng);
    atom.position.y() += dist(rng);
    atom.position.z() += dist(rng);
  }
  return mol;
}

bool is_iqmol_high_symmetry(const char* file) {
  static const char* const files[] = {
      "c60.xyz",
      "cubane.xyz",
      "adamantane.xyz",
      "coronene.xyz",
      "corannulene.xyz",
      "ferrocene.xyz",
      "cyclopropane.xyz",
      "octahedral.xyz",
  };
  return std::any_of(
      std::begin(files), std::end(files), [&](const char* f) { return std::string(f) == file; });
}

void test_xyz_slight_rattle() {
  // IQmol high-symmetry geometries: tiny noise must not change the group.
  unsigned seed = 42;
  for (const XyzCase& c : xyz_corpus()) {
    if (!is_iqmol_high_symmetry(c.file)) {
      continue;
    }
    try {
      const std::string path = std::string(MOLSYMM_TEST_DATA_DIR) + "/" + c.file;
      molsymm::Molecule mol = rattle_molecule(load_xyz(path), 1e-4, seed++);
      molsymm::Tolerances tol;
      tol.primary = c.primary;
      auto a = molsymm::analyze(mol, tol);
      expect(a.point_group == c.expect_pg,
             std::string(c.file) + " slight rattle group " + c.expect_pg + " got " + a.point_group);
      expect(a.order == c.expect_order,
             std::string(c.file) + " slight rattle order " + std::to_string(c.expect_order) +
                 " got " + std::to_string(a.order));
      expect(a.symmetry_number == c.expect_sigma,
             std::string(c.file) + " slight rattle sigma " + std::to_string(c.expect_sigma) +
                 " got " + std::to_string(a.symmetry_number));
      expect(a.rmsd < c.max_rmsd,
             std::string(c.file) + " slight rattle rmsd " + std::to_string(a.rmsd));
    } catch (const std::exception& ex) {
      expect(false, std::string(c.file) + " slight rattle exception: " + ex.what());
    }
  }
}

void test_xyz_large_rattle() {
  // Larger noise: default tolerance drops the ideal group; looser tolerance recovers it.
  struct LargeRattleCase {
    const char* file;
    const char* expect_pg;
  };
  const LargeRattleCase cases[] = {
      {"cubane.xyz", "Oh"},
      {"cyclopropane.xyz", "D3h"},
      {"octahedral.xyz", "Oh"},
  };
  unsigned seed = 99;
  for (const LargeRattleCase& c : cases) {
    try {
      const std::string path = std::string(MOLSYMM_TEST_DATA_DIR) + "/" + c.file;
      molsymm::Molecule mol = rattle_molecule(load_xyz(path), 5e-3, seed++);
      molsymm::Tolerances tight;
      tight.primary = 1e-3;
      auto a_tight = molsymm::analyze(mol, tight);
      expect(a_tight.point_group != c.expect_pg,
             std::string(c.file) + " large rattle at 1e-3 should drop " + c.expect_pg + " got " +
                 a_tight.point_group);
      molsymm::Tolerances loose;
      loose.primary = 2e-2;
      auto a_loose = molsymm::analyze(mol, loose);
      expect(a_loose.point_group == c.expect_pg,
             std::string(c.file) + " large rattle at 2e-2 recovers " + c.expect_pg + " got " +
                 a_loose.point_group);
    } catch (const std::exception& ex) {
      expect(false, std::string(c.file) + " large rattle exception: " + ex.what());
    }
  }
}

} // namespace

int main() {
  test_point_groups();
  test_subgroups_and_table();
  test_displacements_water();
  test_displacements_methane();
  test_xyz_corpus();
  test_xyz_slight_rattle();
  test_xyz_large_rattle();
  if (g_failed) {
    std::cerr << g_failed << " failures\n";
    return 1;
  }
  std::cout << "all tests passed\n";
  return 0;
}
