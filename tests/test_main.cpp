#include "molsymm/molsymm.hpp"

#include <Eigen/Geometry>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#ifndef MOLSYMM_TEST_DATA_DIR
#error "MOLSYMM_TEST_DATA_DIR must be defined"
#endif

namespace {

constexpr double kPi = 3.14159265358979323846;
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
    auto s4 = molsymm::analyze(allene());
    s4.select_subgroup("S4");
    expect(s4.point_group == "S4" && s4.order == 4, "allene select_subgroup S4");
    expect(s4.symmetry_number == 2,
           "S4 symmetry number expect 2 got " + std::to_string(s4.symmetry_number));
  }

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

bool is_pure_cn_or_sn(const std::string& pg) {
  if (pg.size() < 2 || (pg[0] != 'C' && pg[0] != 'S')) {
    return false;
  }
  for (std::size_t i = 1; i < pg.size(); ++i) {
    if (!std::isdigit(static_cast<unsigned char>(pg[i]))) {
      return false;
    }
  }
  return true;
}

void expect_character_table_consistent(const molsymm::Analysis& a, const std::string& label) {
  const auto& ct = a.character_table;
  expect(!ct.irrep_names.empty(), label + " character table non-empty");
  if (a.infinite || a.order <= 0 || ct.characters.empty()) {
    return;
  }
  bool all_one = true;
  for (double chi : ct.characters.front()) {
    if (std::abs(chi - 1.0) > 1e-8) {
      all_one = false;
      break;
    }
  }
  expect(all_one, label + " totally symmetric row is +1 on every class");
  if (is_pure_cn_or_sn(a.point_group) || ct.irrep_dims.empty()) {
    return;
  }
  int dim2 = 0;
  for (int d : ct.irrep_dims) {
    dim2 += d * d;
  }
  expect(dim2 == a.order,
         label + " sum dim^2 expect " + std::to_string(a.order) + " got " + std::to_string(dim2));
}

bool operations_form_finite_group(const molsymm::Analysis& a) {
  if (a.infinite || a.order <= 0) {
    return true;
  }
  const auto& ops = a.operations;
  if (static_cast<int>(ops.size()) != a.order) {
    return false;
  }
  constexpr double kOrthTol = 1e-5;
  constexpr double kCloseTol = 5e-3;
  for (std::size_t i = 0; i < ops.size(); ++i) {
    const molsymm::Mat3& R = ops[i].matrix;
    if ((R.transpose() * R - molsymm::Mat3::Identity()).cwiseAbs().maxCoeff() > kOrthTol) {
      return false;
    }
    if (std::abs(std::abs(R.determinant()) - 1.0) > kOrthTol) {
      return false;
    }
    for (std::size_t j = i + 1; j < ops.size(); ++j) {
      if ((ops[i].matrix - ops[j].matrix).cwiseAbs().maxCoeff() <= 1e-6) {
        return false;
      }
    }
  }
  for (const auto& a_op : ops) {
    for (const auto& b_op : ops) {
      const molsymm::Mat3 P = a_op.matrix * b_op.matrix;
      bool found = false;
      for (const auto& c_op : ops) {
        if ((c_op.matrix - P).cwiseAbs().maxCoeff() <= kCloseTol) {
          found = true;
          break;
        }
      }
      if (!found) {
        return false;
      }
    }
  }
  return true;
}

void test_subgroups_and_table() {
  auto a = molsymm::analyze(benzene());
  expect(has_subgroup(a, "D6") && has_subgroup(a, "C6v") && has_subgroup(a, "C2v"),
         "D6h subgroups include D6, C6v, C2v");
  expect(operations_form_finite_group(a), "benzene closed group");
  expect_character_table_consistent(a, "benzene");

  for (const auto& [name, order] : {std::pair{"D6", 12}, {"C6v", 12}, {"C2", 2}, {"C3", 3}}) {
    auto sub = a;
    sub.select_subgroup(name);
    expect(sub.point_group == name && sub.order == order,
           std::string("benzene select_subgroup ") + name);
    expect_character_table_consistent(sub, std::string("benzene/") + name);
  }

  {
    molsymm::Tolerances tol;
    tol.primary = 1e-2;
    auto d3h = molsymm::analyze(load_xyz(std::string(MOLSYMM_TEST_DATA_DIR) +
                                         "/iqmol/Molecules/Cyclic_Alkanes/Cyclopropane.xyz"),
                                tol);
    expect(d3h.point_group == "D3h", "cyclopropane is D3h got " + d3h.point_group);
    expect(has_subgroup(d3h, "C2v") && !has_subgroup(d3h, "Ci") && !has_subgroup(d3h, "D2h"),
           "D3h subgroups");
    auto dn = d3h;
    dn.select_subgroup("D3");
    expect(dn.point_group == "D3" && dn.order == 6 && has_subgroup(dn, "C2") &&
               !has_subgroup(dn, "C2v"),
           "D3 from D3h");
  }

  auto m = molsymm::analyze(methane());
  m.select_subgroup("T");
  expect(m.point_group == "T" && m.order == 12, "methane select_subgroup T");
  expect_character_table_consistent(m, "methane/T");

  {
    auto w = molsymm::analyze(water());
    expect(operations_form_finite_group(w), "water closed group");
    const auto& ct = w.character_table;
    int i_b1 = -1;
    int i_b2 = -1;
    for (std::size_t i = 0; i < ct.irrep_names.size(); ++i) {
      if (ct.irrep_names[i] == "B1") {
        i_b1 = static_cast<int>(i);
      }
      if (ct.irrep_names[i] == "B2") {
        i_b2 = static_cast<int>(i);
      }
    }
    expect(i_b1 >= 0 && i_b2 >= 0, "water has B1 and B2");
    bool split = false;
    for (std::size_t c = 0; c < ct.characters[static_cast<std::size_t>(i_b1)].size(); ++c) {
      if (std::abs(ct.characters[static_cast<std::size_t>(i_b1)][c] -
                   ct.characters[static_cast<std::size_t>(i_b2)][c]) > 1e-8) {
        split = true;
        break;
      }
    }
    expect(split, "water B1 and B2 differ");
  }

  {
    auto d2d = molsymm::analyze(allene());
    expect(operations_form_finite_group(d2d), "allene closed group");
    const auto& names = d2d.character_table.irrep_names;
    expect(std::find(names.begin(), names.end(), "E") != names.end(), "D2d has E");
    expect(std::none_of(names.begin(),
                        names.end(),
                        [](const std::string& n) {
                          return !n.empty() && (n.back() == 'g' || n.back() == 'u');
                        }),
           "D2d irreps have no g/u");
  }
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

struct ManifestRow {
  std::string rel_path;
  std::string expect_pg;
  int expect_order = 1;
  int expect_sigma = 1;
  double primary = 1e-3;
  double max_rmsd = 1e-3;
};

std::vector<ManifestRow> load_manifest(const std::string& path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("cannot open manifest " + path);
  }
  std::vector<ManifestRow> rows;
  std::string line;
  std::getline(in, line); // header
  while (std::getline(in, line)) {
    if (line.empty()) {
      continue;
    }
    std::istringstream ss(line);
    ManifestRow r;
    std::string order_s;
    std::string sigma_s;
    std::string primary_s;
    std::string rmsd_s;
    std::string source;
    if (!std::getline(ss, r.rel_path, '\t') || !std::getline(ss, r.expect_pg, '\t') ||
        !std::getline(ss, order_s, '\t') || !std::getline(ss, sigma_s, '\t') ||
        !std::getline(ss, primary_s, '\t') || !std::getline(ss, rmsd_s, '\t') ||
        !std::getline(ss, source, '\t')) {
      throw std::runtime_error("bad manifest line: " + line);
    }
    r.expect_order = std::stoi(order_s);
    r.expect_sigma = std::stoi(sigma_s);
    r.primary = std::stod(primary_s);
    r.max_rmsd = std::stod(rmsd_s);
    rows.push_back(std::move(r));
  }
  return rows;
}

bool operations_permute_atoms(const molsymm::Analysis& a, double tol) {
  if (a.operations.empty()) {
    return false;
  }
  if (!a.infinite && a.order > 0 && static_cast<int>(a.operations.size()) != a.order) {
    return false;
  }
  const int n = a.molecule.size();
  const double match = std::max(tol, 1e-4);
  for (const molsymm::Operation& op : a.operations) {
    std::vector<bool> claimed(static_cast<std::size_t>(n), false);
    for (int i = 0; i < n; ++i) {
      const molsymm::Vec3 rp = op.matrix * a.molecule.atoms[static_cast<std::size_t>(i)].position;
      const int zi = a.molecule.atoms[static_cast<std::size_t>(i)].atomic_number;
      int best = -1;
      double best_d = match;
      for (int j = 0; j < n; ++j) {
        if (claimed[static_cast<std::size_t>(j)]) {
          continue;
        }
        if (a.molecule.atoms[static_cast<std::size_t>(j)].atomic_number != zi) {
          continue;
        }
        const double d = (rp - a.molecule.atoms[static_cast<std::size_t>(j)].position).norm();
        if (d <= best_d) {
          best_d = d;
          best = j;
        }
      }
      if (best < 0) {
        return false;
      }
      claimed[static_cast<std::size_t>(best)] = true;
    }
  }
  return true;
}

void check_manifest_row(const ManifestRow& c, const std::string& data_root) {
  molsymm::Tolerances tol;
  tol.primary = c.primary;
  auto a = molsymm::analyze(load_xyz(data_root + "/" + c.rel_path), tol);
  expect(a.point_group == c.expect_pg,
         c.rel_path + " group " + c.expect_pg + " got " + a.point_group);
  expect(a.order == c.expect_order,
         c.rel_path + " order expect " + std::to_string(c.expect_order) + " got " +
             std::to_string(a.order));
  expect(a.symmetry_number == c.expect_sigma,
         c.rel_path + " sigma expect " + std::to_string(c.expect_sigma) + " got " +
             std::to_string(a.symmetry_number));
  expect(a.rmsd < c.max_rmsd,
         c.rel_path + " rmsd " + std::to_string(a.rmsd) + " max " + std::to_string(c.max_rmsd));
  expect(operations_permute_atoms(a, c.primary), c.rel_path + " ops permute atoms");
  if (!a.infinite) {
    expect_character_table_consistent(a, c.rel_path);
    auto d = molsymm::unique_displacements(a, 1e-3);
    const int nd = static_cast<int>(d.displacements.size());
    expect(nd >= 1 && nd <= 3 * a.original.size(), c.rel_path + " unique disp in [1,3N]");
  }
}

void test_high_symmetry_counts() {
  const std::string root = MOLSYMM_TEST_DATA_DIR;
  auto check_nd = [&](const char* rel, double primary, int want) {
    molsymm::Tolerances tol;
    tol.primary = primary;
    auto a = molsymm::analyze(load_xyz(root + "/" + rel), tol);
    const int nd = static_cast<int>(molsymm::unique_displacements(a, 1e-3).displacements.size());
    expect(nd == want,
           std::string(rel) + " unique disp expect " + std::to_string(want) + " got " +
               std::to_string(nd));
  };
  check_nd("iqmol/Molecules/Cyclic_Alkanes/Cubane.xyz", 1e-3, 2);
  check_nd("iqmol/Molecules/Cyclic_Alkanes/Cyclopropane.xyz", 1e-2, 3);
  check_nd("iqmol/Molecules/Fullerenes/Buckminsterfullerene.xyz", 1e-4, 8);

  molsymm::Tolerances tol;
  tol.primary = 1e-4;
  auto c60 = molsymm::analyze(
      load_xyz(root + "/iqmol/Molecules/Fullerenes/Buckminsterfullerene.xyz"), tol);
  int n_id = 0;
  int n_c5 = 0;
  int n_c52 = 0;
  int n_c3 = 0;
  int n_c2 = 0;
  for (const auto& op : c60.operations) {
    if (op.kind == molsymm::OperationKind::Identity) {
      ++n_id;
    } else if (op.kind == molsymm::OperationKind::ProperRotation && op.order == 5) {
      const int p = op.power > 0 ? op.power % 5 : 1;
      ((p == 1 || p == 4) ? n_c5 : n_c52)++;
    } else if (op.kind == molsymm::OperationKind::ProperRotation && op.order == 3) {
      ++n_c3;
    } else if (op.kind == molsymm::OperationKind::ProperRotation && op.order == 2) {
      ++n_c2;
    }
  }
  expect(n_id == 1 && n_c5 == 12 && n_c52 == 12 && n_c3 == 20 && n_c2 == 15,
         "C60 Ih proper class sizes");
}

void test_xyz_corpus() {
  const std::string data_root = MOLSYMM_TEST_DATA_DIR;
  const auto rows = load_manifest(data_root + "/manifest.tsv");
  expect(!rows.empty(), "manifest non-empty");
  for (const ManifestRow& c : rows) {
    try {
      check_manifest_row(c, data_root);
    } catch (const std::exception& ex) {
      expect(false, c.rel_path + " exception: " + ex.what());
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

molsymm::Molecule rotate_molecule(molsymm::Molecule mol, const molsymm::Mat3& R) {
  for (auto& atom : mol.atoms) {
    atom.position = R * atom.position;
  }
  return mol;
}

void test_xyz_slight_rattle() {
  const std::string data_root = MOLSYMM_TEST_DATA_DIR;
  const char* files[] = {
      "iqmol/Molecules/Fullerenes/Buckminsterfullerene.xyz",
      "iqmol/Molecules/Cyclic_Alkanes/Cubane.xyz",
      "iqmol/Molecules/Cyclic_Alkanes/Adamantane.xyz",
      "iqmol/Molecules/Polyaromatic_Hydrocarbons/Coronene.xyz",
      "iqmol/Molecules/Polyaromatic_Hydrocarbons/Corannulene.xyz",
      "iqmol/Molecules/Organometalics/Ferrocene.xyz",
      "iqmol/Molecules/Cyclic_Alkanes/Cyclopropane.xyz",
      "iqmol/Molecules/Coordination/6-Octahedral.xyz",
  };
  const char* expect_pg[] = {"Ih", "Oh", "Td", "D6h", "C5v", "D5d", "D3h", "Oh"};
  const double primary[] = {1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-2, 1e-3};
  unsigned seed = 42;
  for (std::size_t i = 0; i < std::size(files); ++i) {
    try {
      molsymm::Molecule mol = rattle_molecule(load_xyz(data_root + "/" + files[i]), 1e-4, seed++);
      molsymm::Tolerances tol;
      tol.primary = primary[i];
      auto a = molsymm::analyze(mol, tol);
      expect(a.point_group == expect_pg[i],
             std::string(files[i]) + " slight rattle group " + expect_pg[i] + " got " +
                 a.point_group);
    } catch (const std::exception& ex) {
      expect(false, std::string(files[i]) + " slight rattle exception: " + ex.what());
    }
  }
}

void test_xyz_large_rattle() {
  const std::string data_root = MOLSYMM_TEST_DATA_DIR;
  struct LargeRattleCase {
    const char* file;
    const char* expect_pg;
  };
  const LargeRattleCase cases[] = {
      {"iqmol/Molecules/Cyclic_Alkanes/Cubane.xyz", "Oh"},
      {"iqmol/Molecules/Cyclic_Alkanes/Cyclopropane.xyz", "D3h"},
      {"iqmol/Molecules/Coordination/6-Octahedral.xyz", "Oh"},
  };
  unsigned seed = 99;
  for (const LargeRattleCase& c : cases) {
    try {
      molsymm::Molecule mol = rattle_molecule(load_xyz(data_root + "/" + c.file), 5e-3, seed++);
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

bool frame_coords_match(const molsymm::Analysis& a, const molsymm::Analysis& b, double tol) {
  if (a.molecule.size() != b.molecule.size()) {
    return false;
  }
  for (int i = 0; i < a.molecule.size(); ++i) {
    if ((a.molecule.atoms[static_cast<std::size_t>(i)].position -
         b.molecule.atoms[static_cast<std::size_t>(i)].position)
            .norm() > tol) {
      return false;
    }
  }
  if (a.operations.size() != b.operations.size()) {
    return false;
  }
  for (std::size_t i = 0; i < a.operations.size(); ++i) {
    if ((a.operations[i].axis - b.operations[i].axis).norm() > 1e-6 &&
        (a.operations[i].axis + b.operations[i].axis).norm() > 1e-6) {
      // Axis direction may flip for improper ops; require parallel.
      const double c =
          std::abs(a.operations[i].axis.normalized().dot(b.operations[i].axis.normalized()));
      if (c < 0.999) {
        return false;
      }
    }
  }
  return true;
}

void test_inertia_frame_stable() {
  const molsymm::Mat3 R =
      Eigen::AngleAxisd(0.7, molsymm::Vec3(1, 2, 3).normalized()).toRotationMatrix();

  auto check = [&](const molsymm::Molecule& mol, const std::string& label, bool check_frame) {
    auto a0 = molsymm::analyze(mol);
    auto a1 = molsymm::analyze(rotate_molecule(mol, R));
    expect(a0.point_group == a1.point_group, label + " rotated group " + a1.point_group);
    if (check_frame) {
      expect(frame_coords_match(a0, a1, 1e-8), label + " rotated analysis frame matches");
    }
  };

  check(water(), "water asymmetric/C2v", true);
  check(co2(), "co2 linear/Dinfh", true);

  auto axis_on_z = [](const molsymm::Analysis& a, int expect_n) {
    int best_n = 1;
    for (const auto& op : a.operations) {
      if (op.kind == molsymm::OperationKind::ProperRotation) {
        best_n = std::max(best_n, op.order);
      }
    }
    if (best_n != expect_n) {
      return false;
    }
    return std::any_of(a.operations.begin(), a.operations.end(), [&](const molsymm::Operation& op) {
      return op.kind == molsymm::OperationKind::ProperRotation && op.order == best_n &&
             std::abs(op.axis.normalized().z()) > 0.9;
    });
  };
  for (const auto& [mol, label, n] :
       {std::tuple{ammonia(), "ammonia", 3}, std::tuple{methane(), "methane", 3}}) {
    auto a0 = molsymm::analyze(mol);
    auto a1 = molsymm::analyze(rotate_molecule(mol, R));
    expect(a0.point_group == a1.point_group && a0.orientation.determinant() > 0.0 &&
               a1.orientation.determinant() > 0.0 && axis_on_z(a0, n) && axis_on_z(a1, n),
           std::string(label) + " oriented with Cn on z");
  }
}

void expect_near(const molsymm::Vec3& got, const molsymm::Vec3& want, const std::string& msg) {
  expect((got - want).norm() < 1e-12, msg + " |dr|=" + std::to_string((got - want).norm()));
}

template <typename Fn>
bool throws_error(Fn&& fn) {
  try {
    fn();
  } catch (const molsymm::Error&) {
    return true;
  } catch (const std::exception& ex) {
    expect(false, std::string("unexpected exception: ") + ex.what());
    return false;
  }
  return false;
}

void test_zmatrix() {
  using molsymm::ZMatrixAtom;
  const double half_pi = 0.5 * kPi;

  const std::vector<ZMatrixAtom> chain = {
      {6, 0, 0, 0, 0.0, 0.0, 0.0},
      {1, 1, 0, 0, 1.0, 0.0, 0.0},
      {1, 1, 2, 0, 1.0, half_pi, 0.0},
      {1, 1, 2, 3, 2.0, half_pi, half_pi},
  };
  const std::vector<molsymm::Vec3> pos = molsymm::zmatrix_positions(chain);
  expect_near(pos[0], molsymm::Vec3(0.0, 0.0, 0.0), "z-matrix atom 1 at origin");
  expect_near(pos[1], molsymm::Vec3(0.0, 1.0, 0.0), "z-matrix atom 2 along +y");
  expect_near(pos[2], molsymm::Vec3(1.0, 0.0, 0.0), "z-matrix atom 3 in the xy plane");
  expect_near(pos[3], molsymm::Vec3(0.0, 0.0, -2.0), "z-matrix atom 4 dihedral");

  const double water_angle = 104.5 * kPi / 180.0;
  const std::vector<ZMatrixAtom> water = {
      {8, 0, 0, 0, 0.0, 0.0, 0.0},
      {1, 1, 0, 0, 0.957, 0.0, 0.0},
      {1, 1, 2, 0, 0.957, water_angle, 0.0},
  };
  const molsymm::Molecule water_mol = molsymm::molecule_from_zmatrix(water);
  expect(water_mol.atoms[0].atomic_number == 8 && water_mol.atoms[2].atomic_number == 1,
         "z-matrix atomic numbers are kept");
  const molsymm::Analysis water_pg = molsymm::analyze(water_mol);
  expect(water_pg.point_group == "C2v", "water z-matrix is C2v, got " + water_pg.point_group);

  const double tet = std::acos(-1.0 / 3.0);
  const double twist = 2.0 * kPi / 3.0;
  const std::vector<ZMatrixAtom> methane = {
      {6, 0, 0, 0, 0.0, 0.0, 0.0},
      {1, 1, 0, 0, 1.09, 0.0, 0.0},
      {1, 1, 2, 0, 1.09, tet, 0.0},
      {1, 1, 2, 3, 1.09, tet, twist},
      {1, 1, 2, 3, 1.09, tet, -twist},
  };
  const molsymm::Analysis methane_pg = molsymm::analyze(molsymm::molecule_from_zmatrix(methane));
  expect(methane_pg.point_group == "Td", "methane z-matrix is Td, got " + methane_pg.point_group);

  expect(throws_error([] {
           const std::vector<ZMatrixAtom> rows = {
               {6, 0, 0, 0, 0.0, 0.0, 0.0},
               {1, 5, 0, 0, 1.0, 0.0, 0.0},
           };
           (void)molsymm::zmatrix_positions(rows);
         }),
         "out-of-range bond reference throws");
  expect(throws_error([] {
           const std::vector<ZMatrixAtom> rows = {
               {6, 0, 0, 0, 0.0, 0.0, 0.0},
               {1, 1, 0, 0, 1.0, 0.0, 0.0},
               {1, 1, 3, 0, 1.0, 0.5 * kPi, 0.0},
           };
           (void)molsymm::zmatrix_positions(rows);
         }),
         "forward angle reference throws");
  expect(throws_error([] {
           const std::vector<ZMatrixAtom> rows = {
               {6, 0, 0, 0, 0.0, 0.0, 0.0},
               {1, 1, 0, 0, 0.0, 0.0, 0.0},
           };
           (void)molsymm::zmatrix_positions(rows);
         }),
         "zero bond length throws");
  expect(throws_error([] {
           const std::vector<ZMatrixAtom> rows = {
               {6, 0, 0, 0, 0.0, 0.0, 0.0},
               {1, 1, 0, 0, 1.0, 0.0, 0.0},
               {1, 1, 2, 0, 1.0, kPi, 0.0},
               {1, 1, 2, 3, 1.0, 0.5 * kPi, 0.0},
           };
           (void)molsymm::zmatrix_positions(rows);
         }),
         "collinear dihedral references throw");
  expect(throws_error([] {
           const std::vector<ZMatrixAtom> rows = {
               {6, 0, 0, 0, 0.0, 0.0, 0.0},
               {1, 1, 0, 0, 1.0, 0.0, 0.0},
               {1, 1, 2, 0, 1.0, 0.0, 0.0},
               {1, 2, 3, 1, 1.0, 0.5 * kPi, 0.0},
           };
           (void)molsymm::zmatrix_positions(rows);
         }),
         "zero-length bond reference throws");
}

} // namespace

int main() {
  try {
    test_point_groups();
    test_subgroups_and_table();
    test_displacements_water();
    test_displacements_methane();
    test_inertia_frame_stable();
    test_high_symmetry_counts();
    test_xyz_corpus();
    test_xyz_slight_rattle();
    test_xyz_large_rattle();
    test_zmatrix();
  } catch (const std::exception& ex) {
    std::cerr << "exception: " << ex.what() << "\n";
    return 1;
  }
  if (g_failed) {
    std::cerr << g_failed << " failures (" << g_passed << " checks passed)\n";
    return 1;
  }
  std::cout << "all tests passed (" << g_passed << " checks)\n";
  return 0;
}
