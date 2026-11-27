// Deep subgroup walks on corpus geometries (Cubane Oh, C60 Ih) with full
// character-table checks. Complements the synthetic audit's orbit molecules.
#include "molsymm/molsymm.hpp"

#include "fixtures.hpp"
#include "table_checks.hpp"

#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef MOLSYMM_TEST_DATA_DIR
#error "MOLSYMM_TEST_DATA_DIR must be defined"
#endif

using molsymm::Analysis;

namespace {

int g_fail = 0;
int g_checked = 0;

void fail(const std::string& msg) {
  ++g_fail;
  std::cerr << "FAIL: " << msg << "\n";
}

void check_table(const Analysis& a, const std::string& label) {
  ++g_checked;
  molsymm::test::check_character_table(a, label, fail);
}

molsymm::Molecule load_xyz(const std::string& path) {
  static const std::map<std::string, int> z_of{{"H", 1}, {"C", 6}};
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("cannot open " + path);
  }
  int n = 0;
  in >> n;
  std::string line;
  std::getline(in, line);
  std::getline(in, line);
  std::vector<int> z;
  std::vector<double> xyz;
  for (int i = 0; i < n; ++i) {
    std::string sym;
    double x = 0;
    double y = 0;
    double zz = 0;
    in >> sym >> x >> y >> zz;
    z.push_back(z_of.at(sym));
    xyz.insert(xyz.end(), {x, y, zz});
  }
  return molsymm::Molecule(z, xyz);
}

void walk(const Analysis& a, const std::string& path, int depth) {
  for (const std::string& sg : a.subgroups) {
    const std::string label = path + "->" + sg;
    try {
      Analysis s = a;
      s.select_subgroup(sg);
      if (s.point_group != sg) {
        fail(label + " produced " + s.point_group);
      }
      check_table(s, label);
      if (depth < 2 && sg != a.point_group && sg != "C1") {
        walk(s, label, depth + 1);
      }
    } catch (const std::exception& e) {
      fail(label + " threw: " + e.what());
    }
  }
}

} // namespace

void test_subgroup_probe() {
  g_fail = 0;
  g_checked = 0;
  const std::string root = std::string(MOLSYMM_TEST_DATA_DIR) + "/iqmol/Molecules";

  struct Case {
    std::string name{};
    molsymm::Molecule mol{};
    double tol = 0;
  };
  const std::vector<Case> cases{
      {"Cubane", load_xyz(root + "/Cyclic_Alkanes/Cubane.xyz"), 1e-3},
      {"C60", load_xyz(root + "/Fullerenes/Buckminsterfullerene.xyz"), 1e-4},
  };

  for (const Case& c : cases) {
    molsymm::Tolerances tol;
    tol.primary = c.tol;
    Analysis an = molsymm::analyze(c.mol, tol);
    check_table(an, c.name);
    walk(an, c.name, 0);
  }

  // Fixture Oh geometry: targeted axial chains that previously broke tables.
  {
    Analysis an = molsymm::analyze(molsymm::test::sf6());
    for (const std::vector<std::string>& chain :
         std::vector<std::vector<std::string>>{{"C2v"}, {"D4"}, {"Th", "C2v"}, {"Td", "C2v"}}) {
      Analysis s2 = an;
      std::string lab = "SF6";
      try {
        for (const auto& g : chain) {
          s2.select_subgroup(g);
          lab += "->" + g;
        }
        check_table(s2, lab);
      } catch (const std::exception& e) {
        fail(lab + " threw: " + e.what());
      }
    }
  }

  if (g_fail != 0) {
    throw std::runtime_error(std::to_string(g_fail) + " subgroup-probe failures (" +
                             std::to_string(g_checked) + " realizations checked)");
  }
  std::cout << "subgroup-probe checks passed (" << g_checked << " realizations)\n";
}
