#include "molsymm/molsymm.hpp"

#include <algorithm>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
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

const char* kFiniteGroups[] = {
    "C1",  "Cs",  "Ci",  "C2",  "C3",  "C4",  "C5",  "C6",  "C7",  "C8",  "C2v", "C3v",
    "C4v", "C5v", "C6v", "C2h", "C3h", "C4h", "C5h", "C6h", "D2",  "D3",  "D4",  "D5",
    "D6",  "D2h", "D3h", "D4h", "D5h", "D6h", "D2d", "D3d", "D4d", "D5d", "D6d", "S4",
    "S6",  "S8",  "T",   "Td",  "Th",  "O",   "Oh",  "I",   "Ih",
};

void test_catalog_and_generic() {
  for (const char* name : kFiniteGroups) {
    std::vector<molsymm::WyckoffSite> sites;
    try {
      sites = molsymm::wyckoff_positions(name);
    } catch (const std::exception& e) {
      expect(false, std::string(name) + " positions: " + e.what());
      continue;
    }
    expect(!sites.empty(), std::string(name) + " has sites");
    int g_order = 1;
    for (const auto& s : sites) {
      expect(s.multiplicity > 0, std::string(name) + " " + s.label + " mult>0");
      g_order = std::max(g_order, s.multiplicity);
    }
    const molsymm::WyckoffSite* generic = nullptr;
    for (const auto& s : sites) {
      if (s.degrees_of_freedom == 3) {
        generic = &s;
      }
      expect(g_order % s.multiplicity == 0, std::string(name) + " " + s.label + " divides |G|");
    }
    expect(generic != nullptr, std::string(name) + " has generic site");
    if (generic == nullptr) {
      continue;
    }
    try {
      const std::vector<molsymm::Vec3> seeds{
          {0.91, 0.37, 0.53}, {-0.41, 1.07, -0.29}, {0.52, -0.83, 1.19}, {-1.13, 0.29, -0.71}};
      const std::vector<int> zs{6, 7, 8, 9};
      const int norb = (generic->multiplicity <= 12) ? 4 : 2;
      std::vector<int> numbers(zs.begin(), zs.begin() + norb);
      std::vector<molsymm::Vec3> pos(seeds.begin(), seeds.begin() + norb);
      std::vector<std::string> labels(static_cast<std::size_t>(norb), generic->label);
      auto a = molsymm::assemble_from_au(name, numbers, pos, labels);
      expect(a.point_group == name, std::string(name) + " stamped group");
      expect(a.order == generic->multiplicity, std::string(name) + " order");
      expect(a.molecule.size() == norb * generic->multiplicity, std::string(name) + " atom count");
      const auto detected = molsymm::analyze(a.molecule).point_group;
      expect(detected == name, std::string(name) + " analyze->" + detected);
    } catch (const std::exception& e) {
      expect(false, std::string(name) + " generic assemble: " + e.what());
    }
  }
}

void test_oh_plans() {
  auto sf6 = molsymm::wyckoff_plans("Oh", std::vector<int>{16, 9, 9, 9, 9, 9, 9});
  expect(sf6.size() == 1, "Oh SF6 one plan");
  if (!sf6.empty()) {
    std::set<std::string> labels;
    for (const auto& o : sf6[0].orbits) {
      labels.insert(o.label);
    }
    expect(labels.contains("1a") && labels.contains("6a"), "Oh SF6 {1a,6a}");
  }

  auto p12 = molsymm::wyckoff_plans("Oh", 12);
  expect(p12.size() >= 2, "Oh n=12 has >=2 plans");
  bool has_12a = false;
  bool has_6a6a = false;
  for (const auto& p : p12) {
    if (p.orbits.size() == 1 && p.orbits[0].label == "12a") {
      has_12a = true;
    }
    if (p.orbits.size() == 2 && p.orbits[0].label == "6a" && p.orbits[1].label == "6a") {
      has_6a6a = true;
    }
  }
  expect(has_12a, "Oh n=12 includes {12a}");
  expect(has_6a6a, "Oh n=12 includes {6a,6a}");
  expect(molsymm::wyckoff_plans("Oh", 5).empty(), "Oh n=5 no plan");
}

void test_d2h_axes() {
  int n2 = 0;
  for (const auto& s : molsymm::wyckoff_positions("D2h")) {
    if (s.multiplicity == 2 && s.degrees_of_freedom == 1) {
      ++n2;
    }
  }
  expect(n2 == 3, "D2h has three multiplicity-2 axes");
}

void test_assemble_molecules() {
  {
    auto a =
        molsymm::assemble_from_au("C2v",
                                  {8, 1},
                                  {molsymm::Vec3(0, 0, 0.1173), molsymm::Vec3(0, 0.7572, -0.4692)},
                                  {"1a", "2b"});
    expect(a.point_group == "C2v" && a.molecule.size() == 3, "water assemble");
    expect(molsymm::analyze(a.molecule).point_group == "C2v", "water analyze");
    auto frag = molsymm::fragmentation(a);
    expect(frag.representatives.size() == 2, "water 2 reps");
    std::vector<molsymm::Vec3> reps;
    reps.reserve(frag.representatives.size());
    for (int r : frag.representatives) {
      reps.push_back(a.molecule.atoms[static_cast<std::size_t>(r)].position);
    }
    expect(molsymm::reconstruct(frag, reps).size() == 3, "water reconstruct");
    expect(!molsymm::subgroup_assembly(a, "Cs").pieces.empty(), "water subgroup_assembly");
  }

  {
    const molsymm::Vec3 h = molsymm::Vec3(1, 1, 1).normalized() * 1.09;
    auto a = molsymm::assemble_from_au("Td", {6, 1}, {molsymm::Vec3::Zero(), h}, {"1a", "4a"});
    expect(a.molecule.size() == 5, "methane n=5");
    expect(molsymm::analyze(a.molecule).point_group == "Td", "methane analyze");
    auto frag = molsymm::fragmentation(a);
    std::vector<molsymm::Vec3> reps;
    reps.reserve(frag.representatives.size());
    for (int r : frag.representatives) {
      reps.push_back(a.molecule.atoms[static_cast<std::size_t>(r)].position);
    }
    expect(molsymm::reconstruct(frag, reps).size() == 5, "methane reconstruct");
  }

  {
    auto a = molsymm::assemble_from_au(
        "Oh", {16, 9}, {molsymm::Vec3::Zero(), molsymm::Vec3(0, 0, 1.56)}, {"1a", "6a"});
    expect(a.molecule.size() == 7, "SF6 n=7");
    expect(molsymm::analyze(a.molecule).point_group == "Oh", "SF6 analyze");
    auto frag = molsymm::fragmentation(a);
    std::vector<molsymm::Vec3> reps;
    reps.reserve(frag.representatives.size());
    for (int r : frag.representatives) {
      reps.push_back(a.molecule.atoms[static_cast<std::size_t>(r)].position);
    }
    expect(molsymm::reconstruct(frag, reps).size() == 7, "SF6 reconstruct");
    expect(!molsymm::subgroup_assembly(a, "C4v").pieces.empty(), "SF6 subgroup_assembly");
  }
}

void test_errors() {
  bool threw = false;
  try {
    (void)molsymm::wyckoff_positions("Cinfv");
  } catch (const molsymm::Error&) {
    threw = true;
  }
  expect(threw, "Cinfv throws");

  threw = false;
  try {
    (void)molsymm::assemble_from_au(
        "C2v", {8, 1}, {molsymm::Vec3(0, 0, 0.1), molsymm::Vec3(0.5, 0.5, 0.5)}, {"1a", "2a"});
  } catch (const molsymm::Error&) {
    threw = true;
  }
  expect(threw, "off-stratum throws");

  threw = false;
  try {
    (void)molsymm::assemble_from_au(
        "Oh", {9, 9}, {molsymm::Vec3(0, 0, 1.0), molsymm::Vec3(0, 0, 1.0)}, {"6a", "6a"});
  } catch (const molsymm::Error&) {
    threw = true;
  }
  expect(threw, "coincident orbits throw");
}

} // namespace

void test_wyckoff() {
  test_catalog_and_generic();
  test_oh_plans();
  test_d2h_axes();
  test_assemble_molecules();
  test_errors();
  if (g_failed) {
    throw std::runtime_error(std::to_string(g_failed) + " wyckoff failures (" +
                             std::to_string(g_passed) + " checks passed)");
  }
  std::cout << "wyckoff checks passed (" << g_passed << ")\n";
}
