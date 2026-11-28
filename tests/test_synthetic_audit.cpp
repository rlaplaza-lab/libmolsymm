// Synthetic point-group audit: listed vs enumerated subgroups, character tables
// on every realization, cubic orientation robustness, axis-free T/Th detection.
// Also tautology from H2/CO: apply generators stepwise and re-detect each step.
#include "fixtures.hpp"
#include "group_gen.hpp"
#include "molsymm/molsymm.hpp"
#include "table_checks.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

using molsymm::Analysis;
using molsymm::Mat3;
using molsymm::Molecule;
using molsymm::Vec3;
using molsymm::test::close_group;
using molsymm::test::co;
using molsymm::test::h2;
using molsymm::test::kPi;
using molsymm::test::mat_near;
using molsymm::test::molecule_of;
using molsymm::test::refl;
using molsymm::test::rot;

namespace {

int g_fail = 0;

void fail(const std::string& msg) {
  ++g_fail;
  std::cerr << "FAIL: " << msg << "\n";
}

void check_table(const Analysis& a, const std::string& label) {
  molsymm::test::check_character_table(a, label, fail);
}

Molecule axis_orbit(const Molecule& axis, const std::vector<Mat3>& g, std::size_t norb) {
  const std::array<Vec3, 4> seeds{Vec3(1.13, 0.37, 0.71),
                                  Vec3(0.52, 1.41, -0.33),
                                  Vec3(-0.87, 0.29, 1.52),
                                  Vec3(0.41, -1.07, -1.23)};
  const std::array<int, 4> zs{9, 17, 35, 53};
  std::vector<int> z;
  std::vector<double> xyz;
  for (const auto& a : axis.atoms) {
    z.push_back(a.atomic_number);
    xyz.insert(xyz.end(), {a.position.x(), a.position.y(), a.position.z()});
  }
  for (std::size_t s = 0; s < norb; ++s) {
    std::vector<Vec3> orbit;
    for (const Mat3& m : g) {
      const Vec3 p = m * seeds[s];
      if (std::none_of(
              orbit.begin(), orbit.end(), [&](const Vec3& q) { return (q - p).norm() < 1e-8; })) {
        orbit.push_back(p);
      }
    }
    for (const Vec3& p : orbit) {
      z.push_back(zs[s]);
      xyz.insert(xyz.end(), {p.x(), p.y(), p.z()});
    }
  }
  return Molecule(z, xyz);
}

bool matrices_permute(const Molecule& mol, const std::vector<Mat3>& g, double tol) {
  const int n = mol.size();
  for (const Mat3& R : g) {
    std::vector<char> claimed(static_cast<std::size_t>(n), 0);
    for (int i = 0; i < n; ++i) {
      const Vec3 rp = R * mol.atoms[static_cast<std::size_t>(i)].position;
      const int zi = mol.atoms[static_cast<std::size_t>(i)].atomic_number;
      int best = -1;
      double best_d = tol;
      for (int j = 0; j < n; ++j) {
        if (claimed[static_cast<std::size_t>(j)] ||
            mol.atoms[static_cast<std::size_t>(j)].atomic_number != zi) {
          continue;
        }
        const double d = (rp - mol.atoms[static_cast<std::size_t>(j)].position).norm();
        if (d <= best_d) {
          best_d = d;
          best = j;
        }
      }
      if (best < 0) {
        return false;
      }
      claimed[static_cast<std::size_t>(best)] = 1;
    }
  }
  return true;
}

using Bits = std::array<std::uint64_t, 2>;

std::set<Bits> enumerate_subgroups(const std::vector<Mat3>& g) {
  const int h = static_cast<int>(g.size());
  std::vector<std::vector<int>> mul(static_cast<std::size_t>(h), std::vector<int>(h, -1));
  int inv_idx = -1;
  for (int i = 0; i < h; ++i) {
    if (mat_near(g[static_cast<std::size_t>(i)], -Mat3::Identity())) {
      inv_idx = i;
    }
    for (int j = 0; j < h; ++j) {
      const Mat3 p = g[static_cast<std::size_t>(i)] * g[static_cast<std::size_t>(j)];
      for (int k = 0; k < h; ++k) {
        if (mat_near(g[static_cast<std::size_t>(k)], p)) {
          mul[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = k;
          break;
        }
      }
    }
  }
  auto gen = [&](const std::vector<int>& gens) {
    std::vector<char> in(static_cast<std::size_t>(h), 0);
    std::vector<int> elems{0};
    in[0] = 1;
    for (std::size_t q = 0; q < elems.size(); ++q) {
      for (int s : gens) {
        const int k = mul[static_cast<std::size_t>(elems[q])][static_cast<std::size_t>(s)];
        if (!in[static_cast<std::size_t>(k)]) {
          in[static_cast<std::size_t>(k)] = 1;
          elems.push_back(k);
        }
      }
    }
    Bits b{0, 0};
    for (int e : elems) {
      b[static_cast<std::size_t>(e / 64)] |= (std::uint64_t{1} << static_cast<unsigned>(e % 64));
    }
    return b;
  };
  std::set<Bits> out;
  const bool triples = h <= 48;
  for (int i = 0; i < h; ++i) {
    for (int j = i; j < h; ++j) {
      out.insert(gen({i, j}));
      if (triples) {
        for (int k = j; k < h; ++k) {
          out.insert(gen({i, j, k}));
        }
      } else if (inv_idx >= 0) {
        out.insert(gen({i, j, inv_idx}));
      }
    }
  }
  return out;
}

struct Parent {
  std::string name{};
  std::vector<Mat3> gens{};
};

std::vector<Parent> parents() {
  std::vector<Parent> p;
  const Vec3 z = Vec3::UnitZ();
  const Mat3 c2x = rot(Vec3::UnitX(), kPi);
  const Mat3 sh = refl(z);
  const Mat3 sv = refl(Vec3::UnitY());
  const Mat3 inv = -Mat3::Identity();
  p.push_back({"C1", {}});
  p.push_back({"Cs", {sh}});
  p.push_back({"Ci", {inv}});
  for (int n = 2; n <= 8; ++n) {
    const std::string s = std::to_string(n);
    const Mat3 cn = rot(z, 2.0 * kPi / n);
    const Mat3 s2n = sh * rot(z, kPi / n);
    p.push_back({"C" + s, {cn}});
    p.push_back({"C" + s + "v", {cn, sv}});
    p.push_back({"C" + s + "h", {cn, sh}});
    p.push_back({"D" + s, {cn, c2x}});
    p.push_back({"D" + s + "h", {cn, c2x, sh}});
    p.push_back({"D" + s + "d", {s2n, c2x}});
  }
  for (int n : {4, 6, 8, 10, 12}) {
    p.push_back({"S" + std::to_string(n), {sh * rot(z, 2.0 * kPi / n)}});
  }
  const Mat3 c3 = rot(Vec3(1, 1, 1), 2.0 * kPi / 3.0);
  const Mat3 c2z = rot(z, kPi);
  const Mat3 c4z = rot(z, 0.5 * kPi);
  const double phi = 0.5 * (1.0 + std::sqrt(5.0));
  const Mat3 c5 = rot(Vec3(0, 1, phi), 2.0 * kPi / 5.0);
  p.push_back({"T", {c3, c2z}});
  p.push_back({"Th", {c3, c2z, inv}});
  p.push_back({"Td", {c3, sh * c4z}});
  p.push_back({"O", {c3, c4z}});
  p.push_back({"Oh", {c3, c4z, inv}});
  p.push_back({"I", {c3, c5}});
  p.push_back({"Ih", {c3, c5, inv}});
  return p;
}

std::string join(const std::set<std::string>& s) {
  std::string o;
  for (const auto& x : s) {
    o += (o.empty() ? "" : " ") + x;
  }
  return o;
}

void expect_step(const std::string& tag,
                 const Molecule& mol,
                 const std::vector<Mat3>& planted,
                 const std::string& expect) {
  Analysis a;
  try {
    a = molsymm::analyze(mol);
  } catch (const std::exception& e) {
    fail(tag + " analyze threw: " + e.what());
    return;
  }
  if (a.point_group != expect) {
    fail(tag + " detected as " + a.point_group + " (want " + expect + ")");
    return;
  }
  if (!a.infinite && !matrices_permute(mol, planted, 1e-6)) {
    fail(tag + " planted ops do not permute atoms");
  }
  if (!a.infinite && a.order != static_cast<int>(planted.size())) {
    fail(tag + " order " + std::to_string(a.order) + " != planted " +
         std::to_string(planted.size()));
  }
}

void tautology_from_infinity() {
  const Molecule axis_h2 = h2();
  const Molecule axis_co = co();
  if (molsymm::analyze(axis_h2).point_group != "Dinfh") {
    fail("H2 axis not Dinfh");
  }
  if (molsymm::analyze(axis_co).point_group != "Cinfv") {
    fail("CO axis not Cinfv");
  }

  const Vec3 z = Vec3::UnitZ();
  const Mat3 c2x = rot(Vec3::UnitX(), kPi);
  const Mat3 sh = refl(z);
  const Mat3 sv = refl(Vec3::UnitY());
  const Mat3 inv = -Mat3::Identity();

  struct Step {
    Mat3 gen{};
    std::string expect{};
  };
  struct Sc {
    std::string name{};
    const Molecule* axis = nullptr;
    std::vector<Step> steps{};
  };

  std::vector<Sc> scs{
      {"C1", &axis_co, {}}, {"Cs", &axis_h2, {{sh, "Cs"}}}, {"Ci", &axis_h2, {{inv, "Ci"}}}};
  for (int n = 2; n <= 12; ++n) {
    const std::string s = std::to_string(n);
    const Mat3 cn = rot(z, 2.0 * kPi / n);
    const Mat3 s2n = sh * rot(z, kPi / n);
    scs.push_back({"C" + s, &axis_co, {{cn, "C" + s}}});
    scs.push_back({"C" + s + "v", &axis_co, {{cn, "C" + s}, {sv, "C" + s + "v"}}});
    scs.push_back({"C" + s + "h", &axis_h2, {{cn, "C" + s}, {sh, "C" + s + "h"}}});
    scs.push_back({"D" + s, &axis_h2, {{cn, "C" + s}, {c2x, "D" + s}}});
    scs.push_back({"D" + s + "h", &axis_h2, {{cn, "C" + s}, {c2x, "D" + s}, {sh, "D" + s + "h"}}});
    scs.push_back(
        {"D" + s + "d", &axis_h2, {{s2n, "S" + std::to_string(2 * n)}, {c2x, "D" + s + "d"}}});
  }
  for (int n : {4, 6, 8, 10, 12}) {
    scs.push_back({"S" + std::to_string(n),
                   &axis_h2,
                   {{sh * rot(z, 2.0 * kPi / n), "S" + std::to_string(n)}}});
  }

  for (const Sc& sc : scs) {
    std::vector<Mat3> gens;
    if (sc.steps.empty()) {
      const auto planted = close_group(std::vector<Mat3>{});
      expect_step(sc.name, axis_orbit(*sc.axis, planted, 2), planted, sc.name);
      continue;
    }
    for (std::size_t i = 0; i < sc.steps.size(); ++i) {
      gens.push_back(sc.steps[i].gen);
      const auto planted = close_group(gens);
      expect_step(sc.name + "[" + std::to_string(i + 1) + "->" + sc.steps[i].expect + "]",
                  axis_orbit(*sc.axis, planted, planted.size() <= 12 ? 4 : 2),
                  planted,
                  sc.steps[i].expect);
    }
  }
}

} // namespace

void test_synthetic_audit() {
  g_fail = 0;
  std::mt19937 rng(12345); // NOLINT(cert-msc32-c,cert-msc51-cpp) — reproducible audit
  auto arrow = [](std::string a, const std::string& b) {
    a += "->";
    a += b;
    return a;
  };
  for (const Parent& par : parents()) {
    const std::vector<Mat3> g = close_group(par.gens);
    Analysis a;
    try {
      a = molsymm::analyze(molecule_of(g, Mat3::Identity(), nullptr));
    } catch (const std::exception& e) {
      fail(par.name + " analyze threw: " + e.what());
      continue;
    }
    if (a.point_group != par.name) {
      fail(par.name + " detected as " + a.point_group);
      continue;
    }
    check_table(a, par.name);

    const std::set<std::string> listed(a.subgroups.begin(), a.subgroups.end());
    for (const std::string& sg : a.subgroups) {
      try {
        Analysis s = a;
        s.select_subgroup(sg);
        check_table(s, arrow(par.name, sg));
        for (const std::string& sg2 : s.subgroups) {
          try {
            Analysis s2 = s;
            s2.select_subgroup(sg2);
            check_table(s2, arrow(arrow(par.name, sg), sg2));
          } catch (const std::exception& e) {
            fail(arrow(arrow(par.name, sg), sg2) + " threw: " + e.what());
          }
        }
      } catch (const std::exception& e) {
        fail(arrow(par.name, sg) + " threw: " + e.what());
      }
    }

    std::set<std::string> actual;
    for (const Bits& b : enumerate_subgroups(g)) {
      std::vector<Mat3> hsub;
      for (std::size_t e = 0; e < g.size(); ++e) {
        if ((b[e / 64] >> (e % 64)) & 1U) {
          hsub.push_back(g[e]);
        }
      }
      try {
        const Analysis sa = molsymm::analyze(molecule_of(hsub, Mat3::Identity(), nullptr));
        if (sa.order != static_cast<int>(hsub.size())) {
          fail(par.name + " subgroup of order " + std::to_string(hsub.size()) + " detected as " +
               sa.point_group);
        }
        actual.insert(sa.point_group);
      } catch (const std::exception& e) {
        fail(par.name + " subgroup analyze threw: " + e.what());
      }
    }
    std::set<std::string> missing;
    std::set<std::string> extra;
    std::set_difference(actual.begin(),
                        actual.end(),
                        listed.begin(),
                        listed.end(),
                        std::inserter(missing, missing.end()));
    std::set_difference(listed.begin(),
                        listed.end(),
                        actual.begin(),
                        actual.end(),
                        std::inserter(extra, extra.end()));
    if (!missing.empty()) {
      fail(par.name + " list is missing: " + join(missing));
    }
    if (!extra.empty()) {
      fail(par.name + " lists non-subgroups: " + join(extra));
    }

    if (par.name == "T" || par.name == "Th" || par.name == "Td" || par.name == "O" ||
        par.name == "Oh") {
      std::uniform_real_distribution<double> u(-1.0, 1.0);
      for (int trial = 0; trial < 5; ++trial) {
        const Mat3 frame = rot(Vec3(u(rng), u(rng), u(rng)), kPi * u(rng));
        const std::string tag = par.name + "[rot " + std::to_string(trial) + "]";
        try {
          Analysis r = molsymm::analyze(molecule_of(g, frame, &rng));
          if (r.point_group != par.name) {
            fail(tag + " detected as " + r.point_group);
            continue;
          }
          check_table(r, tag);
          for (const char* sg : {"T", "Th", "Td", "C2v", "D2", "D2h", "C2h", "D4"}) {
            if (std::find(r.subgroups.begin(), r.subgroups.end(), sg) == r.subgroups.end()) {
              continue;
            }
            Analysis s = r;
            s.select_subgroup(sg);
            check_table(s, tag + "->" + sg);
          }
        } catch (const std::exception& e) {
          fail(tag + " threw: " + e.what());
        }
      }
    }
  }

  for (const Parent& par : parents()) {
    const std::vector<Mat3> g = close_group(par.gens);
    try {
      const Analysis a = molsymm::analyze(molecule_of(g, Mat3::Identity(), nullptr, false));
      if (a.point_group != par.name) {
        fail("generic-orbit " + par.name + " detected as " + a.point_group);
      }
    } catch (const std::exception& e) {
      fail("generic-orbit " + par.name + " threw: " + e.what());
    }
  }
  {
    std::vector<int> z{26};
    std::vector<double> xyz{0, 0, 0};
    const double ro = 2.0;
    const double a = 2.6;
    const double b = 0.8;
    for (double sgn : {1.0, -1.0}) {
      const std::vector<std::array<double, 3>> o{{ro, 0, 0}, {0, ro, 0}, {0, 0, ro}};
      const std::vector<std::array<double, 3>> hh{
          {a, b, 0}, {a, -b, 0}, {0, a, b}, {0, a, -b}, {b, 0, a}, {-b, 0, a}};
      for (const auto& p : o) {
        z.push_back(8);
        xyz.insert(xyz.end(), {sgn * p[0], sgn * p[1], sgn * p[2]});
      }
      for (const auto& p : hh) {
        z.push_back(1);
        xyz.insert(xyz.end(), {sgn * p[0], sgn * p[1], sgn * p[2]});
      }
    }
    const Analysis a2 = molsymm::analyze(molsymm::Molecule(z, xyz));
    if (a2.point_group != "Th") {
      fail("Th hexaaqua detected as " + a2.point_group);
    }
  }

  tautology_from_infinity();

  if (g_fail != 0) {
    throw std::runtime_error(std::to_string(g_fail) + " synthetic-audit failures");
  }
  std::cout << "synthetic-audit checks passed\n";
}
