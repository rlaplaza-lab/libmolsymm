#pragma once

#include "molsymm/analysis.hpp"

#include <string>
#include <vector>

namespace molsymm {

struct WyckoffSite {
  std::string label;         // "6a", "2b"
  std::string site_symmetry; // stabilizer, e.g. "C4v"
  int multiplicity = 1;
  int degrees_of_freedom = 3; // 0 origin, 1 axis, 2 mirror, 3 generic
  Mat3 projector = Mat3::Identity();
  std::vector<Vec3> chamber; // normals n with n·p > 0 in the fundamental chamber
};

struct WyckoffOrbit {
  int atomic_number = 0; // 0 when the plan came from n_sites only
  std::string label;
  int multiplicity = 1;
  int degrees_of_freedom = 3;
};

struct WyckoffPlan {
  std::string point_group;
  int order = 1;
  std::vector<WyckoffOrbit> orbits;

  [[nodiscard]] Analysis assemble(const std::vector<Vec3>& au_positions) const;
};

std::vector<WyckoffSite> wyckoff_positions(const std::string& point_group);

std::vector<WyckoffPlan> wyckoff_plans(const std::string& point_group,
                                       const std::vector<int>& atomic_numbers);
std::vector<WyckoffPlan> wyckoff_plans(const std::string& point_group, int n_sites);

WyckoffPlan wyckoff_plan(const std::string& point_group,
                         const std::vector<std::string>& site_labels,
                         const std::vector<int>& atomic_numbers = {});

Analysis assemble_from_au(const std::string& point_group,
                          const std::vector<int>& representative_numbers,
                          const std::vector<Vec3>& au_positions,
                          const std::vector<std::string>& site_labels);

} // namespace molsymm
