#pragma once

#include "molsymm/analysis.hpp"

#include <cstdint>
#include <span>
#include <string>
#include <vector>

namespace molsymm {

struct FragmentSite {
  int atom = 0;
  Mat3 projector = Mat3::Identity(); // analysis frame
  std::vector<Vec3> free_axes;       // original-frame +1 eigenspace
};

struct FragmentCopy {
  std::vector<int> atoms; // disjoint partition of the molecule
};

struct Fragmentation {
  std::vector<int> representatives; // one atom per orbit (lowest index)
  std::vector<FragmentSite> sites;
  std::vector<FragmentCopy> copies;
  std::vector<std::uint8_t> mask; // 1 on representatives
  std::vector<int> atomic_numbers;
  std::vector<double> masses;
  Vec3 center_of_mass = Vec3::Zero();
  Mat3 orientation = Mat3::Identity();
  std::vector<Operation> operations;

  [[nodiscard]] Vec3 frame_to_input(const Vec3& p) const {
    return orientation.transpose() * p + center_of_mass;
  }

  [[nodiscard]] Vec3 input_to_frame(const Vec3& p) const {
    return orientation * (p - center_of_mass);
  }
};

struct AssemblyPiece {
  int parent_copy = 0;
  std::vector<int> atoms;
};

struct SubgroupAssembly {
  std::string subgroup;
  std::vector<AssemblyPiece> pieces;
  std::vector<FragmentSite> sites;
  std::vector<int> atomic_numbers;
  std::vector<double> masses;
  Vec3 center_of_mass = Vec3::Zero();
  Mat3 orientation = Mat3::Identity();
  std::vector<Operation> operations;

  [[nodiscard]] Vec3 frame_to_input(const Vec3& p) const {
    return orientation.transpose() * p + center_of_mass;
  }

  [[nodiscard]] Vec3 input_to_frame(const Vec3& p) const {
    return orientation * (p - center_of_mass);
  }
};

Fragmentation fragmentation(const Analysis& analysis);
Molecule reconstruct(const Fragmentation& frag, std::span<const Vec3> representative_positions);

SubgroupAssembly subgroup_assembly(const Analysis& analysis, const std::string& subgroup);
Molecule assemble(const SubgroupAssembly& assembly,
                  const std::vector<std::vector<Vec3>>& piece_positions);

} // namespace molsymm
