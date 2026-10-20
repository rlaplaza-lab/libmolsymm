#include "molsymm/zmatrix.hpp"

#include "molsymm/error.hpp"

#include <Eigen/Geometry>

#include <cmath>
#include <string>

namespace molsymm {
namespace {

constexpr double kLengthEps = 1e-15;

void require_earlier(int ref, int current, int n, const char* role) {
  if (ref < 1 || ref > n) {
    throw Error(std::string(role) + " reference is out of range");
  }
  if (ref >= current) {
    throw Error(std::string(role) + " reference must refer to an earlier atom");
  }
}

const Vec3& placed(const std::vector<Vec3>& pos, int one_based) {
  return pos[static_cast<std::size_t>(one_based - 1)];
}

} // namespace

std::vector<Vec3> zmatrix_positions(std::span<const ZMatrixAtom> atoms) {
  const int n = static_cast<int>(atoms.size());
  if (n <= 0) {
    throw Error("z-matrix is empty");
  }

  std::vector<Vec3> pos(static_cast<std::size_t>(n), Vec3::Zero());
  for (int i = 0; i < n; ++i) {
    const ZMatrixAtom& atom = atoms[static_cast<std::size_t>(i)];
    const int current = i + 1;
    if (i == 0) {
      continue;
    }

    require_earlier(atom.bond_to, current, n, "bond");
    if (!(atom.distance > 0.0)) {
      throw Error("distance must be positive");
    }
    const Vec3& a = placed(pos, atom.bond_to);

    if (i == 1) {
      pos[static_cast<std::size_t>(i)] = Vec3(a.x(), a.y() + atom.distance, a.z());
      continue;
    }

    require_earlier(atom.angle_to, current, n, "angle");
    if (atom.angle_to == atom.bond_to) {
      throw Error("angle reference must differ from the bond reference");
    }
    const Vec3& b = placed(pos, atom.angle_to);

    if (i == 2) {
      const double ysign = (b.y() < a.y()) ? -1.0 : 1.0;
      pos[static_cast<std::size_t>(i)] =
          Vec3(a.x() + atom.distance * std::abs(std::sin(atom.angle)),
               a.y() + ysign * atom.distance * std::cos(atom.angle),
               a.z());
      continue;
    }

    require_earlier(atom.dihedral_to, current, n, "dihedral");
    if (atom.dihedral_to == atom.bond_to || atom.dihedral_to == atom.angle_to) {
      throw Error("dihedral reference must differ from the bond and angle references");
    }
    const Vec3& c = placed(pos, atom.dihedral_to);

    Vec3 ab = b - a;
    const double ab2 = ab.squaredNorm();
    if (ab2 < kLengthEps) {
      throw Error("bond reference has zero length");
    }
    ab /= std::sqrt(ab2);

    const Vec3 bc = b - c;
    Vec3 perp = ab.cross(bc);
    const double perp2 = perp.squaredNorm();
    if (perp2 < kLengthEps) {
      throw Error("dihedral references are collinear");
    }
    perp /= std::sqrt(perp2);

    const Vec3 along = atom.distance * ab;
    const Vec3 bent = Eigen::AngleAxisd(-atom.angle, perp) * along;
    pos[static_cast<std::size_t>(i)] = Eigen::AngleAxisd(atom.dihedral, ab) * bent + a;
  }
  return pos;
}

Molecule molecule_from_zmatrix(std::span<const ZMatrixAtom> atoms) {
  const std::vector<Vec3> pos = zmatrix_positions(atoms);
  Molecule mol;
  mol.atoms.resize(atoms.size());
  for (std::size_t i = 0; i < atoms.size(); ++i) {
    mol.atoms[i].atomic_number = atoms[i].atomic_number;
    mol.atoms[i].position = pos[i];
  }
  mol.resolve_masses();
  return mol;
}

} // namespace molsymm
