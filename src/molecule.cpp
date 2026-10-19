#include "molsymm/molecule.hpp"

#include "molsymm/error.hpp"
#include "molsymm/masses.hpp"

namespace molsymm {

Molecule::Molecule(std::span<const int> atomic_numbers, std::span<const double> xyz_flat) {
  if (xyz_flat.size() != atomic_numbers.size() * 3) {
    throw Error("xyz length must be 3 * number of atoms");
  }
  atoms.resize(atomic_numbers.size());
  for (std::size_t i = 0; i < atomic_numbers.size(); ++i) {
    atoms[i].atomic_number = atomic_numbers[i];
    atoms[i].position = Vec3(xyz_flat[3 * i], xyz_flat[3 * i + 1], xyz_flat[3 * i + 2]);
  }
  resolve_masses();
}

Molecule::Molecule(std::span<const int> atomic_numbers,
                   const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>& xyz) {
  if (xyz.cols() != static_cast<Eigen::Index>(atomic_numbers.size())) {
    throw Error("xyz columns must match number of atoms");
  }
  atoms.resize(atomic_numbers.size());
  for (std::size_t i = 0; i < atomic_numbers.size(); ++i) {
    atoms[i].atomic_number = atomic_numbers[i];
    atoms[i].position = xyz.col(static_cast<Eigen::Index>(i));
  }
  resolve_masses();
}

void Molecule::resolve_masses() {
  for (Atom& a : atoms) {
    if (a.mass <= 0.0) {
      a.mass = atomic_mass(a.atomic_number);
      if (a.mass <= 0.0) {
        a.mass = 1.0;
      }
    }
  }
}

Vec3 center_of_mass(const Molecule& mol) {
  double mtot = 0.0;
  Vec3 com = Vec3::Zero();
  for (const Atom& a : mol.atoms) {
    com += a.mass * a.position;
    mtot += a.mass;
  }
  if (mtot <= 0.0) {
    throw Error("total mass is zero");
  }
  return com / mtot;
}

Mat3 inertia_tensor(const Molecule& mol, const Vec3& origin) {
  Mat3 I = Mat3::Zero();
  for (const Atom& a : mol.atoms) {
    const Vec3 r = a.position - origin;
    I += a.mass * ((r.squaredNorm() * Mat3::Identity()) - r * r.transpose());
  }
  return I;
}

} // namespace molsymm
