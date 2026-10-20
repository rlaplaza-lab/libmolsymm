#include "molsymm/molsymm.hpp"

#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <utility>
#include <vector>

namespace nb = nanobind;

namespace {

constexpr const char* kZEntryHint =
    "z-matrix entry must be (Z), (Z, bond_to, R), (Z, bond_to, R, angle_to, angle), "
    "or (Z, bond_to, R, angle_to, angle, dihedral_to, dihedral)";

molsymm::ZMatrixAtom zatom_from_python(nb::handle entry) {
  if (!nb::isinstance<nb::sequence>(entry) || nb::isinstance<nb::str>(entry)) {
    throw molsymm::Error(kZEntryHint);
  }
  nb::sequence seq = nb::cast<nb::sequence>(entry);
  const std::size_t n = nb::len(seq);
  if (n != 1 && n != 3 && n != 5 && n != 7) {
    throw molsymm::Error(kZEntryHint);
  }
  molsymm::ZMatrixAtom atom;
  atom.atomic_number = nb::cast<int>(seq[0]);
  if (n >= 3) {
    atom.bond_to = nb::cast<int>(seq[1]);
    atom.distance = nb::cast<double>(seq[2]);
  }
  if (n >= 5) {
    atom.angle_to = nb::cast<int>(seq[3]);
    atom.angle = nb::cast<double>(seq[4]);
  }
  if (n >= 7) {
    atom.dihedral_to = nb::cast<int>(seq[5]);
    atom.dihedral = nb::cast<double>(seq[6]);
  }
  return atom;
}

std::vector<molsymm::ZMatrixAtom> zrows_from_python(nb::sequence entries) {
  std::vector<molsymm::ZMatrixAtom> rows;
  rows.reserve(nb::len(entries));
  for (nb::handle entry : entries) {
    rows.push_back(zatom_from_python(entry));
  }
  return rows;
}

Eigen::MatrixXd positions_to_rows(const std::vector<molsymm::Vec3>& pos) {
  Eigen::MatrixXd xyz(static_cast<Eigen::Index>(pos.size()), 3);
  for (std::size_t i = 0; i < pos.size(); ++i) {
    xyz.row(static_cast<Eigen::Index>(i)) = pos[i].transpose();
  }
  return xyz;
}

} // namespace

NB_MODULE(molsymm, m) {
  m.doc() = "Molecular point-group symmetry";

  nb::class_<molsymm::Analysis>(m, "Analysis")
      .def_ro("point_group", &molsymm::Analysis::point_group)
      .def_ro("order", &molsymm::Analysis::order)
      .def_ro("symmetry_number", &molsymm::Analysis::symmetry_number)
      .def_ro("rmsd", &molsymm::Analysis::rmsd)
      .def_ro("subgroups", &molsymm::Analysis::subgroups)
      .def_prop_ro("operations",
                   [](const molsymm::Analysis& a) {
                     std::vector<Eigen::Matrix3d> mats;
                     mats.reserve(a.operations.size());
                     for (const auto& op : a.operations) {
                       mats.push_back(op.matrix);
                     }
                     return mats;
                   })
      .def_prop_ro("irrep_names",
                   [](const molsymm::Analysis& a) { return a.character_table.irrep_names; })
      .def("symmetrize",
           [](molsymm::Analysis& a) {
             const double rms = a.symmetrize();
             const int n = a.molecule.size();
             Eigen::MatrixXd xyz(n, 3);
             for (int i = 0; i < n; ++i) {
               xyz.row(i) = a.frame_to_input(a.molecule.atoms[static_cast<std::size_t>(i)].position)
                                .transpose();
             }
             return nb::make_tuple(rms, xyz);
           })
      .def("select_subgroup", &molsymm::Analysis::select_subgroup)
      .def(
          "displacements",
          [](const molsymm::Analysis& a, double step) {
            return molsymm::unique_displacements(a, step);
          },
          nb::arg("step") = 1e-3);

  auto mols_to_xyz = [](const std::vector<molsymm::Molecule>& mols) {
    std::vector<Eigen::MatrixXd> out;
    out.reserve(mols.size());
    for (const auto& mol : mols) {
      Eigen::MatrixXd xyz(mol.size(), 3);
      for (int i = 0; i < mol.size(); ++i) {
        xyz.row(i) = mol.atoms[static_cast<std::size_t>(i)].position.transpose();
      }
      out.push_back(std::move(xyz));
    }
    return out;
  };

  nb::class_<molsymm::DisplacementSet>(m, "DisplacementSet")
      .def_ro("step", &molsymm::DisplacementSet::step)
      .def_prop_ro("n",
                   [](const molsymm::DisplacementSet& d) {
                     return static_cast<int>(d.displacements.size());
                   })
      .def_prop_ro("plus",
                   [mols_to_xyz](const molsymm::DisplacementSet& d) { return mols_to_xyz(d.plus); })
      .def_prop_ro(
          "minus",
          [mols_to_xyz](const molsymm::DisplacementSet& d) { return mols_to_xyz(d.minus); })
      .def("gradient_from_energies",
           [](const molsymm::DisplacementSet& d,
              const std::vector<double>& ep,
              const std::vector<double>& em) { return d.gradient_from_energies(ep, em); })
      .def("hessian_from_gradients",
           [](const molsymm::DisplacementSet& d,
              const std::vector<Eigen::VectorXd>& gp,
              const std::vector<Eigen::VectorXd>& gm) { return d.hessian_from_gradients(gp, gm); });

  m.def(
      "analyze",
      [](const std::vector<int>& numbers, const Eigen::MatrixXd& positions, double primary) {
        if (positions.cols() != 3 && positions.rows() != 3) {
          throw molsymm::Error("positions must be (n,3) or (3,n)");
        }
        molsymm::Tolerances tol;
        tol.primary = primary;
        if (positions.cols() == 3) {
          Eigen::MatrixXd xyz = positions.transpose();
          if (static_cast<int>(numbers.size()) != xyz.cols()) {
            throw molsymm::Error("numbers/positions size mismatch");
          }
          molsymm::Molecule mol(numbers, xyz);
          return molsymm::analyze(mol, tol);
        }
        if (static_cast<int>(numbers.size()) != positions.cols()) {
          throw molsymm::Error("numbers/positions size mismatch");
        }
        molsymm::Molecule mol(numbers, positions);
        return molsymm::analyze(mol, tol);
      },
      nb::arg("numbers"),
      nb::arg("positions"),
      nb::arg("primary") = 1e-3);

  m.def(
      "zmatrix_to_xyz",
      [](nb::sequence entries) {
        return positions_to_rows(molsymm::zmatrix_positions(zrows_from_python(entries)));
      },
      nb::arg("entries"));

  m.def(
      "analyze_zmatrix",
      [](nb::sequence entries, double primary) {
        molsymm::Tolerances tol;
        tol.primary = primary;
        molsymm::Molecule mol = molsymm::molecule_from_zmatrix(zrows_from_python(entries));
        return molsymm::analyze(mol, tol);
      },
      nb::arg("entries"),
      nb::arg("primary") = 1e-3);
}
