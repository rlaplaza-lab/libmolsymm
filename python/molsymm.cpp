#include "molsymm/molsymm.hpp"

#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <utility>
#include <vector>

namespace nb = nanobind;

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
        Eigen::MatrixXd xyz;
        if (positions.cols() == 3) {
          xyz = positions.transpose(); // 3 x n
        } else {
          xyz = positions;
        }
        if (static_cast<int>(numbers.size()) != xyz.cols()) {
          throw molsymm::Error("numbers/positions size mismatch");
        }
        molsymm::Molecule mol(numbers, xyz);
        molsymm::Tolerances tol;
        tol.primary = primary;
        return molsymm::analyze(mol, tol);
      },
      nb::arg("numbers"),
      nb::arg("positions"),
      nb::arg("primary") = 1e-3);
}
