#include "molsymm/molsymm.hpp"

#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <utility>
#include <vector>

namespace nb = nanobind;

namespace {

using NpMatrix = nb::ndarray<nb::numpy, double, nb::ndim<2>, nb::c_contig, nb::device::cpu>;
using NpVector = nb::ndarray<nb::numpy, double, nb::ndim<1>, nb::c_contig, nb::device::cpu>;

constexpr const char* kZEntryHint =
    "z-matrix entry must be (Z), (Z, bond_to, R), (Z, bond_to, R, angle_to, angle), "
    "or (Z, bond_to, R, angle_to, angle, dihedral_to, dihedral)";

molsymm::ZMatrixAtom zatom_from_python(nb::handle entry) {
  if (!nb::isinstance<nb::sequence>(entry) || nb::isinstance<nb::str>(entry)) {
    throw molsymm::Error(kZEntryHint);
  }
  auto seq = nb::cast<nb::sequence>(entry);
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

std::vector<molsymm::ZMatrixAtom> zrows_from_python(const nb::sequence& entries) {
  std::vector<molsymm::ZMatrixAtom> rows;
  rows.reserve(nb::len(entries));
  for (nb::handle entry : entries) {
    rows.push_back(zatom_from_python(entry));
  }
  return rows;
}

// Owned C-contiguous buffers; Eigen map/view returns corrupted columns on some builds.
nb::object cast_move(NpMatrix m) {
  return m.cast(nb::rv_policy::move);
}
nb::object cast_move(NpVector v) {
  return v.cast(nb::rv_policy::move);
}

NpMatrix owned_matrix(std::size_t rows, std::size_t cols) {
  auto* data = new double[rows * cols];
  nb::capsule owner(data, [](void* p) noexcept { delete[] static_cast<double*>(p); });
  return NpMatrix(data, {rows, cols}, owner);
}

NpVector owned_vector(std::size_t n) {
  auto* data = new double[n];
  nb::capsule owner(data, [](void* p) noexcept { delete[] static_cast<double*>(p); });
  return NpVector(data, {n}, owner);
}

NpMatrix positions_to_rows(const std::vector<molsymm::Vec3>& pos) {
  const auto n = pos.size();
  NpMatrix xyz = owned_matrix(n, 3);
  double* data = xyz.data();
  for (std::size_t i = 0; i < n; ++i) {
    data[3 * i + 0] = pos[i].x();
    data[3 * i + 1] = pos[i].y();
    data[3 * i + 2] = pos[i].z();
  }
  return xyz;
}

NpMatrix molecule_to_rows(const molsymm::Molecule& mol) {
  const auto n = static_cast<std::size_t>(mol.size());
  NpMatrix xyz = owned_matrix(n, 3);
  double* data = xyz.data();
  for (std::size_t i = 0; i < n; ++i) {
    const molsymm::Vec3& p = mol.atoms[i].position;
    data[3 * i + 0] = p.x();
    data[3 * i + 1] = p.y();
    data[3 * i + 2] = p.z();
  }
  return xyz;
}

nb::list molecules_to_list(const std::vector<molsymm::Molecule>& mols) {
  nb::list out;
  for (const auto& mol : mols) {
    out.append(cast_move(molecule_to_rows(mol)));
  }
  return out;
}

template <typename Derived>
NpMatrix matrix_to_numpy(const Eigen::MatrixBase<Derived>& m) {
  const auto rows = static_cast<std::size_t>(m.rows());
  const auto cols = static_cast<std::size_t>(m.cols());
  NpMatrix out = owned_matrix(rows, cols);
  double* data = out.data();
  for (std::size_t r = 0; r < rows; ++r) {
    for (std::size_t c = 0; c < cols; ++c) {
      data[r * cols + c] = m(static_cast<Eigen::Index>(r), static_cast<Eigen::Index>(c));
    }
  }
  return out;
}

NpVector vector_to_numpy(const Eigen::VectorXd& v) {
  const auto n = static_cast<std::size_t>(v.size());
  NpVector out = owned_vector(n);
  double* data = out.data();
  for (std::size_t i = 0; i < n; ++i) {
    data[i] = v(static_cast<Eigen::Index>(i));
  }
  return out;
}

NpMatrix forces_to_rows(const Eigen::VectorXd& forces) {
  const auto n = static_cast<int>(forces.size() / 3);
  NpMatrix f = owned_matrix(static_cast<std::size_t>(n), 3);
  double* data = f.data();
  for (int i = 0; i < n; ++i) {
    data[static_cast<std::size_t>(3 * i + 0)] = forces(3 * i);
    data[static_cast<std::size_t>(3 * i + 1)] = forces(3 * i + 1);
    data[static_cast<std::size_t>(3 * i + 2)] = forces(3 * i + 2);
  }
  return f;
}

std::vector<molsymm::Vec3> parse_positions(const Eigen::MatrixXd& positions, int expected) {
  Eigen::MatrixXd rows;
  if (positions.cols() == 3 && positions.rows() == expected) {
    rows = positions;
  } else if (positions.rows() == 3 && positions.cols() == expected) {
    rows = positions.transpose();
  } else {
    throw molsymm::Error("positions must be (n,3) or (3,n)");
  }
  std::vector<molsymm::Vec3> out(static_cast<std::size_t>(expected));
  for (int i = 0; i < expected; ++i) {
    out[static_cast<std::size_t>(i)] = rows.row(i).transpose();
  }
  return out;
}

std::vector<std::vector<std::vector<double>>>
free_axes_prop(const std::vector<molsymm::FragmentSite>& sites) {
  std::vector<std::vector<std::vector<double>>> out;
  out.reserve(sites.size());
  for (const auto& site : sites) {
    std::vector<std::vector<double>> axes;
    axes.reserve(site.free_axes.size());
    for (const auto& ax : site.free_axes) {
      axes.push_back({ax.x(), ax.y(), ax.z()});
    }
    out.push_back(std::move(axes));
  }
  return out;
}

} // namespace

NB_MODULE(molsymm, m) {
  m.doc() = "Molecular point-group symmetry, including Z-matrix input";

  nb::class_<molsymm::Analysis>(m, "Analysis")
      .def_ro("point_group", &molsymm::Analysis::point_group)
      .def_ro("order", &molsymm::Analysis::order)
      .def_ro("symmetry_number", &molsymm::Analysis::symmetry_number)
      .def_ro("rmsd", &molsymm::Analysis::rmsd)
      .def_ro("subgroups", &molsymm::Analysis::subgroups)
      .def_prop_ro("atomic_numbers",
                   [](const molsymm::Analysis& a) {
                     std::vector<int> z;
                     z.reserve(static_cast<std::size_t>(a.molecule.size()));
                     for (const auto& atom : a.molecule.atoms) {
                       z.push_back(atom.atomic_number);
                     }
                     return z;
                   })
      .def_prop_ro(
          "frame_xyz",
          [](const molsymm::Analysis& a) { return cast_move(molecule_to_rows(a.molecule)); })
      .def_prop_ro("operations",
                   [](const molsymm::Analysis& a) {
                     nb::list mats;
                     for (const auto& op : a.operations) {
                       mats.append(cast_move(matrix_to_numpy(op.matrix)));
                     }
                     return mats;
                   })
      .def_prop_ro("irrep_names",
                   [](const molsymm::Analysis& a) { return a.character_table.irrep_names; })
      .def("symmetrize",
           [](molsymm::Analysis& a) {
             const double rms = a.symmetrize();
             std::vector<molsymm::Vec3> pos(static_cast<std::size_t>(a.molecule.size()));
             for (int i = 0; i < a.molecule.size(); ++i) {
               pos[static_cast<std::size_t>(i)] =
                   a.frame_to_input(a.molecule.atoms[static_cast<std::size_t>(i)].position);
             }
             return nb::make_tuple(rms, cast_move(positions_to_rows(pos)));
           })
      .def("select_subgroup", &molsymm::Analysis::select_subgroup)
      .def(
          "displacements",
          [](const molsymm::Analysis& a, double step) {
            return molsymm::unique_displacements(a, step);
          },
          nb::arg("step") = 1e-3)
      .def("fragmentation", [](const molsymm::Analysis& a) { return molsymm::fragmentation(a); })
      .def(
          "subgroup_assembly",
          [](const molsymm::Analysis& a, const std::string& name) {
            return molsymm::subgroup_assembly(a, name);
          },
          nb::arg("name"));

  nb::class_<molsymm::DisplacementSet>(m, "DisplacementSet")
      .def_ro("step", &molsymm::DisplacementSet::step)
      .def_prop_ro("n",
                   [](const molsymm::DisplacementSet& d) {
                     return static_cast<int>(d.displacements.size());
                   })
      .def_prop_ro("plus",
                   [](const molsymm::DisplacementSet& d) { return molecules_to_list(d.plus); })
      .def_prop_ro("minus",
                   [](const molsymm::DisplacementSet& d) { return molecules_to_list(d.minus); })
      .def("gradient_from_energies",
           [](const molsymm::DisplacementSet& d,
              const std::vector<double>& ep,
              const std::vector<double>& em) {
             return cast_move(vector_to_numpy(d.gradient_from_energies(ep, em)));
           })
      .def("hessian_from_gradients",
           [](const molsymm::DisplacementSet& d,
              const std::vector<Eigen::VectorXd>& gp,
              const std::vector<Eigen::VectorXd>& gm) {
             return cast_move(matrix_to_numpy(d.hessian_from_gradients(gp, gm)));
           });

  nb::class_<molsymm::Fragmentation>(m, "Fragmentation")
      .def_ro("representatives", &molsymm::Fragmentation::representatives)
      .def_prop_ro("copies",
                   [](const molsymm::Fragmentation& f) {
                     std::vector<std::vector<int>> out;
                     out.reserve(f.copies.size());
                     for (const auto& c : f.copies) {
                       out.push_back(c.atoms);
                     }
                     return out;
                   })
      .def_prop_ro("mask",
                   [](const molsymm::Fragmentation& f) {
                     return std::vector<int>(f.mask.begin(), f.mask.end());
                   })
      .def_prop_ro("free_axes",
                   [](const molsymm::Fragmentation& f) { return free_axes_prop(f.sites); })
      .def(
          "reconstruct",
          [](const molsymm::Fragmentation& f, const Eigen::MatrixXd& positions) {
            return cast_move(molecule_to_rows(molsymm::reconstruct(
                f, parse_positions(positions, static_cast<int>(f.representatives.size())))));
          },
          nb::arg("positions"));

  nb::class_<molsymm::SubgroupAssembly>(m, "SubgroupAssembly")
      .def_ro("subgroup", &molsymm::SubgroupAssembly::subgroup)
      .def_prop_ro("pieces",
                   [](const molsymm::SubgroupAssembly& a) {
                     std::vector<std::vector<int>> out;
                     out.reserve(a.pieces.size());
                     for (const auto& p : a.pieces) {
                       out.push_back(p.atoms);
                     }
                     return out;
                   })
      .def_prop_ro("free_axes",
                   [](const molsymm::SubgroupAssembly& a) { return free_axes_prop(a.sites); })
      .def(
          "assemble",
          [](const molsymm::SubgroupAssembly& a, const nb::sequence& pieces) {
            if (static_cast<int>(nb::len(pieces)) != static_cast<int>(a.pieces.size())) {
              throw molsymm::Error("piece count mismatch");
            }
            std::vector<std::vector<molsymm::Vec3>> piece_pos(a.pieces.size());
            for (std::size_t p = 0; p < a.pieces.size(); ++p) {
              piece_pos[p] = parse_positions(nb::cast<Eigen::MatrixXd>(pieces[p]),
                                             static_cast<int>(a.pieces[p].atoms.size()));
            }
            return cast_move(molecule_to_rows(molsymm::assemble(a, piece_pos)));
          },
          nb::arg("pieces"));

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

  nb::class_<molsymm::Symmetrized>(m, "Symmetrized")
      .def_ro("rmsd", &molsymm::Symmetrized::rmsd)
      .def_prop_ro(
          "positions",
          [](const molsymm::Symmetrized& s) { return cast_move(positions_to_rows(s.positions)); })
      .def_prop_ro("forces",
                   [](const molsymm::Symmetrized& s) -> nb::object {
                     if (s.forces.size() == 0) {
                       return nb::none();
                     }
                     return cast_move(forces_to_rows(s.forces));
                   })
      .def_prop_ro("hessian", [](const molsymm::Symmetrized& s) -> nb::object {
        if (s.hessian.rows() == 0 || s.hessian.cols() == 0) {
          return nb::none();
        }
        return cast_move(matrix_to_numpy(s.hessian));
      });

  m.def(
      "symmetrize",
      [](const std::vector<int>& numbers,
         const Eigen::MatrixXd& positions,
         const std::string& point_group,
         const nb::object& forces,
         const nb::object& hessian,
         double primary) {
        if (positions.cols() != 3 && positions.rows() != 3) {
          throw molsymm::Error("positions must be (n,3) or (3,n)");
        }
        const auto n = static_cast<int>(numbers.size());
        molsymm::Molecule mol;
        if (positions.cols() == 3) {
          if (positions.rows() != n) {
            throw molsymm::Error("numbers/positions size mismatch");
          }
          mol = molsymm::Molecule(numbers, positions.transpose());
        } else {
          if (positions.cols() != n) {
            throw molsymm::Error("numbers/positions size mismatch");
          }
          mol = molsymm::Molecule(numbers, positions);
        }

        Eigen::VectorXd forces_vec;
        const Eigen::VectorXd* forces_ptr = nullptr;
        if (!forces.is_none()) {
          const bool flat = nb::hasattr(forces, "ndim") && nb::cast<int>(forces.attr("ndim")) == 1;
          if (flat) {
            forces_vec = nb::cast<Eigen::VectorXd>(forces);
            if (forces_vec.size() != 3 * n) {
              throw molsymm::Error("forces must be (n,3), (3,n), or flat (3n,)");
            }
          } else {
            auto fmat = nb::cast<Eigen::MatrixXd>(forces);
            forces_vec.resize(3 * n);
            if (fmat.cols() == 3 && fmat.rows() == n) {
              for (int i = 0; i < n; ++i) {
                forces_vec.segment<3>(3 * i) = fmat.row(i).transpose();
              }
            } else if (fmat.rows() == 3 && fmat.cols() == n) {
              for (int i = 0; i < n; ++i) {
                forces_vec.segment<3>(3 * i) = fmat.col(i);
              }
            } else if (fmat.size() == 3 * n) {
              for (int k = 0; k < 3 * n; ++k) {
                forces_vec(k) = fmat.data()[k];
              }
            } else {
              throw molsymm::Error("forces must be (n,3), (3,n), or flat (3n,)");
            }
          }
          forces_ptr = &forces_vec;
        }

        Eigen::MatrixXd hess_mat;
        const Eigen::MatrixXd* hess_ptr = nullptr;
        if (!hessian.is_none()) {
          hess_mat = nb::cast<Eigen::MatrixXd>(hessian);
          if (hess_mat.rows() != 3 * n || hess_mat.cols() != 3 * n) {
            throw molsymm::Error("hessian must be (3n,3n)");
          }
          hess_ptr = &hess_mat;
        }

        molsymm::Tolerances tol;
        tol.primary = primary;
        return molsymm::symmetrize(mol, point_group, forces_ptr, hess_ptr, tol);
      },
      nb::arg("numbers"),
      nb::arg("positions"),
      nb::arg("point_group"),
      nb::arg("forces") = nb::none(),
      nb::arg("hessian") = nb::none(),
      nb::arg("primary") = 1e-3);

  m.def(
      "zmatrix_to_xyz",
      [](const nb::sequence& entries) {
        return cast_move(positions_to_rows(molsymm::zmatrix_positions(zrows_from_python(entries))));
      },
      nb::arg("entries"));

  m.def(
      "analyze_zmatrix",
      [](const nb::sequence& entries, double primary) {
        molsymm::Tolerances tol;
        tol.primary = primary;
        molsymm::Molecule mol = molsymm::molecule_from_zmatrix(zrows_from_python(entries));
        return molsymm::analyze(mol, tol);
      },
      nb::arg("entries"),
      nb::arg("primary") = 1e-3);
}
