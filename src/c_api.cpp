#include "molsymm/molsymm.h"

#include "molsymm/molsymm.hpp"

#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

struct molsymm_analysis {
  molsymm::Analysis a;
};

struct molsymm_displacements {
  molsymm::DisplacementSet set;
};

namespace {

thread_local std::string g_error;

void set_error(const std::string& msg) {
  g_error = msg;
}

template <typename F>
int guard(F&& f) {
  try {
    f();
    return 0;
  } catch (const std::exception& ex) {
    set_error(ex.what());
    return -1;
  }
}

} // namespace

extern "C" {

const char* molsymm_last_error(void) {
  return g_error.c_str();
}

int molsymm_analyze(const int* atomic_numbers,
                    const double* xyz,
                    int n,
                    double primary_tol,
                    molsymm_analysis** out) {
  return guard([&] {
    if (!atomic_numbers || !xyz || !out || n <= 0) {
      throw molsymm::Error("invalid arguments");
    }
    std::vector<int> z(atomic_numbers, atomic_numbers + n);
    std::vector<double> coords(xyz, xyz + 3 * n);
    molsymm::Molecule mol(z, coords);
    molsymm::Tolerances tol;
    tol.primary = primary_tol;
    auto* handle = new molsymm_analysis;
    handle->a = molsymm::analyze(mol, tol);
    *out = handle;
  });
}

void molsymm_free(molsymm_analysis* a) {
  delete a;
}

int molsymm_point_group(const molsymm_analysis* a, char* buf, int buflen) {
  return guard([&] {
    if (!a || !buf || buflen <= 0) {
      throw molsymm::Error("invalid arguments");
    }
    std::snprintf(buf, static_cast<std::size_t>(buflen), "%s", a->a.point_group.c_str());
  });
}

int molsymm_symmetry_number(const molsymm_analysis* a, int* out) {
  return guard([&] {
    if (!a || !out) {
      throw molsymm::Error("invalid arguments");
    }
    *out = a->a.symmetry_number;
  });
}

int molsymm_order(const molsymm_analysis* a, int* out) {
  return guard([&] {
    if (!a || !out) {
      throw molsymm::Error("invalid arguments");
    }
    *out = a->a.order;
  });
}

int molsymm_rmsd(const molsymm_analysis* a, double* out) {
  return guard([&] {
    if (!a || !out) {
      throw molsymm::Error("invalid arguments");
    }
    *out = a->a.rmsd;
  });
}

int molsymm_n_atoms(const molsymm_analysis* a, int* out) {
  return guard([&] {
    if (!a || !out) {
      throw molsymm::Error("invalid arguments");
    }
    *out = a->a.original.size();
  });
}

int molsymm_n_operations(const molsymm_analysis* a, int* out) {
  return guard([&] {
    if (!a || !out) {
      throw molsymm::Error("invalid arguments");
    }
    *out = static_cast<int>(a->a.operations.size());
  });
}

int molsymm_operation_matrix(const molsymm_analysis* a, int index, double* mat9) {
  return guard([&] {
    if (!a || !mat9) {
      throw molsymm::Error("invalid arguments");
    }
    if (index < 0 || index >= static_cast<int>(a->a.operations.size())) {
      throw molsymm::Error("operation index out of range");
    }
    const auto& M = a->a.operations[static_cast<std::size_t>(index)].matrix;
    for (int r = 0; r < 3; ++r) {
      for (int c = 0; c < 3; ++c) {
        mat9[3 * r + c] = M(r, c);
      }
    }
  });
}

int molsymm_n_irreps(const molsymm_analysis* a, int* out) {
  return guard([&] {
    if (!a || !out) {
      throw molsymm::Error("invalid arguments");
    }
    *out = static_cast<int>(a->a.character_table.irrep_names.size());
  });
}

int molsymm_irrep_name(const molsymm_analysis* a, int index, char* buf, int buflen) {
  return guard([&] {
    if (!a || !buf || buflen <= 0) {
      throw molsymm::Error("invalid arguments");
    }
    if (index < 0 || index >= static_cast<int>(a->a.character_table.irrep_names.size())) {
      throw molsymm::Error("irrep index out of range");
    }
    std::snprintf(buf,
                  static_cast<std::size_t>(buflen),
                  "%s",
                  a->a.character_table.irrep_names[static_cast<std::size_t>(index)].c_str());
  });
}

int molsymm_n_subgroups(const molsymm_analysis* a, int* out) {
  return guard([&] {
    if (!a || !out) {
      throw molsymm::Error("invalid arguments");
    }
    *out = static_cast<int>(a->a.subgroups.size());
  });
}

int molsymm_subgroup_name(const molsymm_analysis* a, int index, char* buf, int buflen) {
  return guard([&] {
    if (!a || !buf || buflen <= 0) {
      throw molsymm::Error("invalid arguments");
    }
    if (index < 0 || index >= static_cast<int>(a->a.subgroups.size())) {
      throw molsymm::Error("subgroup index out of range");
    }
    std::snprintf(buf,
                  static_cast<std::size_t>(buflen),
                  "%s",
                  a->a.subgroups[static_cast<std::size_t>(index)].c_str());
  });
}

int molsymm_select_subgroup(molsymm_analysis* a, const char* name) {
  return guard([&] {
    if (!a || !name) {
      throw molsymm::Error("invalid arguments");
    }
    a->a.select_subgroup(name);
  });
}

int molsymm_symmetrize(molsymm_analysis* a, double* rmsd, double* xyz_out) {
  return guard([&] {
    if (!a) {
      throw molsymm::Error("invalid arguments");
    }
    const double r = a->a.symmetrize();
    if (rmsd) {
      *rmsd = r;
    }
    if (xyz_out) {
      for (int i = 0; i < a->a.molecule.size(); ++i) {
        const molsymm::Vec3 p =
            a->a.frame_to_input(a->a.molecule.atoms[static_cast<std::size_t>(i)].position);
        xyz_out[3 * i] = p.x();
        xyz_out[3 * i + 1] = p.y();
        xyz_out[3 * i + 2] = p.z();
      }
    }
  });
}

int molsymm_unique_displacements(const molsymm_analysis* a,
                                 double step,
                                 molsymm_displacements** out) {
  return guard([&] {
    if (!a || !out) {
      throw molsymm::Error("invalid arguments");
    }
    auto* handle = new molsymm_displacements;
    handle->set = molsymm::unique_displacements(a->a, step);
    *out = handle;
  });
}

void molsymm_free_displacements(molsymm_displacements* d) {
  delete d;
}

int molsymm_n_displacements(const molsymm_displacements* d, int* out) {
  return guard([&] {
    if (!d || !out) {
      throw molsymm::Error("invalid arguments");
    }
    *out = static_cast<int>(d->set.displacements.size());
  });
}

int molsymm_displacement_atom(const molsymm_displacements* d, int index, int* atom) {
  return guard([&] {
    if (!d || !atom) {
      throw molsymm::Error("invalid arguments");
    }
    if (index < 0 || index >= static_cast<int>(d->set.displacements.size())) {
      throw molsymm::Error("index out of range");
    }
    *atom = d->set.displacements[static_cast<std::size_t>(index)].atom;
  });
}

int molsymm_displacement_direction(const molsymm_displacements* d, int index, double* xyz) {
  return guard([&] {
    if (!d || !xyz) {
      throw molsymm::Error("invalid arguments");
    }
    if (index < 0 || index >= static_cast<int>(d->set.displacements.size())) {
      throw molsymm::Error("index out of range");
    }
    const auto& v = d->set.displacements[static_cast<std::size_t>(index)].direction;
    xyz[0] = v.x();
    xyz[1] = v.y();
    xyz[2] = v.z();
  });
}

int molsymm_displacement_preserved(const molsymm_displacements* d,
                                   int index,
                                   char* buf,
                                   int buflen) {
  return guard([&] {
    if (!d || !buf || buflen <= 0) {
      throw molsymm::Error("invalid arguments");
    }
    if (index < 0 || index >= static_cast<int>(d->set.displacements.size())) {
      throw molsymm::Error("index out of range");
    }
    std::snprintf(buf,
                  static_cast<std::size_t>(buflen),
                  "%s",
                  d->set.displacements[static_cast<std::size_t>(index)].preserved_subgroup.c_str());
  });
}

int molsymm_displacement_geometry(const molsymm_displacements* d,
                                  int index,
                                  int plus,
                                  double* xyz_out) {
  return guard([&] {
    if (!d || !xyz_out) {
      throw molsymm::Error("invalid arguments");
    }
    if (index < 0 || index >= static_cast<int>(d->set.displacements.size())) {
      throw molsymm::Error("index out of range");
    }
    const molsymm::Molecule& m = plus ? d->set.plus[static_cast<std::size_t>(index)]
                                      : d->set.minus[static_cast<std::size_t>(index)];
    for (int i = 0; i < m.size(); ++i) {
      xyz_out[3 * i] = m.atoms[static_cast<std::size_t>(i)].position.x();
      xyz_out[3 * i + 1] = m.atoms[static_cast<std::size_t>(i)].position.y();
      xyz_out[3 * i + 2] = m.atoms[static_cast<std::size_t>(i)].position.z();
    }
  });
}

int molsymm_gradient_from_energies(const molsymm_displacements* d,
                                   const double* e_plus,
                                   const double* e_minus,
                                   double* grad_out) {
  return guard([&] {
    if (!d || !e_plus || !e_minus || !grad_out) {
      throw molsymm::Error("invalid arguments");
    }
    const int nu = static_cast<int>(d->set.displacements.size());
    std::vector<double> ep(e_plus, e_plus + nu);
    std::vector<double> em(e_minus, e_minus + nu);
    Eigen::VectorXd g = d->set.gradient_from_energies(ep, em);
    for (int i = 0; i < g.size(); ++i) {
      grad_out[i] = g(i);
    }
  });
}

int molsymm_hessian_from_gradients(const molsymm_displacements* d,
                                   const double* g_plus,
                                   const double* g_minus,
                                   double* hess_out) {
  return guard([&] {
    if (!d || !g_plus || !g_minus || !hess_out) {
      throw molsymm::Error("invalid arguments");
    }
    const int nu = static_cast<int>(d->set.displacements.size());
    const int n3 = 3 * d->set.analysis.original.size();
    std::vector<Eigen::VectorXd> gp(nu);
    std::vector<Eigen::VectorXd> gm(nu);
    for (int u = 0; u < nu; ++u) {
      gp[static_cast<std::size_t>(u)] = Eigen::VectorXd::Zero(n3);
      gm[static_cast<std::size_t>(u)] = Eigen::VectorXd::Zero(n3);
      for (int i = 0; i < n3; ++i) {
        gp[static_cast<std::size_t>(u)](i) = g_plus[u * n3 + i];
        gm[static_cast<std::size_t>(u)](i) = g_minus[u * n3 + i];
      }
    }
    Eigen::MatrixXd H = d->set.hessian_from_gradients(gp, gm);
    for (int r = 0; r < n3; ++r) {
      for (int c = 0; c < n3; ++c) {
        hess_out[r * n3 + c] = H(r, c);
      }
    }
  });
}

} // extern "C"
