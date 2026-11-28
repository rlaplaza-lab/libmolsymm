#include "molsymm/molsymm.h"

#include "molsymm/molsymm.hpp"

#include <Eigen/Dense>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <span>
#include <string>
#include <utility>
#include <vector>

static_assert(sizeof(molsymm_zatom) == sizeof(molsymm::ZMatrixAtom),
              "C and C++ Z-matrix rows must have the same size");
static_assert(offsetof(molsymm_zatom, atomic_number) ==
                  offsetof(molsymm::ZMatrixAtom, atomic_number),
              "atomic_number offset");
static_assert(offsetof(molsymm_zatom, bond_to) == offsetof(molsymm::ZMatrixAtom, bond_to),
              "bond_to offset");
static_assert(offsetof(molsymm_zatom, angle_to) == offsetof(molsymm::ZMatrixAtom, angle_to),
              "angle_to offset");
static_assert(offsetof(molsymm_zatom, dihedral_to) == offsetof(molsymm::ZMatrixAtom, dihedral_to),
              "dihedral_to offset");
static_assert(offsetof(molsymm_zatom, distance) == offsetof(molsymm::ZMatrixAtom, distance),
              "distance offset");
static_assert(offsetof(molsymm_zatom, angle) == offsetof(molsymm::ZMatrixAtom, angle),
              "angle offset");
static_assert(offsetof(molsymm_zatom, dihedral) == offsetof(molsymm::ZMatrixAtom, dihedral),
              "dihedral offset");

struct molsymm_analysis {
  molsymm::Analysis a;
};

struct molsymm_displacements {
  molsymm::DisplacementSet set;
};

struct molsymm_fragmentation {
  molsymm::Fragmentation frag;
};

struct molsymm_assembly {
  molsymm::SubgroupAssembly assembly;
};

struct molsymm_wyckoff_table {
  std::vector<molsymm::WyckoffSite> sites;
};

struct molsymm_wyckoff_plans {
  std::vector<molsymm::WyckoffPlan> plans;
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

std::vector<molsymm::ZMatrixAtom> zrows_from_c(const molsymm_zatom* atoms, int n) {
  std::vector<molsymm::ZMatrixAtom> rows(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    molsymm::ZMatrixAtom& row = rows[static_cast<std::size_t>(i)];
    row.atomic_number = atoms[i].atomic_number;
    row.bond_to = atoms[i].bond_to;
    row.angle_to = atoms[i].angle_to;
    row.dihedral_to = atoms[i].dihedral_to;
    row.distance = atoms[i].distance;
    row.angle = atoms[i].angle;
    row.dihedral = atoms[i].dihedral;
  }
  return rows;
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
    molsymm::Molecule mol(std::span<const int>{atomic_numbers, static_cast<std::size_t>(n)},
                          std::span<const double>{xyz, static_cast<std::size_t>(3 * n)});
    molsymm::Tolerances tol;
    tol.primary = primary_tol;
    auto* handle = new molsymm_analysis;
    handle->a = molsymm::analyze(mol, tol);
    *out = handle;
  });
}

int molsymm_zmatrix_to_xyz(const molsymm_zatom* atoms, int n, double* xyz_out) {
  return guard([&] {
    if (!atoms || !xyz_out || n <= 0) {
      throw molsymm::Error("invalid arguments");
    }
    const std::vector<molsymm::Vec3> pos = molsymm::zmatrix_positions(zrows_from_c(atoms, n));
    for (int i = 0; i < n; ++i) {
      xyz_out[3 * i] = pos[static_cast<std::size_t>(i)].x();
      xyz_out[3 * i + 1] = pos[static_cast<std::size_t>(i)].y();
      xyz_out[3 * i + 2] = pos[static_cast<std::size_t>(i)].z();
    }
  });
}

int molsymm_analyze_zmatrix(const molsymm_zatom* atoms,
                            int n,
                            double primary_tol,
                            molsymm_analysis** out) {
  return guard([&] {
    if (!atoms || !out || n <= 0) {
      throw molsymm::Error("invalid arguments");
    }
    molsymm::Molecule mol = molsymm::molecule_from_zmatrix(zrows_from_c(atoms, n));
    molsymm::Tolerances tol;
    tol.primary = primary_tol;
    molsymm::Analysis analysis = molsymm::analyze(mol, tol);
    auto* handle = new molsymm_analysis;
    handle->a = std::move(analysis);
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

int molsymm_symmetrize_data(const int* atomic_numbers,
                            const double* xyz,
                            int n,
                            const char* point_group,
                            double primary_tol,
                            const double* forces,
                            const double* hessian,
                            double* rmsd,
                            double* xyz_out,
                            double* forces_out,
                            double* hessian_out) {
  return guard([&] {
    if (!atomic_numbers || !xyz || n <= 0 || !point_group || !xyz_out) {
      throw molsymm::Error("invalid arguments");
    }
    molsymm::Molecule mol(std::span<const int>{atomic_numbers, static_cast<std::size_t>(n)},
                          std::span<const double>{xyz, static_cast<std::size_t>(3 * n)});
    molsymm::Tolerances tol;
    tol.primary = primary_tol;

    Eigen::VectorXd forces_vec;
    const Eigen::VectorXd* forces_ptr = nullptr;
    if (forces) {
      forces_vec = Eigen::VectorXd::Map(forces, 3 * n);
      forces_ptr = &forces_vec;
    }
    Eigen::MatrixXd hess_mat;
    const Eigen::MatrixXd* hess_ptr = nullptr;
    if (hessian) {
      hess_mat =
          Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
              hessian, 3 * n, 3 * n);
      hess_ptr = &hess_mat;
    }

    molsymm::Symmetrized out = molsymm::symmetrize(mol, point_group, forces_ptr, hess_ptr, tol);
    if (rmsd) {
      *rmsd = out.rmsd;
    }
    for (int i = 0; i < n; ++i) {
      const molsymm::Vec3& p = out.positions[static_cast<std::size_t>(i)];
      xyz_out[3 * i] = p.x();
      xyz_out[3 * i + 1] = p.y();
      xyz_out[3 * i + 2] = p.z();
    }
    if (forces && forces_out) {
      Eigen::Map<Eigen::VectorXd>(forces_out, 3 * n) = out.forces;
    }
    if (hessian && hessian_out) {
      Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
          hessian_out, 3 * n, 3 * n) = out.hessian;
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
    Eigen::VectorXd g = d->set.gradient_from_energies(
        std::span<const double>{e_plus, static_cast<std::size_t>(nu)},
        std::span<const double>{e_minus, static_cast<std::size_t>(nu)});
    Eigen::Map<Eigen::VectorXd>(grad_out, g.size()) = g;
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
    std::vector<Eigen::VectorXd> gp(static_cast<std::size_t>(nu));
    std::vector<Eigen::VectorXd> gm(static_cast<std::size_t>(nu));
    for (int u = 0; u < nu; ++u) {
      gp[static_cast<std::size_t>(u)] = Eigen::VectorXd::Map(g_plus + u * n3, n3);
      gm[static_cast<std::size_t>(u)] = Eigen::VectorXd::Map(g_minus + u * n3, n3);
    }
    Eigen::MatrixXd H = d->set.hessian_from_gradients(gp, gm);
    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
        hess_out, n3, n3) = H;
  });
}

int molsymm_make_fragmentation(const molsymm_analysis* a, molsymm_fragmentation** out) {
  return guard([&] {
    if (!a || !out) {
      throw molsymm::Error("invalid arguments");
    }
    molsymm::Fragmentation frag = molsymm::fragmentation(a->a);
    auto* handle = new molsymm_fragmentation;
    handle->frag = std::move(frag);
    *out = handle;
  });
}

void molsymm_free_fragmentation(molsymm_fragmentation* f) {
  delete f;
}

int molsymm_n_representatives(const molsymm_fragmentation* f, int* out) {
  return guard([&] {
    if (!f || !out) {
      throw molsymm::Error("invalid arguments");
    }
    *out = static_cast<int>(f->frag.representatives.size());
  });
}

int molsymm_representative_atom(const molsymm_fragmentation* f, int index, int* atom) {
  return guard([&] {
    if (!f || !atom) {
      throw molsymm::Error("invalid arguments");
    }
    if (index < 0 || index >= static_cast<int>(f->frag.representatives.size())) {
      throw molsymm::Error("index out of range");
    }
    *atom = f->frag.representatives[static_cast<std::size_t>(index)];
  });
}

int molsymm_n_copies(const molsymm_fragmentation* f, int* out) {
  return guard([&] {
    if (!f || !out) {
      throw molsymm::Error("invalid arguments");
    }
    *out = static_cast<int>(f->frag.copies.size());
  });
}

int molsymm_copy_n_atoms(const molsymm_fragmentation* f, int copy, int* out) {
  return guard([&] {
    if (!f || !out) {
      throw molsymm::Error("invalid arguments");
    }
    if (copy < 0 || copy >= static_cast<int>(f->frag.copies.size())) {
      throw molsymm::Error("copy index out of range");
    }
    *out = static_cast<int>(f->frag.copies[static_cast<std::size_t>(copy)].atoms.size());
  });
}

int molsymm_copy_atoms(const molsymm_fragmentation* f, int copy, int* atoms_out) {
  return guard([&] {
    if (!f || !atoms_out) {
      throw molsymm::Error("invalid arguments");
    }
    if (copy < 0 || copy >= static_cast<int>(f->frag.copies.size())) {
      throw molsymm::Error("copy index out of range");
    }
    const auto& atoms = f->frag.copies[static_cast<std::size_t>(copy)].atoms;
    for (std::size_t i = 0; i < atoms.size(); ++i) {
      atoms_out[i] = atoms[i];
    }
  });
}

int molsymm_fragment_mask(const molsymm_fragmentation* f, unsigned char* mask_out) {
  return guard([&] {
    if (!f || !mask_out) {
      throw molsymm::Error("invalid arguments");
    }
    for (std::size_t i = 0; i < f->frag.mask.size(); ++i) {
      mask_out[i] = f->frag.mask[i];
    }
  });
}

int molsymm_n_free_axes(const molsymm_fragmentation* f, int site, int* out) {
  return guard([&] {
    if (!f || !out) {
      throw molsymm::Error("invalid arguments");
    }
    if (site < 0 || site >= static_cast<int>(f->frag.sites.size())) {
      throw molsymm::Error("site index out of range");
    }
    *out = static_cast<int>(f->frag.sites[static_cast<std::size_t>(site)].free_axes.size());
  });
}

int molsymm_free_axis(const molsymm_fragmentation* f, int site, int axis, double* xyz) {
  return guard([&] {
    if (!f || !xyz) {
      throw molsymm::Error("invalid arguments");
    }
    if (site < 0 || site >= static_cast<int>(f->frag.sites.size())) {
      throw molsymm::Error("site index out of range");
    }
    const auto& axes = f->frag.sites[static_cast<std::size_t>(site)].free_axes;
    if (axis < 0 || axis >= static_cast<int>(axes.size())) {
      throw molsymm::Error("axis index out of range");
    }
    const auto& v = axes[static_cast<std::size_t>(axis)];
    xyz[0] = v.x();
    xyz[1] = v.y();
    xyz[2] = v.z();
  });
}

int molsymm_reconstruct(const molsymm_fragmentation* f, const double* rep_xyz, double* xyz_out) {
  return guard([&] {
    if (!f || !rep_xyz || !xyz_out) {
      throw molsymm::Error("invalid arguments");
    }
    const int nrep = static_cast<int>(f->frag.representatives.size());
    std::vector<molsymm::Vec3> reps(static_cast<std::size_t>(nrep));
    for (int i = 0; i < nrep; ++i) {
      reps[static_cast<std::size_t>(i)] =
          molsymm::Vec3(rep_xyz[3 * i], rep_xyz[3 * i + 1], rep_xyz[3 * i + 2]);
    }
    molsymm::Molecule mol = molsymm::reconstruct(f->frag, reps);
    for (int i = 0; i < mol.size(); ++i) {
      const auto& p = mol.atoms[static_cast<std::size_t>(i)].position;
      xyz_out[3 * i] = p.x();
      xyz_out[3 * i + 1] = p.y();
      xyz_out[3 * i + 2] = p.z();
    }
  });
}

int molsymm_subgroup_assembly(const molsymm_analysis* a,
                              const char* subgroup,
                              molsymm_assembly** out) {
  return guard([&] {
    if (!a || !subgroup || !out) {
      throw molsymm::Error("invalid arguments");
    }
    molsymm::SubgroupAssembly assembly = molsymm::subgroup_assembly(a->a, subgroup);
    auto* handle = new molsymm_assembly;
    handle->assembly = std::move(assembly);
    *out = handle;
  });
}

void molsymm_free_assembly(molsymm_assembly* asm_) {
  delete asm_;
}

int molsymm_assembly_n_pieces(const molsymm_assembly* asm_, int* out) {
  return guard([&] {
    if (!asm_ || !out) {
      throw molsymm::Error("invalid arguments");
    }
    *out = static_cast<int>(asm_->assembly.pieces.size());
  });
}

int molsymm_assembly_piece_n_atoms(const molsymm_assembly* asm_, int piece, int* out) {
  return guard([&] {
    if (!asm_ || !out) {
      throw molsymm::Error("invalid arguments");
    }
    if (piece < 0 || piece >= static_cast<int>(asm_->assembly.pieces.size())) {
      throw molsymm::Error("piece index out of range");
    }
    *out = static_cast<int>(asm_->assembly.pieces[static_cast<std::size_t>(piece)].atoms.size());
  });
}

int molsymm_assembly_piece_atoms(const molsymm_assembly* asm_, int piece, int* atoms_out) {
  return guard([&] {
    if (!asm_ || !atoms_out) {
      throw molsymm::Error("invalid arguments");
    }
    if (piece < 0 || piece >= static_cast<int>(asm_->assembly.pieces.size())) {
      throw molsymm::Error("piece index out of range");
    }
    const auto& atoms = asm_->assembly.pieces[static_cast<std::size_t>(piece)].atoms;
    for (std::size_t i = 0; i < atoms.size(); ++i) {
      atoms_out[i] = atoms[i];
    }
  });
}

int molsymm_assembly_n_free_axes(const molsymm_assembly* asm_, int site, int* out) {
  return guard([&] {
    if (!asm_ || !out) {
      throw molsymm::Error("invalid arguments");
    }
    if (site < 0 || site >= static_cast<int>(asm_->assembly.sites.size())) {
      throw molsymm::Error("site index out of range");
    }
    *out = static_cast<int>(asm_->assembly.sites[static_cast<std::size_t>(site)].free_axes.size());
  });
}

int molsymm_assembly_free_axis(const molsymm_assembly* asm_, int site, int axis, double* xyz) {
  return guard([&] {
    if (!asm_ || !xyz) {
      throw molsymm::Error("invalid arguments");
    }
    if (site < 0 || site >= static_cast<int>(asm_->assembly.sites.size())) {
      throw molsymm::Error("site index out of range");
    }
    const auto& axes = asm_->assembly.sites[static_cast<std::size_t>(site)].free_axes;
    if (axis < 0 || axis >= static_cast<int>(axes.size())) {
      throw molsymm::Error("axis index out of range");
    }
    const auto& v = axes[static_cast<std::size_t>(axis)];
    xyz[0] = v.x();
    xyz[1] = v.y();
    xyz[2] = v.z();
  });
}

int molsymm_assemble(const molsymm_assembly* asm_, const double* piece_xyz, double* xyz_out) {
  return guard([&] {
    if (!asm_ || !piece_xyz || !xyz_out) {
      throw molsymm::Error("invalid arguments");
    }
    std::vector<std::vector<molsymm::Vec3>> pieces(asm_->assembly.pieces.size());
    int offset = 0;
    for (std::size_t p = 0; p < asm_->assembly.pieces.size(); ++p) {
      const int na = static_cast<int>(asm_->assembly.pieces[p].atoms.size());
      pieces[p].resize(static_cast<std::size_t>(na));
      for (int i = 0; i < na; ++i) {
        pieces[p][static_cast<std::size_t>(i)] = molsymm::Vec3(piece_xyz[3 * (offset + i)],
                                                               piece_xyz[3 * (offset + i) + 1],
                                                               piece_xyz[3 * (offset + i) + 2]);
      }
      offset += na;
    }
    molsymm::Molecule mol = molsymm::assemble(asm_->assembly, pieces);
    for (int i = 0; i < mol.size(); ++i) {
      const auto& p = mol.atoms[static_cast<std::size_t>(i)].position;
      xyz_out[3 * i] = p.x();
      xyz_out[3 * i + 1] = p.y();
      xyz_out[3 * i + 2] = p.z();
    }
  });
}

int molsymm_wyckoff_positions(const char* point_group, molsymm_wyckoff_table** out) {
  return guard([&] {
    if (!point_group || !out) {
      throw molsymm::Error("invalid arguments");
    }
    auto* handle = new molsymm_wyckoff_table;
    handle->sites = molsymm::wyckoff_positions(point_group);
    *out = handle;
  });
}

void molsymm_free_wyckoff_table(molsymm_wyckoff_table* t) {
  delete t;
}

int molsymm_wyckoff_n_sites(const molsymm_wyckoff_table* t, int* out) {
  return guard([&] {
    if (!t || !out) {
      throw molsymm::Error("invalid arguments");
    }
    *out = static_cast<int>(t->sites.size());
  });
}

int molsymm_wyckoff_site_label(const molsymm_wyckoff_table* t, int index, char* buf, int buflen) {
  return guard([&] {
    if (!t || !buf || buflen <= 0) {
      throw molsymm::Error("invalid arguments");
    }
    if (index < 0 || index >= static_cast<int>(t->sites.size())) {
      throw molsymm::Error("index out of range");
    }
    const std::string& s = t->sites[static_cast<std::size_t>(index)].label;
    if (static_cast<int>(s.size()) + 1 > buflen) {
      throw molsymm::Error("buffer too small");
    }
    std::snprintf(buf, static_cast<std::size_t>(buflen), "%s", s.c_str());
  });
}

int molsymm_wyckoff_site_symmetry(const molsymm_wyckoff_table* t,
                                  int index,
                                  char* buf,
                                  int buflen) {
  return guard([&] {
    if (!t || !buf || buflen <= 0) {
      throw molsymm::Error("invalid arguments");
    }
    if (index < 0 || index >= static_cast<int>(t->sites.size())) {
      throw molsymm::Error("index out of range");
    }
    const std::string& s = t->sites[static_cast<std::size_t>(index)].site_symmetry;
    if (static_cast<int>(s.size()) + 1 > buflen) {
      throw molsymm::Error("buffer too small");
    }
    std::snprintf(buf, static_cast<std::size_t>(buflen), "%s", s.c_str());
  });
}

int molsymm_wyckoff_site_multiplicity(const molsymm_wyckoff_table* t, int index, int* out) {
  return guard([&] {
    if (!t || !out) {
      throw molsymm::Error("invalid arguments");
    }
    if (index < 0 || index >= static_cast<int>(t->sites.size())) {
      throw molsymm::Error("index out of range");
    }
    *out = t->sites[static_cast<std::size_t>(index)].multiplicity;
  });
}

int molsymm_wyckoff_site_dof(const molsymm_wyckoff_table* t, int index, int* out) {
  return guard([&] {
    if (!t || !out) {
      throw molsymm::Error("invalid arguments");
    }
    if (index < 0 || index >= static_cast<int>(t->sites.size())) {
      throw molsymm::Error("index out of range");
    }
    *out = t->sites[static_cast<std::size_t>(index)].degrees_of_freedom;
  });
}

int molsymm_wyckoff_plans_composition(const char* point_group,
                                      const int* atomic_numbers,
                                      int n,
                                      molsymm_wyckoff_plans** out) {
  return guard([&] {
    if (!point_group || !out || n < 0 || (n > 0 && !atomic_numbers)) {
      throw molsymm::Error("invalid arguments");
    }
    std::vector<int> zs(atomic_numbers, atomic_numbers + n);
    auto* handle = new molsymm_wyckoff_plans;
    handle->plans = molsymm::wyckoff_plans(point_group, zs);
    *out = handle;
  });
}

int molsymm_wyckoff_plans_nsites(const char* point_group,
                                 int n_sites,
                                 molsymm_wyckoff_plans** out) {
  return guard([&] {
    if (!point_group || !out) {
      throw molsymm::Error("invalid arguments");
    }
    auto* handle = new molsymm_wyckoff_plans;
    handle->plans = molsymm::wyckoff_plans(point_group, n_sites);
    *out = handle;
  });
}

void molsymm_free_wyckoff_plans(molsymm_wyckoff_plans* p) {
  delete p;
}

int molsymm_wyckoff_n_plans(const molsymm_wyckoff_plans* p, int* out) {
  return guard([&] {
    if (!p || !out) {
      throw molsymm::Error("invalid arguments");
    }
    *out = static_cast<int>(p->plans.size());
  });
}

int molsymm_wyckoff_plan_n_orbits(const molsymm_wyckoff_plans* p, int plan, int* out) {
  return guard([&] {
    if (!p || !out) {
      throw molsymm::Error("invalid arguments");
    }
    if (plan < 0 || plan >= static_cast<int>(p->plans.size())) {
      throw molsymm::Error("plan index out of range");
    }
    *out = static_cast<int>(p->plans[static_cast<std::size_t>(plan)].orbits.size());
  });
}

int molsymm_wyckoff_plan_orbit_label(
    const molsymm_wyckoff_plans* p, int plan, int orbit, char* buf, int buflen) {
  return guard([&] {
    if (!p || !buf || buflen <= 0) {
      throw molsymm::Error("invalid arguments");
    }
    if (plan < 0 || plan >= static_cast<int>(p->plans.size())) {
      throw molsymm::Error("plan index out of range");
    }
    const auto& orbits = p->plans[static_cast<std::size_t>(plan)].orbits;
    if (orbit < 0 || orbit >= static_cast<int>(orbits.size())) {
      throw molsymm::Error("orbit index out of range");
    }
    const std::string& s = orbits[static_cast<std::size_t>(orbit)].label;
    if (static_cast<int>(s.size()) + 1 > buflen) {
      throw molsymm::Error("buffer too small");
    }
    std::snprintf(buf, static_cast<std::size_t>(buflen), "%s", s.c_str());
  });
}

int molsymm_wyckoff_plan_orbit_z(const molsymm_wyckoff_plans* p, int plan, int orbit, int* out) {
  return guard([&] {
    if (!p || !out) {
      throw molsymm::Error("invalid arguments");
    }
    if (plan < 0 || plan >= static_cast<int>(p->plans.size())) {
      throw molsymm::Error("plan index out of range");
    }
    const auto& orbits = p->plans[static_cast<std::size_t>(plan)].orbits;
    if (orbit < 0 || orbit >= static_cast<int>(orbits.size())) {
      throw molsymm::Error("orbit index out of range");
    }
    *out = orbits[static_cast<std::size_t>(orbit)].atomic_number;
  });
}

int molsymm_assemble_from_au(const char* point_group,
                             const int* representative_numbers,
                             const double* au_xyz,
                             int n_orbits,
                             const char** site_labels,
                             molsymm_analysis** out) {
  return guard([&] {
    if (!point_group || !out || n_orbits < 0 ||
        (n_orbits > 0 && (!representative_numbers || !au_xyz || !site_labels))) {
      throw molsymm::Error("invalid arguments");
    }
    std::vector<int> zs(representative_numbers, representative_numbers + n_orbits);
    std::vector<molsymm::Vec3> pos(static_cast<std::size_t>(n_orbits));
    std::vector<std::string> labels(static_cast<std::size_t>(n_orbits));
    for (int i = 0; i < n_orbits; ++i) {
      pos[static_cast<std::size_t>(i)] =
          molsymm::Vec3(au_xyz[3 * i], au_xyz[3 * i + 1], au_xyz[3 * i + 2]);
      if (!site_labels[i]) {
        throw molsymm::Error("null site label");
      }
      labels[static_cast<std::size_t>(i)] = site_labels[i];
    }
    auto* handle = new molsymm_analysis;
    handle->a = molsymm::assemble_from_au(point_group, zs, pos, labels);
    *out = handle;
  });
}

} // extern "C"
