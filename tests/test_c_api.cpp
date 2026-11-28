#include "molsymm/molsymm.h"

#include <cmath>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

int g_failed = 0;
int g_passed = 0;

void expect(bool cond, const std::string& msg) {
  if (!cond) {
    std::cerr << "FAIL: " << msg << "\n";
    ++g_failed;
  } else {
    ++g_passed;
  }
}

constexpr int kWaterZ[3] = {8, 1, 1};
constexpr double kWaterXyz[9] = {0.0, 0.0, 0.1173, 0.0, 0.7572, -0.4692, 0.0, -0.7572, -0.4692};

void test_c_api_nulls() {
  expect(molsymm_analyze(nullptr, kWaterXyz, 3, 1e-3, nullptr) == -1, "analyze null");
  expect(std::strlen(molsymm_last_error()) > 0, "last_error set");
  expect(molsymm_analyze(kWaterZ, kWaterXyz, 0, 1e-3, nullptr) == -1, "analyze n=0");
  expect(molsymm_point_group(nullptr, nullptr, 0) == -1, "point_group null");
  expect(molsymm_symmetry_number(nullptr, nullptr) == -1, "symmetry_number null");
  expect(molsymm_order(nullptr, nullptr) == -1, "order null");
  expect(molsymm_rmsd(nullptr, nullptr) == -1, "rmsd null");
  expect(molsymm_n_atoms(nullptr, nullptr) == -1, "n_atoms null");
  expect(molsymm_n_operations(nullptr, nullptr) == -1, "n_operations null");
  expect(molsymm_operation_matrix(nullptr, 0, nullptr) == -1, "operation_matrix null");
  expect(molsymm_n_irreps(nullptr, nullptr) == -1, "n_irreps null");
  expect(molsymm_irrep_name(nullptr, 0, nullptr, 0) == -1, "irrep_name null");
  expect(molsymm_n_subgroups(nullptr, nullptr) == -1, "n_subgroups null");
  expect(molsymm_subgroup_name(nullptr, 0, nullptr, 0) == -1, "subgroup_name null");
  expect(molsymm_select_subgroup(nullptr, nullptr) == -1, "select_subgroup null");
  expect(molsymm_symmetrize(nullptr, nullptr, nullptr) == -1, "symmetrize null");
  expect(molsymm_symmetrize_data(nullptr,
                                 nullptr,
                                 0,
                                 nullptr,
                                 1e-3,
                                 nullptr,
                                 nullptr,
                                 nullptr,
                                 nullptr,
                                 nullptr,
                                 nullptr) == -1,
         "symmetrize_data null");
  expect(molsymm_unique_displacements(nullptr, 1e-3, nullptr) == -1, "unique_displacements null");
  expect(molsymm_n_displacements(nullptr, nullptr) == -1, "n_displacements null");
  expect(molsymm_displacement_atom(nullptr, 0, nullptr) == -1, "displacement_atom null");
  expect(molsymm_displacement_direction(nullptr, 0, nullptr) == -1, "displacement_direction null");
  expect(molsymm_displacement_preserved(nullptr, 0, nullptr, 0) == -1,
         "displacement_preserved null");
  expect(molsymm_displacement_geometry(nullptr, 0, 1, nullptr) == -1, "displacement_geometry null");
  expect(molsymm_gradient_from_energies(nullptr, nullptr, nullptr, nullptr) == -1,
         "gradient_from_energies null");
  expect(molsymm_hessian_from_gradients(nullptr, nullptr, nullptr, nullptr) == -1,
         "hessian_from_gradients null");
  expect(molsymm_make_fragmentation(nullptr, nullptr) == -1, "make_fragmentation null");
  expect(molsymm_n_representatives(nullptr, nullptr) == -1, "n_representatives null");
  expect(molsymm_representative_atom(nullptr, 0, nullptr) == -1, "representative_atom null");
  expect(molsymm_n_copies(nullptr, nullptr) == -1, "n_copies null");
  expect(molsymm_copy_n_atoms(nullptr, 0, nullptr) == -1, "copy_n_atoms null");
  expect(molsymm_copy_atoms(nullptr, 0, nullptr) == -1, "copy_atoms null");
  expect(molsymm_fragment_mask(nullptr, nullptr) == -1, "fragment_mask null");
  expect(molsymm_n_free_axes(nullptr, 0, nullptr) == -1, "n_free_axes null");
  expect(molsymm_free_axis(nullptr, 0, 0, nullptr) == -1, "free_axis null");
  expect(molsymm_reconstruct(nullptr, nullptr, nullptr) == -1, "reconstruct null");
  expect(molsymm_subgroup_assembly(nullptr, nullptr, nullptr) == -1, "subgroup_assembly null");
  expect(molsymm_assembly_n_pieces(nullptr, nullptr) == -1, "assembly_n_pieces null");
  expect(molsymm_assembly_piece_n_atoms(nullptr, 0, nullptr) == -1, "assembly_piece_n_atoms null");
  expect(molsymm_assembly_piece_atoms(nullptr, 0, nullptr) == -1, "assembly_piece_atoms null");
  expect(molsymm_assembly_n_free_axes(nullptr, 0, nullptr) == -1, "assembly_n_free_axes null");
  expect(molsymm_assembly_free_axis(nullptr, 0, 0, nullptr) == -1, "assembly_free_axis null");
  expect(molsymm_assemble(nullptr, nullptr, nullptr) == -1, "assemble null");
  expect(molsymm_zmatrix_to_xyz(nullptr, 0, nullptr) == -1, "zmatrix_to_xyz null");
  expect(molsymm_analyze_zmatrix(nullptr, 0, 1e-3, nullptr) == -1, "analyze_zmatrix null");
  expect(molsymm_wyckoff_positions(nullptr, nullptr) == -1, "wyckoff_positions null");
  expect(molsymm_wyckoff_n_sites(nullptr, nullptr) == -1, "wyckoff_n_sites null");
  expect(molsymm_wyckoff_plans_nsites(nullptr, 0, nullptr) == -1, "wyckoff_plans_nsites null");
  expect(molsymm_assemble_from_au(nullptr, nullptr, nullptr, 0, nullptr, nullptr) == -1,
         "assemble_from_au null");
  molsymm_free(nullptr);
  molsymm_free_displacements(nullptr);
  molsymm_free_fragmentation(nullptr);
  molsymm_free_assembly(nullptr);
  molsymm_free_wyckoff_table(nullptr);
  molsymm_free_wyckoff_plans(nullptr);
}

void test_c_api_water() {
  molsymm_analysis* a = nullptr;
  expect(molsymm_analyze(kWaterZ, kWaterXyz, 3, 1e-3, &a) == 0 && a, "analyze water");
  if (!a) {
    return;
  }

  char pg[32] = {};
  int sigma = 0;
  int order = 0;
  int nat = 0;
  int nops = 0;
  int nirr = 0;
  int nsub = 0;
  double rmsd = 1.0;
  expect(molsymm_point_group(a, pg, static_cast<int>(sizeof(pg))) == 0 &&
             std::strcmp(pg, "C2v") == 0,
         "water C2v");
  expect(molsymm_point_group(a, pg, 0) == -1, "point_group buflen=0");
  expect(molsymm_symmetry_number(a, &sigma) == 0 && sigma == 2, "water sigma");
  expect(molsymm_order(a, &order) == 0 && order == 4, "water order");
  expect(molsymm_rmsd(a, &rmsd) == 0 && rmsd < 1e-6, "water rmsd");
  expect(molsymm_n_atoms(a, &nat) == 0 && nat == 3, "water n_atoms");
  expect(molsymm_n_operations(a, &nops) == 0 && nops == 4, "water n_operations");
  expect(molsymm_n_irreps(a, &nirr) == 0 && nirr == 4, "water n_irreps");
  expect(molsymm_n_subgroups(a, &nsub) == 0 && nsub > 0, "water n_subgroups");

  double mat9[9] = {};
  char buf[32] = {};
  expect(molsymm_operation_matrix(a, 0, mat9) == 0, "operation_matrix");
  {
    const double det = mat9[0] * (mat9[4] * mat9[8] - mat9[5] * mat9[7]) -
                       mat9[1] * (mat9[3] * mat9[8] - mat9[5] * mat9[6]) +
                       mat9[2] * (mat9[3] * mat9[7] - mat9[4] * mat9[6]);
    expect(std::abs(std::abs(det) - 1.0) < 1e-8, "operation_matrix det ±1");
  }
  expect(molsymm_operation_matrix(a, 99, mat9) == -1, "operation_matrix OOR");
  {
    const char* want[] = {"A1", "A2", "B1", "B2"};
    for (int i = 0; i < 4; ++i) {
      expect(molsymm_irrep_name(a, i, buf, static_cast<int>(sizeof(buf))) == 0 &&
                 std::strcmp(buf, want[i]) == 0,
             std::string("irrep_name ") + want[i]);
    }
  }
  expect(molsymm_irrep_name(a, 99, buf, static_cast<int>(sizeof(buf))) == -1, "irrep_name OOR");
  expect(molsymm_subgroup_name(a, 0, buf, static_cast<int>(sizeof(buf))) == 0, "subgroup_name");
  expect(molsymm_subgroup_name(a, 99, buf, static_cast<int>(sizeof(buf))) == -1,
         "subgroup_name OOR");
  expect(molsymm_symmetrize(a, &rmsd, mat9) == 0, "symmetrize");
  expect(molsymm_symmetrize(a, nullptr, nullptr) == 0, "symmetrize optional outs");

  molsymm_displacements* d = nullptr;
  expect(molsymm_unique_displacements(a, 1e-3, &d) == 0 && d, "unique_displacements");
  if (d) {
    int nd = 0;
    int atom = -1;
    double dir[3] = {};
    double geom[9] = {};
    expect(molsymm_n_displacements(d, &nd) == 0 && nd == 6, "n_displacements");
    expect(molsymm_displacement_atom(d, 0, &atom) == 0 && atom >= 0 && atom < 3,
           "displacement_atom in range");
    expect(molsymm_displacement_direction(d, 0, dir) == 0, "displacement_direction");
    expect(std::abs(std::hypot(dir[0], dir[1], dir[2]) - 1.0) < 1e-8,
           "displacement_direction unit");
    expect(molsymm_displacement_preserved(d, 0, buf, static_cast<int>(sizeof(buf))) == 0 &&
               buf[0] != '\0',
           "displacement_preserved");
    expect(molsymm_displacement_geometry(d, 0, 1, geom) == 0, "displacement_geometry plus");
    expect(molsymm_displacement_geometry(d, 0, 0, geom) == 0, "displacement_geometry minus");
    expect(molsymm_displacement_atom(d, 99, &atom) == -1, "displacement_atom OOR");
    expect(molsymm_displacement_direction(d, 99, dir) == -1, "displacement_direction OOR");
    expect(molsymm_displacement_preserved(d, 99, buf, static_cast<int>(sizeof(buf))) == -1,
           "displacement_preserved OOR");
    expect(molsymm_displacement_geometry(d, 99, 1, geom) == -1, "displacement_geometry OOR");

    std::vector<double> ep(static_cast<std::size_t>(nd), 0.0);
    std::vector<double> em(static_cast<std::size_t>(nd), 0.0);
    ep[0] = 2e-6;
    std::vector<double> grad(9, 0.0);
    expect(molsymm_gradient_from_energies(d, ep.data(), em.data(), grad.data()) == 0,
           "gradient_from_energies");
    double gabs = 0.0;
    for (double v : grad) {
      gabs += std::abs(v);
    }
    expect(gabs > 0.0, "gradient_from_energies non-zero");
    std::vector<double> g(static_cast<std::size_t>(nd * 9), 0.0);
    std::vector<double> hess(81, 0.0);
    expect(molsymm_hessian_from_gradients(d, g.data(), g.data(), hess.data()) == 0,
           "hessian_from_gradients");
    molsymm_free_displacements(d);
  }

  molsymm_fragmentation* f = nullptr;
  expect(molsymm_make_fragmentation(a, &f) == 0 && f, "make_fragmentation");
  if (f) {
    int nrep = 0;
    int ncopies = 0;
    int atom = -1;
    int copy_n = 0;
    int naxes = 0;
    unsigned char mask[3] = {};
    double axis[3] = {};
    const double rep_xyz[6] = {
        kWaterXyz[0], kWaterXyz[1], kWaterXyz[2], kWaterXyz[3], kWaterXyz[4], kWaterXyz[5]};
    double xyz_out[9] = {};
    expect(molsymm_n_representatives(f, &nrep) == 0 && nrep == 2, "n_representatives");
    expect(molsymm_n_copies(f, &ncopies) == 0 && ncopies > 0, "n_copies");
    expect(molsymm_representative_atom(f, 0, &atom) == 0 && atom == 0, "representative_atom");
    expect(molsymm_representative_atom(f, 99, &atom) == -1, "representative_atom OOR");
    expect(molsymm_copy_n_atoms(f, 0, &copy_n) == 0 && copy_n > 0, "copy_n_atoms");
    expect(molsymm_copy_n_atoms(f, 99, &copy_n) == -1, "copy_n_atoms OOR");
    std::vector<int> atoms(static_cast<std::size_t>(copy_n));
    expect(molsymm_copy_atoms(f, 0, atoms.data()) == 0, "copy_atoms");
    expect(molsymm_copy_atoms(f, 99, atoms.data()) == -1, "copy_atoms OOR");
    expect(molsymm_fragment_mask(f, mask) == 0 && mask[0] && mask[1] && !mask[2], "fragment_mask");
    expect(molsymm_n_free_axes(f, 0, &naxes) == 0, "n_free_axes");
    expect(molsymm_n_free_axes(f, 99, &naxes) == -1, "n_free_axes OOR");
    if (naxes > 0) {
      expect(molsymm_free_axis(f, 0, 0, axis) == 0, "free_axis");
    }
    expect(molsymm_free_axis(f, 0, 99, axis) == -1, "free_axis OOR");
    expect(molsymm_reconstruct(f, rep_xyz, xyz_out) == 0, "reconstruct");
    molsymm_free_fragmentation(f);
  }

  molsymm_assembly* asm_ = nullptr;
  expect(molsymm_subgroup_assembly(a, "C2v", &asm_) == 0 && asm_, "subgroup_assembly");
  if (asm_) {
    int npieces = 0;
    int piece_n = 0;
    int naxes = 0;
    double axis[3] = {};
    const double piece_xyz[6] = {
        kWaterXyz[0], kWaterXyz[1], kWaterXyz[2], kWaterXyz[3], kWaterXyz[4], kWaterXyz[5]};
    double xyz_out[9] = {};
    expect(molsymm_assembly_n_pieces(asm_, &npieces) == 0 && npieces == 1, "assembly_n_pieces");
    expect(molsymm_assembly_piece_n_atoms(asm_, 0, &piece_n) == 0 && piece_n == 2,
           "assembly_piece_n_atoms");
    expect(molsymm_assembly_piece_n_atoms(asm_, 99, &piece_n) == -1, "assembly_piece_n_atoms OOR");
    std::vector<int> atoms(static_cast<std::size_t>(piece_n));
    expect(molsymm_assembly_piece_atoms(asm_, 0, atoms.data()) == 0, "assembly_piece_atoms");
    expect(molsymm_assembly_piece_atoms(asm_, 99, atoms.data()) == -1, "assembly_piece_atoms OOR");
    expect(molsymm_assembly_n_free_axes(asm_, 0, &naxes) == 0, "assembly_n_free_axes");
    expect(molsymm_assembly_n_free_axes(asm_, 99, &naxes) == -1, "assembly_n_free_axes OOR");
    if (naxes > 0) {
      expect(molsymm_assembly_free_axis(asm_, 0, 0, axis) == 0, "assembly_free_axis");
    }
    expect(molsymm_assembly_free_axis(asm_, 0, 99, axis) == -1, "assembly_free_axis OOR");
    expect(molsymm_assemble(asm_, piece_xyz, xyz_out) == 0, "assemble");
    molsymm_free_assembly(asm_);
  }
  expect(molsymm_subgroup_assembly(a, "NotAGroup", &asm_) == -1, "subgroup_assembly unknown");
  expect(molsymm_select_subgroup(a, "Cs") == 0, "select_subgroup Cs");
  expect(molsymm_select_subgroup(a, "NotAGroup") == -1, "select_subgroup unknown");
  molsymm_free(a);
}

void test_c_api_symmetrize_and_zmatrix() {
  const double rattled[9] = {0.0, 0.0, 0.1173, 0.02, 0.7572, -0.4692, 0.0, -0.7572, -0.4692};
  const double forces[9] = {0.1, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, -0.2, 0.0};
  double hess[81] = {};
  for (int i = 0; i < 9; ++i) {
    hess[i * 9 + i] = 1.0;
  }
  double rmsd = 0.0;
  double xyz_sym[9] = {};
  double forces_out[9] = {};
  double hess_out[81] = {};
  expect(
      molsymm_symmetrize_data(
          kWaterZ, rattled, 3, "C2v", 1e-2, forces, hess, &rmsd, xyz_sym, forces_out, hess_out) ==
          0,
      "symmetrize_data full");
  expect(
      molsymm_symmetrize_data(
          kWaterZ, rattled, 3, "C2v", 1e-2, nullptr, nullptr, nullptr, xyz_sym, nullptr, nullptr) ==
          0,
      "symmetrize_data geometry");
  expect(
      molsymm_symmetrize_data(
          kWaterZ, rattled, 3, "C2v", 1e-2, forces, nullptr, &rmsd, xyz_sym, forces_out, nullptr) ==
          0,
      "symmetrize_data forces only");
  expect(molsymm_symmetrize_data(
             kWaterZ, rattled, 3, "C2v", 1e-2, nullptr, hess, &rmsd, xyz_sym, nullptr, hess_out) ==
             0,
         "symmetrize_data hessian only");
  expect(
      molsymm_symmetrize_data(
          kWaterZ, rattled, 3, "Td", 1e-2, nullptr, nullptr, nullptr, xyz_sym, nullptr, nullptr) ==
          -1,
      "symmetrize_data bad subgroup");

  constexpr double kPi = 3.14159265358979323846;
  molsymm_zatom atoms[3] = {};
  atoms[0].atomic_number = 8;
  atoms[1] = {1, 1, 0, 0, 0.957, 0.0, 0.0};
  atoms[2] = {1, 1, 2, 0, 0.957, 104.5 * kPi / 180.0, 0.0};
  double xyz[9] = {};
  expect(molsymm_zmatrix_to_xyz(atoms, 3, xyz) == 0, "zmatrix_to_xyz");
  molsymm_analysis* a = nullptr;
  expect(molsymm_analyze_zmatrix(atoms, 3, 1e-3, &a) == 0 && a, "analyze_zmatrix");
  if (a) {
    char pg[16] = {};
    expect(molsymm_point_group(a, pg, static_cast<int>(sizeof(pg))) == 0 &&
               std::strcmp(pg, "C2v") == 0,
           "zmatrix water C2v");
    molsymm_free(a);
  }
}

void test_c_api_wyckoff() {
  molsymm_wyckoff_table* table = nullptr;
  expect(molsymm_wyckoff_positions("Oh", &table) == 0 && table, "wyckoff_positions Oh");
  if (table) {
    int ns = 0;
    expect(molsymm_wyckoff_n_sites(table, &ns) == 0 && ns >= 5, "wyckoff_n_sites");
    char label[16] = {};
    char sym[16] = {};
    int mult = 0;
    int dof = 0;
    expect(molsymm_wyckoff_site_label(table, 0, label, 16) == 0, "wyckoff_site_label");
    expect(molsymm_wyckoff_site_symmetry(table, 0, sym, 16) == 0, "wyckoff_site_symmetry");
    expect(molsymm_wyckoff_site_multiplicity(table, 0, &mult) == 0 && mult == 1,
           "wyckoff_site_multiplicity origin");
    expect(molsymm_wyckoff_site_dof(table, 0, &dof) == 0 && dof == 0, "wyckoff_site_dof origin");
    molsymm_free_wyckoff_table(table);
  }
  expect(molsymm_wyckoff_positions("Cinfv", &table) == -1, "wyckoff_positions infinite");

  molsymm_wyckoff_plans* plans = nullptr;
  expect(molsymm_wyckoff_plans_nsites("Oh", 12, &plans) == 0 && plans, "wyckoff_plans_nsites");
  if (plans) {
    int np = 0;
    expect(molsymm_wyckoff_n_plans(plans, &np) == 0 && np >= 2, "wyckoff_n_plans");
    int norb = 0;
    expect(molsymm_wyckoff_plan_n_orbits(plans, 0, &norb) == 0 && norb >= 1,
           "wyckoff_plan_n_orbits");
    char olabel[16] = {};
    expect(molsymm_wyckoff_plan_orbit_label(plans, 0, 0, olabel, 16) == 0,
           "wyckoff_plan_orbit_label");
    molsymm_free_wyckoff_plans(plans);
  }

  const int z[2] = {16, 9};
  const double au[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 1.56};
  const char* labels[2] = {"1a", "6a"};
  molsymm_analysis* a = nullptr;
  expect(molsymm_assemble_from_au("Oh", z, au, 2, labels, &a) == 0 && a, "assemble_from_au SF6");
  if (a) {
    char pg[16] = {};
    int n = 0;
    expect(molsymm_point_group(a, pg, 16) == 0 && std::strcmp(pg, "Oh") == 0,
           "assemble_from_au Oh");
    expect(molsymm_n_atoms(a, &n) == 0 && n == 7, "assemble_from_au n=7");
    molsymm_free(a);
  }
}

} // namespace

void test_c_api() {
  test_c_api_nulls();
  test_c_api_water();
  test_c_api_symmetrize_and_zmatrix();
  test_c_api_wyckoff();
  if (g_failed) {
    throw std::runtime_error(std::to_string(g_failed) + " c-api failures (" +
                             std::to_string(g_passed) + " checks passed)");
  }
  std::cout << "c-api checks passed (" << g_passed << ")\n";
}
