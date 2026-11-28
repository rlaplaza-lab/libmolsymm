#ifndef MOLSYMM_H
#define MOLSYMM_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct molsymm_analysis molsymm_analysis;
typedef struct molsymm_displacements molsymm_displacements;
typedef struct molsymm_fragmentation molsymm_fragmentation;
typedef struct molsymm_assembly molsymm_assembly;
typedef struct molsymm_wyckoff_table molsymm_wyckoff_table;
typedef struct molsymm_wyckoff_plans molsymm_wyckoff_plans;

const char* molsymm_last_error(void);

int molsymm_analyze(const int* atomic_numbers,
                    const double* xyz,
                    int n,
                    double primary_tol,
                    molsymm_analysis** out);

/* One Z-matrix row. Indices are 1-based (0 = unused). Distances are Å, angles radians.
   Field order matches molsymm::ZMatrixAtom. */
typedef struct molsymm_zatom {
  int atomic_number;
  int bond_to;
  int angle_to;
  int dihedral_to;
  double distance;
  double angle;
  double dihedral;
} molsymm_zatom;

int molsymm_zmatrix_to_xyz(const molsymm_zatom* atoms, int n, double* xyz_out);
int molsymm_analyze_zmatrix(const molsymm_zatom* atoms,
                            int n,
                            double primary_tol,
                            molsymm_analysis** out);

void molsymm_free(molsymm_analysis* a);

int molsymm_point_group(const molsymm_analysis* a, char* buf, int buflen);
int molsymm_symmetry_number(const molsymm_analysis* a, int* out);
int molsymm_order(const molsymm_analysis* a, int* out);
int molsymm_rmsd(const molsymm_analysis* a, double* out);
int molsymm_n_atoms(const molsymm_analysis* a, int* out);
int molsymm_n_operations(const molsymm_analysis* a, int* out);
int molsymm_operation_matrix(const molsymm_analysis* a, int index, double* mat9);
int molsymm_n_irreps(const molsymm_analysis* a, int* out);
int molsymm_irrep_name(const molsymm_analysis* a, int index, char* buf, int buflen);
int molsymm_n_subgroups(const molsymm_analysis* a, int* out);
int molsymm_subgroup_name(const molsymm_analysis* a, int index, char* buf, int buflen);
int molsymm_select_subgroup(molsymm_analysis* a, const char* name);
int molsymm_symmetrize(molsymm_analysis* a, double* rmsd, double* xyz_out);

/* Project geometry and optional forces / Hessian onto point_group.
   xyz is length 3*n. forces (if non-NULL) is length 3*n. hessian (if non-NULL)
   is row-major 3n×3n. Matching outputs may be NULL when the input is NULL. */
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
                            double* hessian_out);

int molsymm_unique_displacements(const molsymm_analysis* a,
                                 double step,
                                 molsymm_displacements** out);
void molsymm_free_displacements(molsymm_displacements* d);
int molsymm_n_displacements(const molsymm_displacements* d, int* out);
int molsymm_displacement_atom(const molsymm_displacements* d, int index, int* atom);
int molsymm_displacement_direction(const molsymm_displacements* d, int index, double* xyz);
int molsymm_displacement_preserved(const molsymm_displacements* d,
                                   int index,
                                   char* buf,
                                   int buflen);
int molsymm_displacement_geometry(const molsymm_displacements* d,
                                  int index,
                                  int plus,
                                  double* xyz_out);
int molsymm_gradient_from_energies(const molsymm_displacements* d,
                                   const double* e_plus,
                                   const double* e_minus,
                                   double* grad_out);
int molsymm_hessian_from_gradients(const molsymm_displacements* d,
                                   const double* g_plus,
                                   const double* g_minus,
                                   double* hess_out);

int molsymm_make_fragmentation(const molsymm_analysis* a, molsymm_fragmentation** out);
void molsymm_free_fragmentation(molsymm_fragmentation* f);
int molsymm_n_representatives(const molsymm_fragmentation* f, int* out);
int molsymm_representative_atom(const molsymm_fragmentation* f, int index, int* atom);
int molsymm_n_copies(const molsymm_fragmentation* f, int* out);
int molsymm_copy_n_atoms(const molsymm_fragmentation* f, int copy, int* out);
int molsymm_copy_atoms(const molsymm_fragmentation* f, int copy, int* atoms_out);
int molsymm_fragment_mask(const molsymm_fragmentation* f, unsigned char* mask_out);
int molsymm_n_free_axes(const molsymm_fragmentation* f, int site, int* out);
int molsymm_free_axis(const molsymm_fragmentation* f, int site, int axis, double* xyz);
int molsymm_reconstruct(const molsymm_fragmentation* f, const double* rep_xyz, double* xyz_out);

int molsymm_subgroup_assembly(const molsymm_analysis* a,
                              const char* subgroup,
                              molsymm_assembly** out);
void molsymm_free_assembly(molsymm_assembly* asm_);
int molsymm_assembly_n_pieces(const molsymm_assembly* asm_, int* out);
int molsymm_assembly_piece_n_atoms(const molsymm_assembly* asm_, int piece, int* out);
int molsymm_assembly_piece_atoms(const molsymm_assembly* asm_, int piece, int* atoms_out);
int molsymm_assembly_n_free_axes(const molsymm_assembly* asm_, int site, int* out);
int molsymm_assembly_free_axis(const molsymm_assembly* asm_, int site, int axis, double* xyz);
int molsymm_assemble(const molsymm_assembly* asm_, const double* piece_xyz, double* xyz_out);

/* Constructive Wyckoff positions / plans / assemble (finite groups only). */
int molsymm_wyckoff_positions(const char* point_group, molsymm_wyckoff_table** out);
void molsymm_free_wyckoff_table(molsymm_wyckoff_table* t);
int molsymm_wyckoff_n_sites(const molsymm_wyckoff_table* t, int* out);
int molsymm_wyckoff_site_label(const molsymm_wyckoff_table* t, int index, char* buf, int buflen);
int molsymm_wyckoff_site_symmetry(const molsymm_wyckoff_table* t, int index, char* buf, int buflen);
int molsymm_wyckoff_site_multiplicity(const molsymm_wyckoff_table* t, int index, int* out);
int molsymm_wyckoff_site_dof(const molsymm_wyckoff_table* t, int index, int* out);

int molsymm_wyckoff_plans_composition(const char* point_group,
                                      const int* atomic_numbers,
                                      int n,
                                      molsymm_wyckoff_plans** out);
int molsymm_wyckoff_plans_nsites(const char* point_group, int n_sites, molsymm_wyckoff_plans** out);
void molsymm_free_wyckoff_plans(molsymm_wyckoff_plans* p);
int molsymm_wyckoff_n_plans(const molsymm_wyckoff_plans* p, int* out);
int molsymm_wyckoff_plan_n_orbits(const molsymm_wyckoff_plans* p, int plan, int* out);
int molsymm_wyckoff_plan_orbit_label(
    const molsymm_wyckoff_plans* p, int plan, int orbit, char* buf, int buflen);
int molsymm_wyckoff_plan_orbit_z(const molsymm_wyckoff_plans* p, int plan, int orbit, int* out);

/* site_labels[i] is the Wyckoff label for orbit i (length n_orbits). */
int molsymm_assemble_from_au(const char* point_group,
                             const int* representative_numbers,
                             const double* au_xyz,
                             int n_orbits,
                             const char** site_labels,
                             molsymm_analysis** out);

#ifdef __cplusplus
}
#endif

#endif /* MOLSYMM_H */
