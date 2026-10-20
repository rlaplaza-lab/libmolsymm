#ifndef MOLSYMM_H
#define MOLSYMM_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct molsymm_analysis molsymm_analysis;
typedef struct molsymm_displacements molsymm_displacements;

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

#ifdef __cplusplus
}
#endif

#endif /* MOLSYMM_H */
