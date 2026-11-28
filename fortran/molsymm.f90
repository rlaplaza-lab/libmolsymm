module molsymm
  use, intrinsic :: iso_c_binding
  implicit none
  private
  public :: ms_analyze, ms_release, ms_point_group_name, ms_symmetry_number, ms_order, ms_rmsd
  public :: ms_symmetrize, ms_symmetrize_data, ms_unique_displacements, ms_n_displacements
  public :: ms_displacement_geometry, ms_hessian_from_gradients, ms_release_displacements
  public :: ms_last_error, ms_zatom, ms_zmatrix_to_xyz, ms_analyze_zmatrix
  public :: ms_fragmentation, ms_release_fragmentation, ms_n_representatives
  public :: ms_representative_atom, ms_n_copies, ms_copy_atoms, ms_fragment_mask
  public :: ms_reconstruct, ms_subgroup_assembly, ms_release_assembly
  public :: ms_assembly_n_pieces, ms_assembly_piece_atoms, ms_assemble
  public :: ms_wyckoff_positions, ms_release_wyckoff_table, ms_wyckoff_n_sites
  public :: ms_wyckoff_site_label, ms_wyckoff_site_multiplicity, ms_wyckoff_site_dof
  public :: ms_wyckoff_plans_nsites, ms_release_wyckoff_plans, ms_wyckoff_n_plans
  public :: ms_wyckoff_plan_n_orbits, ms_wyckoff_plan_orbit_label
  public :: ms_assemble_from_au

  type, bind(c) :: ms_zatom
    integer(c_int) :: atomic_number
    integer(c_int) :: bond_to
    integer(c_int) :: angle_to
    integer(c_int) :: dihedral_to
    real(c_double) :: distance
    real(c_double) :: angle
    real(c_double) :: dihedral
  end type

interface
  function molsymm_last_error() bind(c, name="molsymm_last_error")
    import :: c_ptr
    type(c_ptr) :: molsymm_last_error
  end function

  function molsymm_analyze(z, xyz, n, primary, out) bind(c, name="molsymm_analyze")
    import :: c_int, c_double, c_ptr
    integer(c_int), intent(in) :: z(*)
    real(c_double), intent(in) :: xyz(*)
    integer(c_int), value :: n
    real(c_double), value :: primary
    type(c_ptr), intent(out) :: out
    integer(c_int) :: molsymm_analyze
  end function

  function molsymm_zmatrix_to_xyz(atoms, n, xyz_out) bind(c, name="molsymm_zmatrix_to_xyz")
    import :: c_int, c_double, ms_zatom
    type(ms_zatom), intent(in) :: atoms(*)
    integer(c_int), value :: n
    real(c_double), intent(out) :: xyz_out(*)
    integer(c_int) :: molsymm_zmatrix_to_xyz
  end function

  function molsymm_analyze_zmatrix(atoms, n, primary, out) bind(c, name="molsymm_analyze_zmatrix")
    import :: c_int, c_double, c_ptr, ms_zatom
    type(ms_zatom), intent(in) :: atoms(*)
    integer(c_int), value :: n
    real(c_double), value :: primary
    type(c_ptr), intent(out) :: out
    integer(c_int) :: molsymm_analyze_zmatrix
  end function

  subroutine molsymm_free(a) bind(c, name="molsymm_free")
    import :: c_ptr
    type(c_ptr), value :: a
  end subroutine

  function molsymm_point_group(a, buf, buflen) bind(c, name="molsymm_point_group")
    import :: c_int, c_ptr, c_char
    type(c_ptr), value :: a
    character(kind=c_char), intent(out) :: buf(*)
    integer(c_int), value :: buflen
    integer(c_int) :: molsymm_point_group
  end function

  function molsymm_symmetry_number(a, out) bind(c, name="molsymm_symmetry_number")
    import :: c_int, c_ptr
    type(c_ptr), value :: a
    integer(c_int), intent(out) :: out
    integer(c_int) :: molsymm_symmetry_number
  end function

  function molsymm_order(a, out) bind(c, name="molsymm_order")
    import :: c_int, c_ptr
    type(c_ptr), value :: a
    integer(c_int), intent(out) :: out
    integer(c_int) :: molsymm_order
  end function

  function molsymm_rmsd(a, out) bind(c, name="molsymm_rmsd")
    import :: c_int, c_ptr, c_double
    type(c_ptr), value :: a
    real(c_double), intent(out) :: out
    integer(c_int) :: molsymm_rmsd
  end function

  function molsymm_symmetrize(a, rmsd, xyz_out) bind(c, name="molsymm_symmetrize")
    import :: c_int, c_ptr, c_double
    type(c_ptr), value :: a
    real(c_double), intent(out) :: rmsd
    real(c_double), intent(out) :: xyz_out(*)
    integer(c_int) :: molsymm_symmetrize
  end function

  function molsymm_symmetrize_data(z, xyz, n, point_group, primary, forces, hessian, &
      rmsd, xyz_out, forces_out, hessian_out) bind(c, name="molsymm_symmetrize_data")
    import :: c_int, c_double, c_ptr, c_char
    integer(c_int), intent(in) :: z(*)
    real(c_double), intent(in) :: xyz(*)
    integer(c_int), value :: n
    character(kind=c_char), intent(in) :: point_group(*)
    real(c_double), value :: primary
    type(c_ptr), value :: forces, hessian
    real(c_double), intent(out) :: rmsd
    real(c_double), intent(out) :: xyz_out(*)
    type(c_ptr), value :: forces_out, hessian_out
    integer(c_int) :: molsymm_symmetrize_data
  end function

  function molsymm_unique_displacements(a, step, out) bind(c, name="molsymm_unique_displacements")
    import :: c_int, c_ptr, c_double
    type(c_ptr), value :: a
    real(c_double), value :: step
    type(c_ptr), intent(out) :: out
    integer(c_int) :: molsymm_unique_displacements
  end function

  subroutine molsymm_free_displacements(d) bind(c, name="molsymm_free_displacements")
    import :: c_ptr
    type(c_ptr), value :: d
  end subroutine

  function molsymm_n_displacements(d, out) bind(c, name="molsymm_n_displacements")
    import :: c_int, c_ptr
    type(c_ptr), value :: d
    integer(c_int), intent(out) :: out
    integer(c_int) :: molsymm_n_displacements
  end function

  function molsymm_displacement_geometry(d, index, plus, xyz_out) &
      bind(c, name="molsymm_displacement_geometry")
    import :: c_int, c_ptr, c_double
    type(c_ptr), value :: d
    integer(c_int), value :: index, plus
    real(c_double), intent(out) :: xyz_out(*)
    integer(c_int) :: molsymm_displacement_geometry
  end function

  function molsymm_hessian_from_gradients(d, g_plus, g_minus, hess_out) &
      bind(c, name="molsymm_hessian_from_gradients")
    import :: c_int, c_ptr, c_double
    type(c_ptr), value :: d
    real(c_double), intent(in) :: g_plus(*), g_minus(*)
    real(c_double), intent(out) :: hess_out(*)
    integer(c_int) :: molsymm_hessian_from_gradients
  end function

  function molsymm_make_fragmentation(a, out) bind(c, name="molsymm_make_fragmentation")
    import :: c_int, c_ptr
    type(c_ptr), value :: a
    type(c_ptr), intent(out) :: out
    integer(c_int) :: molsymm_make_fragmentation
  end function

  subroutine molsymm_free_fragmentation(f) bind(c, name="molsymm_free_fragmentation")
    import :: c_ptr
    type(c_ptr), value :: f
  end subroutine

  function molsymm_n_representatives(f, out) bind(c, name="molsymm_n_representatives")
    import :: c_int, c_ptr
    type(c_ptr), value :: f
    integer(c_int), intent(out) :: out
    integer(c_int) :: molsymm_n_representatives
  end function

  function molsymm_representative_atom(f, index, atom) bind(c, name="molsymm_representative_atom")
    import :: c_int, c_ptr
    type(c_ptr), value :: f
    integer(c_int), value :: index
    integer(c_int), intent(out) :: atom
    integer(c_int) :: molsymm_representative_atom
  end function

  function molsymm_n_copies(f, out) bind(c, name="molsymm_n_copies")
    import :: c_int, c_ptr
    type(c_ptr), value :: f
    integer(c_int), intent(out) :: out
    integer(c_int) :: molsymm_n_copies
  end function

  function molsymm_copy_n_atoms(f, copy, out) bind(c, name="molsymm_copy_n_atoms")
    import :: c_int, c_ptr
    type(c_ptr), value :: f
    integer(c_int), value :: copy
    integer(c_int), intent(out) :: out
    integer(c_int) :: molsymm_copy_n_atoms
  end function

  function molsymm_copy_atoms(f, copy, atoms_out) bind(c, name="molsymm_copy_atoms")
    import :: c_int, c_ptr
    type(c_ptr), value :: f
    integer(c_int), value :: copy
    integer(c_int), intent(out) :: atoms_out(*)
    integer(c_int) :: molsymm_copy_atoms
  end function

  function molsymm_fragment_mask(f, mask_out) bind(c, name="molsymm_fragment_mask")
    import :: c_int, c_ptr, c_signed_char
    type(c_ptr), value :: f
    integer(c_signed_char), intent(out) :: mask_out(*)
    integer(c_int) :: molsymm_fragment_mask
  end function

  function molsymm_reconstruct(f, rep_xyz, xyz_out) bind(c, name="molsymm_reconstruct")
    import :: c_int, c_ptr, c_double
    type(c_ptr), value :: f
    real(c_double), intent(in) :: rep_xyz(*)
    real(c_double), intent(out) :: xyz_out(*)
    integer(c_int) :: molsymm_reconstruct
  end function

  function molsymm_subgroup_assembly(a, subgroup, out) bind(c, name="molsymm_subgroup_assembly")
    import :: c_int, c_ptr, c_char
    type(c_ptr), value :: a
    character(kind=c_char), intent(in) :: subgroup(*)
    type(c_ptr), intent(out) :: out
    integer(c_int) :: molsymm_subgroup_assembly
  end function

  subroutine molsymm_free_assembly(asm) bind(c, name="molsymm_free_assembly")
    import :: c_ptr
    type(c_ptr), value :: asm
  end subroutine

  function molsymm_assembly_n_pieces(asm, out) bind(c, name="molsymm_assembly_n_pieces")
    import :: c_int, c_ptr
    type(c_ptr), value :: asm
    integer(c_int), intent(out) :: out
    integer(c_int) :: molsymm_assembly_n_pieces
  end function

  function molsymm_assembly_piece_n_atoms(asm, piece, out) &
      bind(c, name="molsymm_assembly_piece_n_atoms")
    import :: c_int, c_ptr
    type(c_ptr), value :: asm
    integer(c_int), value :: piece
    integer(c_int), intent(out) :: out
    integer(c_int) :: molsymm_assembly_piece_n_atoms
  end function

  function molsymm_assembly_piece_atoms(asm, piece, atoms_out) &
      bind(c, name="molsymm_assembly_piece_atoms")
    import :: c_int, c_ptr
    type(c_ptr), value :: asm
    integer(c_int), value :: piece
    integer(c_int), intent(out) :: atoms_out(*)
    integer(c_int) :: molsymm_assembly_piece_atoms
  end function

  function molsymm_assemble(asm, piece_xyz, xyz_out) bind(c, name="molsymm_assemble")
    import :: c_int, c_ptr, c_double
    type(c_ptr), value :: asm
    real(c_double), intent(in) :: piece_xyz(*)
    real(c_double), intent(out) :: xyz_out(*)
    integer(c_int) :: molsymm_assemble
  end function

  function molsymm_wyckoff_positions(pg, out) bind(c, name="molsymm_wyckoff_positions")
    import :: c_int, c_ptr, c_char
    character(kind=c_char), intent(in) :: pg(*)
    type(c_ptr), intent(out) :: out
    integer(c_int) :: molsymm_wyckoff_positions
  end function

  subroutine molsymm_free_wyckoff_table(t) bind(c, name="molsymm_free_wyckoff_table")
    import :: c_ptr
    type(c_ptr), value :: t
  end subroutine

  function molsymm_wyckoff_n_sites(t, out) bind(c, name="molsymm_wyckoff_n_sites")
    import :: c_int, c_ptr
    type(c_ptr), value :: t
    integer(c_int), intent(out) :: out
    integer(c_int) :: molsymm_wyckoff_n_sites
  end function

  function molsymm_wyckoff_site_label(t, index, buf, buflen) &
      bind(c, name="molsymm_wyckoff_site_label")
    import :: c_int, c_ptr, c_char
    type(c_ptr), value :: t
    integer(c_int), value :: index, buflen
    character(kind=c_char), intent(out) :: buf(*)
    integer(c_int) :: molsymm_wyckoff_site_label
  end function

  function molsymm_wyckoff_site_multiplicity(t, index, out) &
      bind(c, name="molsymm_wyckoff_site_multiplicity")
    import :: c_int, c_ptr
    type(c_ptr), value :: t
    integer(c_int), value :: index
    integer(c_int), intent(out) :: out
    integer(c_int) :: molsymm_wyckoff_site_multiplicity
  end function

  function molsymm_wyckoff_site_dof(t, index, out) bind(c, name="molsymm_wyckoff_site_dof")
    import :: c_int, c_ptr
    type(c_ptr), value :: t
    integer(c_int), value :: index
    integer(c_int), intent(out) :: out
    integer(c_int) :: molsymm_wyckoff_site_dof
  end function

  function molsymm_wyckoff_plans_nsites(pg, n_sites, out) &
      bind(c, name="molsymm_wyckoff_plans_nsites")
    import :: c_int, c_ptr, c_char
    character(kind=c_char), intent(in) :: pg(*)
    integer(c_int), value :: n_sites
    type(c_ptr), intent(out) :: out
    integer(c_int) :: molsymm_wyckoff_plans_nsites
  end function

  subroutine molsymm_free_wyckoff_plans(p) bind(c, name="molsymm_free_wyckoff_plans")
    import :: c_ptr
    type(c_ptr), value :: p
  end subroutine

  function molsymm_wyckoff_n_plans(p, out) bind(c, name="molsymm_wyckoff_n_plans")
    import :: c_int, c_ptr
    type(c_ptr), value :: p
    integer(c_int), intent(out) :: out
    integer(c_int) :: molsymm_wyckoff_n_plans
  end function

  function molsymm_wyckoff_plan_n_orbits(p, plan, out) &
      bind(c, name="molsymm_wyckoff_plan_n_orbits")
    import :: c_int, c_ptr
    type(c_ptr), value :: p
    integer(c_int), value :: plan
    integer(c_int), intent(out) :: out
    integer(c_int) :: molsymm_wyckoff_plan_n_orbits
  end function

  function molsymm_wyckoff_plan_orbit_label(p, plan, orbit, buf, buflen) &
      bind(c, name="molsymm_wyckoff_plan_orbit_label")
    import :: c_int, c_ptr, c_char
    type(c_ptr), value :: p
    integer(c_int), value :: plan, orbit, buflen
    character(kind=c_char), intent(out) :: buf(*)
    integer(c_int) :: molsymm_wyckoff_plan_orbit_label
  end function

  function molsymm_assemble_from_au(pg, z, au_xyz, n_orbits, labels, out) &
      bind(c, name="molsymm_assemble_from_au")
    import :: c_int, c_ptr, c_char, c_double
    character(kind=c_char), intent(in) :: pg(*)
    integer(c_int), intent(in) :: z(*)
    real(c_double), intent(in) :: au_xyz(*)
    integer(c_int), value :: n_orbits
    type(c_ptr), intent(in) :: labels(*)
    type(c_ptr), intent(out) :: out
    integer(c_int) :: molsymm_assemble_from_au
  end function
end interface

contains

  function ms_last_error() result(msg)
    character(len=:), allocatable :: msg
    type(c_ptr) :: p
    character(kind=c_char), pointer :: fptr(:)
    integer :: i, n
    p = molsymm_last_error()
    if (.not. c_associated(p)) then
      msg = ""
      return
    end if
    call c_f_pointer(p, fptr, [4096])
    n = 0
    do i = 1, 4096
      if (fptr(i) == c_null_char) exit
      n = i
    end do
    allocate(character(len=n) :: msg)
    do i = 1, n
      msg(i:i) = fptr(i)
    end do
  end function

  subroutine ms_analyze(atomic_numbers, xyz, handle, status, primary)
    integer(c_int), intent(in) :: atomic_numbers(:)
    real(c_double), intent(in) :: xyz(:,:)  ! (3,n)
    type(c_ptr), intent(out) :: handle
    integer(c_int), intent(out) :: status
    real(c_double), intent(in), optional :: primary
    real(c_double) :: prim
    real(c_double), allocatable :: flat(:)
    integer :: n, i
    prim = 1.0d-3
    if (present(primary)) prim = primary
    n = size(atomic_numbers)
    allocate(flat(3*n))
    do i = 1, n
      flat(3*(i-1)+1:3*(i-1)+3) = xyz(:, i)
    end do
    status = molsymm_analyze(atomic_numbers, flat, int(n, c_int), prim, handle)
  end subroutine

  subroutine ms_zmatrix_to_xyz(atoms, xyz, status)
    type(ms_zatom), intent(in) :: atoms(:)
    real(c_double), intent(out) :: xyz(:,:)  ! (3,n)
    integer(c_int), intent(out) :: status
    type(ms_zatom), allocatable :: packed(:)
    real(c_double), allocatable :: flat(:)
    integer :: n, i
    n = size(atoms)
    allocate(packed(n), flat(3*n))
    packed = atoms
    status = molsymm_zmatrix_to_xyz(packed, int(n, c_int), flat)
    do i = 1, n
      xyz(:, i) = flat(3*(i-1)+1:3*(i-1)+3)
    end do
  end subroutine

  subroutine ms_analyze_zmatrix(atoms, handle, status, primary)
    type(ms_zatom), intent(in) :: atoms(:)
    type(c_ptr), intent(out) :: handle
    integer(c_int), intent(out) :: status
    real(c_double), intent(in), optional :: primary
    type(ms_zatom), allocatable :: packed(:)
    real(c_double) :: prim
    integer :: n
    prim = 1.0d-3
    if (present(primary)) prim = primary
    n = size(atoms)
    allocate(packed(n))
    packed = atoms
    status = molsymm_analyze_zmatrix(packed, int(n, c_int), prim, handle)
  end subroutine

  subroutine ms_release(handle)
    type(c_ptr), intent(inout) :: handle
    if (c_associated(handle)) call molsymm_free(handle)
    handle = c_null_ptr
  end subroutine

  subroutine ms_point_group_name(handle, name, status)
    type(c_ptr), intent(in) :: handle
    character(len=*), intent(out) :: name
    integer(c_int), intent(out) :: status
    character(kind=c_char) :: buf(64)
    integer :: i
    status = molsymm_point_group(handle, buf, 64_c_int)
    name = ""
    do i = 1, 64
      if (buf(i) == c_null_char) exit
      name(i:i) = buf(i)
    end do
  end subroutine

  subroutine ms_symmetry_number(handle, sigma, status)
    type(c_ptr), intent(in) :: handle
    integer(c_int), intent(out) :: sigma, status
    status = molsymm_symmetry_number(handle, sigma)
  end subroutine

  subroutine ms_order(handle, order, status)
    type(c_ptr), intent(in) :: handle
    integer(c_int), intent(out) :: order, status
    status = molsymm_order(handle, order)
  end subroutine

  subroutine ms_rmsd(handle, rmsd, status)
    type(c_ptr), intent(in) :: handle
    real(c_double), intent(out) :: rmsd
    integer(c_int), intent(out) :: status
    status = molsymm_rmsd(handle, rmsd)
  end subroutine

  subroutine ms_symmetrize(handle, xyz_out, rmsd, status)
    type(c_ptr), intent(in) :: handle
    real(c_double), intent(out) :: xyz_out(:,:)
    real(c_double), intent(out) :: rmsd
    integer(c_int), intent(out) :: status
    real(c_double), allocatable :: flat(:)
    integer :: n, i
    n = size(xyz_out, 2)
    allocate(flat(3*n))
    status = molsymm_symmetrize(handle, rmsd, flat)
    do i = 1, n
      xyz_out(:, i) = flat(3*(i-1)+1:3*(i-1)+3)
    end do
  end subroutine

  subroutine ms_symmetrize_data(atomic_numbers, xyz, point_group, xyz_out, rmsd, status, &
                                forces, hessian, forces_out, hessian_out, primary)
    integer(c_int), intent(in) :: atomic_numbers(:)
    real(c_double), intent(in) :: xyz(:,:)  ! (3,n)
    character(len=*), intent(in) :: point_group
    real(c_double), intent(out) :: xyz_out(:,:)  ! (3,n)
    real(c_double), intent(out) :: rmsd
    integer(c_int), intent(out) :: status
    real(c_double), intent(in), optional :: forces(:,:)   ! (3,n)
    real(c_double), intent(in), optional :: hessian(:,:)  ! (3n,3n)
    real(c_double), intent(out), optional :: forces_out(:,:)
    real(c_double), intent(out), optional :: hessian_out(:,:)
    real(c_double), intent(in), optional :: primary
    real(c_double) :: prim
    real(c_double), allocatable, target :: flat(:), flat_out(:)
    real(c_double), allocatable, target :: f_in(:), f_out(:), h_in(:), h_out(:)
    character(kind=c_char), allocatable :: pg(:)
    type(c_ptr) :: forces_ptr, hessian_ptr, forces_out_ptr, hessian_out_ptr
    integer :: n, i, n3, j, ntrim
    prim = 1.0d-3
    if (present(primary)) prim = primary
    n = size(atomic_numbers)
    n3 = 3 * n
    ntrim = len_trim(point_group)
    allocate(flat(n3), flat_out(n3), pg(ntrim + 1))
    do i = 1, n
      flat(3*(i-1)+1:3*(i-1)+3) = xyz(:, i)
    end do
    do i = 1, ntrim
      pg(i) = point_group(i:i)
    end do
    pg(ntrim + 1) = c_null_char

    forces_ptr = c_null_ptr
    hessian_ptr = c_null_ptr
    forces_out_ptr = c_null_ptr
    hessian_out_ptr = c_null_ptr
    if (present(forces)) then
      allocate(f_in(n3))
      do i = 1, n
        f_in(3*(i-1)+1:3*(i-1)+3) = forces(:, i)
      end do
      forces_ptr = c_loc(f_in(1))
      if (present(forces_out)) then
        allocate(f_out(n3))
        forces_out_ptr = c_loc(f_out(1))
      end if
    end if
    if (present(hessian)) then
      allocate(h_in(n3 * n3))
      do i = 1, n3
        do j = 1, n3
          h_in((i - 1) * n3 + j) = hessian(i, j)
        end do
      end do
      hessian_ptr = c_loc(h_in(1))
      if (present(hessian_out)) then
        allocate(h_out(n3 * n3))
        hessian_out_ptr = c_loc(h_out(1))
      end if
    end if

    status = molsymm_symmetrize_data(atomic_numbers, flat, int(n, c_int), pg, prim, &
                                     forces_ptr, hessian_ptr, rmsd, flat_out, &
                                     forces_out_ptr, hessian_out_ptr)
    do i = 1, n
      xyz_out(:, i) = flat_out(3*(i-1)+1:3*(i-1)+3)
    end do
    if (present(forces) .and. present(forces_out)) then
      do i = 1, n
        forces_out(:, i) = f_out(3*(i-1)+1:3*(i-1)+3)
      end do
    end if
    if (present(hessian) .and. present(hessian_out)) then
      do i = 1, n3
        do j = 1, n3
          hessian_out(i, j) = h_out((i - 1) * n3 + j)
        end do
      end do
    end if
  end subroutine

  subroutine ms_unique_displacements(handle, step, dhandle, status)
    type(c_ptr), intent(in) :: handle
    real(c_double), intent(in) :: step
    type(c_ptr), intent(out) :: dhandle
    integer(c_int), intent(out) :: status
    status = molsymm_unique_displacements(handle, step, dhandle)
  end subroutine

  subroutine ms_release_displacements(dhandle)
    type(c_ptr), intent(inout) :: dhandle
    if (c_associated(dhandle)) call molsymm_free_displacements(dhandle)
    dhandle = c_null_ptr
  end subroutine

  subroutine ms_n_displacements(dhandle, n, status)
    type(c_ptr), intent(in) :: dhandle
    integer(c_int), intent(out) :: n, status
    status = molsymm_n_displacements(dhandle, n)
  end subroutine

  subroutine ms_displacement_geometry(dhandle, index, plus, xyz_out, status)
    type(c_ptr), intent(in) :: dhandle
    integer(c_int), intent(in) :: index, plus
    real(c_double), intent(out) :: xyz_out(:,:)
    integer(c_int), intent(out) :: status
    real(c_double), allocatable :: flat(:)
    integer :: n, i
    n = size(xyz_out, 2)
    allocate(flat(3*n))
    status = molsymm_displacement_geometry(dhandle, index, plus, flat)
    do i = 1, n
      xyz_out(:, i) = flat(3*(i-1)+1:3*(i-1)+3)
    end do
  end subroutine

  subroutine ms_hessian_from_gradients(dhandle, g_plus, g_minus, hess_out, status)
    type(c_ptr), intent(in) :: dhandle
    real(c_double), intent(in) :: g_plus(:), g_minus(:)
    real(c_double), intent(out) :: hess_out(:)
    integer(c_int), intent(out) :: status
    status = molsymm_hessian_from_gradients(dhandle, g_plus, g_minus, hess_out)
  end subroutine

  subroutine ms_fragmentation(handle, fhandle, status)
    type(c_ptr), intent(in) :: handle
    type(c_ptr), intent(out) :: fhandle
    integer(c_int), intent(out) :: status
    status = molsymm_make_fragmentation(handle, fhandle)
  end subroutine

  subroutine ms_release_fragmentation(fhandle)
    type(c_ptr), intent(inout) :: fhandle
    if (c_associated(fhandle)) call molsymm_free_fragmentation(fhandle)
    fhandle = c_null_ptr
  end subroutine

  subroutine ms_n_representatives(fhandle, n, status)
    type(c_ptr), intent(in) :: fhandle
    integer(c_int), intent(out) :: n, status
    status = molsymm_n_representatives(fhandle, n)
  end subroutine

  subroutine ms_representative_atom(fhandle, index, atom, status)
    type(c_ptr), intent(in) :: fhandle
    integer(c_int), intent(in) :: index
    integer(c_int), intent(out) :: atom, status
    ! Fortran indices are 1-based; C is 0-based.
    status = molsymm_representative_atom(fhandle, index - 1_c_int, atom)
    if (status == 0) atom = atom + 1_c_int
  end subroutine

  subroutine ms_n_copies(fhandle, n, status)
    type(c_ptr), intent(in) :: fhandle
    integer(c_int), intent(out) :: n, status
    status = molsymm_n_copies(fhandle, n)
  end subroutine

  subroutine ms_copy_atoms(fhandle, copy, atoms, status)
    type(c_ptr), intent(in) :: fhandle
    integer(c_int), intent(in) :: copy
    integer(c_int), intent(out) :: atoms(:)
    integer(c_int), intent(out) :: status
    integer(c_int) :: n, i
    status = molsymm_copy_n_atoms(fhandle, copy - 1_c_int, n)
    if (status /= 0) return
    if (size(atoms) < n) then
      status = -1
      return
    end if
    status = molsymm_copy_atoms(fhandle, copy - 1_c_int, atoms)
    if (status == 0) then
      do i = 1, n
        atoms(i) = atoms(i) + 1_c_int
      end do
    end if
  end subroutine

  subroutine ms_fragment_mask(fhandle, mask, status)
    type(c_ptr), intent(in) :: fhandle
    integer(c_signed_char), intent(out) :: mask(:)
    integer(c_int), intent(out) :: status
    if (size(mask) <= 0) then
      status = -1
      return
    end if
    status = molsymm_fragment_mask(fhandle, mask)
  end subroutine

  subroutine ms_reconstruct(fhandle, rep_xyz, xyz_out, status)
    type(c_ptr), intent(in) :: fhandle
    real(c_double), intent(in) :: rep_xyz(:,:)   ! (3, n_rep)
    real(c_double), intent(out) :: xyz_out(:,:)  ! (3, n)
    integer(c_int), intent(out) :: status
    real(c_double), allocatable :: flat_in(:), flat_out(:)
    integer :: nrep, n, i
    integer(c_int) :: nrep_expect
    nrep = size(rep_xyz, 2)
    n = size(xyz_out, 2)
    status = molsymm_n_representatives(fhandle, nrep_expect)
    if (status /= 0) return
    if (nrep /= nrep_expect .or. n <= 0) then
      status = -1
      return
    end if
    allocate(flat_in(3*nrep), flat_out(3*n))
    do i = 1, nrep
      flat_in(3*(i-1)+1:3*(i-1)+3) = rep_xyz(:, i)
    end do
    status = molsymm_reconstruct(fhandle, flat_in, flat_out)
    if (status /= 0) return
    do i = 1, n
      xyz_out(:, i) = flat_out(3*(i-1)+1:3*(i-1)+3)
    end do
  end subroutine

  subroutine ms_subgroup_assembly(handle, subgroup, ahandle, status)
    type(c_ptr), intent(in) :: handle
    character(len=*), intent(in) :: subgroup
    type(c_ptr), intent(out) :: ahandle
    integer(c_int), intent(out) :: status
    character(kind=c_char), allocatable :: cname(:)
    integer :: i, n
    n = len_trim(subgroup)
    allocate(cname(n+1))
    do i = 1, n
      cname(i) = subgroup(i:i)
    end do
    cname(n+1) = c_null_char
    status = molsymm_subgroup_assembly(handle, cname, ahandle)
  end subroutine

  subroutine ms_release_assembly(ahandle)
    type(c_ptr), intent(inout) :: ahandle
    if (c_associated(ahandle)) call molsymm_free_assembly(ahandle)
    ahandle = c_null_ptr
  end subroutine

  subroutine ms_assembly_n_pieces(ahandle, n, status)
    type(c_ptr), intent(in) :: ahandle
    integer(c_int), intent(out) :: n, status
    status = molsymm_assembly_n_pieces(ahandle, n)
  end subroutine

  subroutine ms_assembly_piece_atoms(ahandle, piece, atoms, status)
    type(c_ptr), intent(in) :: ahandle
    integer(c_int), intent(in) :: piece
    integer(c_int), intent(out) :: atoms(:)
    integer(c_int), intent(out) :: status
    integer(c_int) :: n, i
    status = molsymm_assembly_piece_n_atoms(ahandle, piece - 1_c_int, n)
    if (status /= 0) return
    if (size(atoms) < n) then
      status = -1
      return
    end if
    status = molsymm_assembly_piece_atoms(ahandle, piece - 1_c_int, atoms)
    if (status == 0) then
      do i = 1, n
        atoms(i) = atoms(i) + 1_c_int
      end do
    end if
  end subroutine

  subroutine ms_assemble(ahandle, piece_xyz, xyz_out, status)
    type(c_ptr), intent(in) :: ahandle
    real(c_double), intent(in) :: piece_xyz(:,:)  ! (3, n_asym)
    real(c_double), intent(out) :: xyz_out(:,:)   ! (3, n)
    integer(c_int), intent(out) :: status
    real(c_double), allocatable :: flat_in(:), flat_out(:)
    integer :: nasym, n, i, p
    integer(c_int) :: npieces, n_piece, nasym_expect
    nasym = size(piece_xyz, 2)
    n = size(xyz_out, 2)
    status = molsymm_assembly_n_pieces(ahandle, npieces)
    if (status /= 0) return
    nasym_expect = 0_c_int
    do p = 0, npieces - 1
      status = molsymm_assembly_piece_n_atoms(ahandle, int(p, c_int), n_piece)
      if (status /= 0) return
      nasym_expect = nasym_expect + n_piece
    end do
    if (nasym /= nasym_expect .or. n <= 0) then
      status = -1
      return
    end if
    allocate(flat_in(3*nasym), flat_out(3*n))
    do i = 1, nasym
      flat_in(3*(i-1)+1:3*(i-1)+3) = piece_xyz(:, i)
    end do
    status = molsymm_assemble(ahandle, flat_in, flat_out)
    if (status /= 0) return
    do i = 1, n
      xyz_out(:, i) = flat_out(3*(i-1)+1:3*(i-1)+3)
    end do
  end subroutine

  subroutine ms_wyckoff_positions(point_group, thandle, status)
    character(len=*), intent(in) :: point_group
    type(c_ptr), intent(out) :: thandle
    integer(c_int), intent(out) :: status
    character(kind=c_char), allocatable :: cname(:)
    integer :: i, n
    n = len_trim(point_group)
    allocate(cname(n+1))
    do i = 1, n
      cname(i) = point_group(i:i)
    end do
    cname(n+1) = c_null_char
    status = molsymm_wyckoff_positions(cname, thandle)
  end subroutine

  subroutine ms_release_wyckoff_table(thandle)
    type(c_ptr), intent(inout) :: thandle
    call molsymm_free_wyckoff_table(thandle)
    thandle = c_null_ptr
  end subroutine

  subroutine ms_wyckoff_n_sites(thandle, n, status)
    type(c_ptr), intent(in) :: thandle
    integer(c_int), intent(out) :: n, status
    status = molsymm_wyckoff_n_sites(thandle, n)
  end subroutine

  subroutine ms_wyckoff_site_label(thandle, index, label, status)
    type(c_ptr), intent(in) :: thandle
    integer(c_int), intent(in) :: index
    character(len=*), intent(out) :: label
    integer(c_int), intent(out) :: status
    character(kind=c_char) :: buf(64)
    integer :: i
    status = molsymm_wyckoff_site_label(thandle, index - 1_c_int, buf, 64_c_int)
    label = ""
    if (status /= 0) return
    do i = 1, 64
      if (buf(i) == c_null_char) exit
      label(i:i) = buf(i)
    end do
  end subroutine

  subroutine ms_wyckoff_site_multiplicity(thandle, index, mult, status)
    type(c_ptr), intent(in) :: thandle
    integer(c_int), intent(in) :: index
    integer(c_int), intent(out) :: mult, status
    status = molsymm_wyckoff_site_multiplicity(thandle, index - 1_c_int, mult)
  end subroutine

  subroutine ms_wyckoff_site_dof(thandle, index, dof, status)
    type(c_ptr), intent(in) :: thandle
    integer(c_int), intent(in) :: index
    integer(c_int), intent(out) :: dof, status
    status = molsymm_wyckoff_site_dof(thandle, index - 1_c_int, dof)
  end subroutine

  subroutine ms_wyckoff_plans_nsites(point_group, n_sites, phandle, status)
    character(len=*), intent(in) :: point_group
    integer(c_int), intent(in) :: n_sites
    type(c_ptr), intent(out) :: phandle
    integer(c_int), intent(out) :: status
    character(kind=c_char), allocatable :: cname(:)
    integer :: i, n
    n = len_trim(point_group)
    allocate(cname(n+1))
    do i = 1, n
      cname(i) = point_group(i:i)
    end do
    cname(n+1) = c_null_char
    status = molsymm_wyckoff_plans_nsites(cname, n_sites, phandle)
  end subroutine

  subroutine ms_release_wyckoff_plans(phandle)
    type(c_ptr), intent(inout) :: phandle
    call molsymm_free_wyckoff_plans(phandle)
    phandle = c_null_ptr
  end subroutine

  subroutine ms_wyckoff_n_plans(phandle, n, status)
    type(c_ptr), intent(in) :: phandle
    integer(c_int), intent(out) :: n, status
    status = molsymm_wyckoff_n_plans(phandle, n)
  end subroutine

  subroutine ms_wyckoff_plan_n_orbits(phandle, plan, n, status)
    type(c_ptr), intent(in) :: phandle
    integer(c_int), intent(in) :: plan
    integer(c_int), intent(out) :: n, status
    status = molsymm_wyckoff_plan_n_orbits(phandle, plan - 1_c_int, n)
  end subroutine

  subroutine ms_wyckoff_plan_orbit_label(phandle, plan, orbit, label, status)
    type(c_ptr), intent(in) :: phandle
    integer(c_int), intent(in) :: plan, orbit
    character(len=*), intent(out) :: label
    integer(c_int), intent(out) :: status
    character(kind=c_char) :: buf(64)
    integer :: i
    status = molsymm_wyckoff_plan_orbit_label(phandle, plan - 1_c_int, orbit - 1_c_int, buf, 64_c_int)
    label = ""
    if (status /= 0) return
    do i = 1, 64
      if (buf(i) == c_null_char) exit
      label(i:i) = buf(i)
    end do
  end subroutine

  subroutine ms_assemble_from_au(point_group, z, au_xyz, labels, handle, status)
    character(len=*), intent(in) :: point_group
    integer(c_int), intent(in) :: z(:)
    real(c_double), intent(in) :: au_xyz(:,:)  ! (3, n_orbits)
    character(len=*), intent(in) :: labels(:)
    type(c_ptr), intent(out) :: handle
    integer(c_int), intent(out) :: status
    character(kind=c_char), allocatable :: cname(:)
    character(kind=c_char), allocatable, target :: clabels(:,:)
    type(c_ptr), allocatable :: label_ptrs(:)
    real(c_double), allocatable :: flat(:)
    integer :: i, j, n, norb, maxlen
    norb = size(z)
    if (size(au_xyz, 2) /= norb .or. size(labels) /= norb) then
      status = -1_c_int
      return
    end if
    n = len_trim(point_group)
    allocate(cname(n+1))
    do i = 1, n
      cname(i) = point_group(i:i)
    end do
    cname(n+1) = c_null_char
    maxlen = 0
    do i = 1, norb
      maxlen = max(maxlen, len_trim(labels(i)))
    end do
    allocate(clabels(maxlen+1, norb), label_ptrs(norb), flat(3*norb))
    do i = 1, norb
      clabels(:, i) = c_null_char
      do j = 1, len_trim(labels(i))
        clabels(j, i) = labels(i)(j:j)
      end do
      label_ptrs(i) = c_loc(clabels(1, i))
      flat(3*(i-1)+1:3*(i-1)+3) = au_xyz(:, i)
    end do
    status = molsymm_assemble_from_au(cname, z, flat, int(norb, c_int), label_ptrs, handle)
  end subroutine

end module molsymm
