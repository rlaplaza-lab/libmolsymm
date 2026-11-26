module molsymm
  use, intrinsic :: iso_c_binding
  implicit none
  private
  public :: ms_analyze, ms_release, ms_point_group_name, ms_symmetry_number, ms_order
  public :: ms_symmetrize, ms_unique_displacements, ms_n_displacements
  public :: ms_displacement_geometry, ms_hessian_from_gradients, ms_release_displacements
  public :: ms_last_error, ms_zatom, ms_zmatrix_to_xyz, ms_analyze_zmatrix
  public :: ms_fragmentation, ms_release_fragmentation, ms_n_representatives
  public :: ms_representative_atom, ms_n_copies, ms_copy_atoms, ms_fragment_mask
  public :: ms_reconstruct, ms_subgroup_assembly, ms_release_assembly
  public :: ms_assembly_n_pieces, ms_assembly_piece_atoms, ms_assemble

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

  function molsymm_symmetrize(a, rmsd, xyz_out) bind(c, name="molsymm_symmetrize")
    import :: c_int, c_ptr, c_double
    type(c_ptr), value :: a
    real(c_double), intent(out) :: rmsd
    real(c_double), intent(out) :: xyz_out(*)
    integer(c_int) :: molsymm_symmetrize
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
    status = molsymm_representative_atom(fhandle, index, atom)
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
    integer(c_int) :: n
    status = molsymm_copy_n_atoms(fhandle, copy, n)
    if (status /= 0) return
    if (size(atoms) < n) then
      status = -1
      return
    end if
    status = molsymm_copy_atoms(fhandle, copy, atoms)
  end subroutine

  subroutine ms_fragment_mask(fhandle, mask, status)
    type(c_ptr), intent(in) :: fhandle
    integer(c_signed_char), intent(out) :: mask(:)
    integer(c_int), intent(out) :: status
    status = molsymm_fragment_mask(fhandle, mask)
  end subroutine

  subroutine ms_reconstruct(fhandle, rep_xyz, xyz_out, status)
    type(c_ptr), intent(in) :: fhandle
    real(c_double), intent(in) :: rep_xyz(:,:)   ! (3, n_rep)
    real(c_double), intent(out) :: xyz_out(:,:)  ! (3, n)
    integer(c_int), intent(out) :: status
    real(c_double), allocatable :: flat_in(:), flat_out(:)
    integer :: nrep, n, i
    nrep = size(rep_xyz, 2)
    n = size(xyz_out, 2)
    allocate(flat_in(3*nrep), flat_out(3*n))
    do i = 1, nrep
      flat_in(3*(i-1)+1:3*(i-1)+3) = rep_xyz(:, i)
    end do
    status = molsymm_reconstruct(fhandle, flat_in, flat_out)
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
    integer(c_int) :: n
    status = molsymm_assembly_piece_n_atoms(ahandle, piece, n)
    if (status /= 0) return
    if (size(atoms) < n) then
      status = -1
      return
    end if
    status = molsymm_assembly_piece_atoms(ahandle, piece, atoms)
  end subroutine

  subroutine ms_assemble(ahandle, piece_xyz, xyz_out, status)
    type(c_ptr), intent(in) :: ahandle
    real(c_double), intent(in) :: piece_xyz(:,:)  ! (3, n_asym)
    real(c_double), intent(out) :: xyz_out(:,:)   ! (3, n)
    integer(c_int), intent(out) :: status
    real(c_double), allocatable :: flat_in(:), flat_out(:)
    integer :: nasym, n, i
    nasym = size(piece_xyz, 2)
    n = size(xyz_out, 2)
    allocate(flat_in(3*nasym), flat_out(3*n))
    do i = 1, nasym
      flat_in(3*(i-1)+1:3*(i-1)+3) = piece_xyz(:, i)
    end do
    status = molsymm_assemble(ahandle, flat_in, flat_out)
    do i = 1, n
      xyz_out(:, i) = flat_out(3*(i-1)+1:3*(i-1)+3)
    end do
  end subroutine

end module molsymm
