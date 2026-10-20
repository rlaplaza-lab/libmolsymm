program test_fortran
  use molsymm
  use, intrinsic :: iso_c_binding
  implicit none
  integer(c_int) :: z(3), status, sigma, nd
  real(c_double) :: xyz(3,3)
  type(c_ptr) :: handle, dhandle
  character(len=32) :: name
  type(ms_zatom) :: zatoms(3)
  real(c_double) :: ang

  z = [8_c_int, 1_c_int, 1_c_int]
  xyz(:,1) = [0.0d0, 0.0d0, 0.1173d0]
  xyz(:,2) = [0.0d0, 0.7572d0, -0.4692d0]
  xyz(:,3) = [0.0d0, -0.7572d0, -0.4692d0]

  call ms_analyze(z, xyz, handle, status)
  if (status /= 0) then
    print *, "analyze failed: ", ms_last_error()
    stop 1
  end if
  call ms_point_group_name(handle, name, status)
  call ms_symmetry_number(handle, sigma, status)
  print *, "point group ", trim(name), " sigma ", sigma
  if (trim(name) /= "C2v" .or. sigma /= 2) stop 1

  call ms_unique_displacements(handle, 1.0d-3, dhandle, status)
  call ms_n_displacements(dhandle, nd, status)
  print *, "n displacements ", nd
  if (nd /= 6) stop 1

  call ms_release_displacements(dhandle)
  call ms_release(handle)

  ang = 104.5d0 * acos(-1.0d0) / 180.0d0
  zatoms = ms_zatom(0, 0, 0, 0, 0.0d0, 0.0d0, 0.0d0)
  zatoms(1)%atomic_number = 8
  zatoms(2)%atomic_number = 1
  zatoms(2)%bond_to = 1
  zatoms(2)%distance = 0.957d0
  zatoms(3)%atomic_number = 1
  zatoms(3)%bond_to = 1
  zatoms(3)%angle_to = 2
  zatoms(3)%distance = 0.957d0
  zatoms(3)%angle = ang
  call ms_analyze_zmatrix(zatoms, handle, status)
  if (status /= 0) then
    print *, "z-matrix analyze failed: ", ms_last_error()
    stop 1
  end if
  call ms_point_group_name(handle, name, status)
  print *, "z-matrix point group ", trim(name)
  if (trim(name) /= "C2v") stop 1
  call ms_release(handle)
  print *, "fortran smoke ok"
end program test_fortran
