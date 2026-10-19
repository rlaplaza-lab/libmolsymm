program test_fortran
  use molsymm
  use, intrinsic :: iso_c_binding
  implicit none
  integer(c_int) :: z(3), status, sigma, nd
  real(c_double) :: xyz(3,3)
  type(c_ptr) :: handle, dhandle
  character(len=32) :: name

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
  print *, "fortran smoke ok"
end program test_fortran
