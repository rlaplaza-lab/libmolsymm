program test_fortran
  use molsymm
  use, intrinsic :: iso_c_binding
  implicit none

  integer(c_int) :: z(3), status, sigma, order, nd, ios, nfail
  integer(c_int) :: expect_order, expect_sigma, got_order, got_sigma
  real(c_double) :: xyz(3,3), primary, max_rmsd, got_rmsd
  type(c_ptr) :: handle, dhandle
  character(len=32) :: name
  character(len=256) :: data_dir, manifest, xyz_path, rel, pg, source
  character(len=512) :: line
  type(ms_zatom) :: zatoms(3)
  real(c_double) :: ang
  integer :: u, nat
  integer(c_int), allocatable :: z_all(:)
  real(c_double), allocatable :: xyz_all(:,:)
  character(len=8) :: sym
  real(c_double) :: xx, yy, zz

  nfail = 0

  ! --- water smoke ---
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
  call ms_order(handle, order, status)
  print *, "point group ", trim(name), " sigma ", sigma, " order ", order
  if (trim(name) /= "C2v" .or. sigma /= 2 .or. order /= 4) then
    nfail = nfail + 1
  end if

  call ms_unique_displacements(handle, 1.0d-3, dhandle, status)
  call ms_n_displacements(dhandle, nd, status)
  print *, "n displacements ", nd
  if (nd /= 6) nfail = nfail + 1

  call ms_release_displacements(dhandle)

  block
    type(c_ptr) :: fhandle, ahandle
    integer(c_int) :: nrep, npieces, atom0
    real(c_double) :: rep_xyz(3,2), xyz_out(3,3)
    call ms_fragmentation(handle, fhandle, status)
    if (status /= 0) then
      print *, "fragmentation failed: ", ms_last_error()
      nfail = nfail + 1
    else
      call ms_n_representatives(fhandle, nrep, status)
      if (nrep /= 2) nfail = nfail + 1
      call ms_representative_atom(fhandle, 1_c_int, atom0, status)
      if (status /= 0 .or. atom0 /= 1_c_int) then
        print *, "representative atom 1-based index expected 1, got ", atom0
        nfail = nfail + 1
      end if
      rep_xyz = xyz(:, 1:2)
      call ms_reconstruct(fhandle, rep_xyz, xyz_out, status)
      if (status /= 0 .or. maxval(abs(xyz_out - xyz)) > 1.0d-6) nfail = nfail + 1
      call ms_release_fragmentation(fhandle)
    end if
    call ms_subgroup_assembly(handle, "C2v", ahandle, status)
    if (status /= 0) then
      print *, "subgroup_assembly failed: ", ms_last_error()
      nfail = nfail + 1
    else
      call ms_assembly_n_pieces(ahandle, npieces, status)
      if (npieces /= 1) nfail = nfail + 1
      call ms_assemble(ahandle, xyz(:, 1:2), xyz_out, status)
      if (status /= 0 .or. maxval(abs(xyz_out - xyz)) > 1.0d-6) nfail = nfail + 1
      call ms_release_assembly(ahandle)
    end if
  end block

  call ms_release(handle)

  ! --- z-matrix smoke ---
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
  if (trim(name) /= "C2v") nfail = nfail + 1
  call ms_release(handle)

  ! --- smoke manifest ---
  call get_environment_variable("MOLSYMM_TEST_DATA_DIR", data_dir, status=ios)
  if (len_trim(data_dir) == 0) then
    data_dir = "tests/data"
  end if
  manifest = trim(data_dir)//"/manifest_smoke.tsv"
  open(newunit=u, file=trim(manifest), status="old", action="read", iostat=ios)
  if (ios /= 0) then
    print *, "cannot open manifest ", trim(manifest)
    stop 1
  end if
  read(u, "(A)", iostat=ios) line  ! header
  do
    read(u, "(A)", iostat=ios) line
    if (ios /= 0) exit
    if (len_trim(line) == 0) cycle
    call parse_manifest_line(line, rel, pg, expect_order, expect_sigma, primary, max_rmsd, source)
    xyz_path = trim(data_dir)//"/"//trim(rel)
    call load_xyz_file(xyz_path, z_all, xyz_all, nat, ios)
    if (ios /= 0) then
      print *, "FAIL load ", trim(xyz_path)
      nfail = nfail + 1
      cycle
    end if
    call ms_analyze(z_all, xyz_all, handle, status, primary=primary)
    if (status /= 0) then
      print *, "FAIL analyze ", trim(rel), " ", ms_last_error()
      nfail = nfail + 1
      if (allocated(z_all)) deallocate(z_all)
      if (allocated(xyz_all)) deallocate(xyz_all)
      cycle
    end if
    call ms_point_group_name(handle, name, status)
    call ms_symmetry_number(handle, got_sigma, status)
    call ms_order(handle, got_order, status)
    call ms_rmsd(handle, got_rmsd, status)
    if (trim(name) /= trim(pg) .or. got_order /= expect_order .or. got_sigma /= expect_sigma) then
      print *, "FAIL ", trim(rel), " expect ", trim(pg), expect_order, expect_sigma, &
               " got ", trim(name), got_order, got_sigma
      nfail = nfail + 1
    else if (got_rmsd >= max_rmsd) then
      print *, "FAIL ", trim(rel), " rmsd ", got_rmsd, " max ", max_rmsd
      nfail = nfail + 1
    end if
    call ms_release(handle)
    if (allocated(z_all)) deallocate(z_all)
    if (allocated(xyz_all)) deallocate(xyz_all)
  end do
  close(u)

  if (nfail /= 0) then
    print *, nfail, " fortran failures"
    stop 1
  end if
  print *, "fortran smoke ok"
contains

  subroutine parse_manifest_line(line, rel, pg, order, sigma, primary, max_rmsd, source)
    character(len=*), intent(in) :: line
    character(len=*), intent(out) :: rel, pg, source
    integer(c_int), intent(out) :: order, sigma
    real(c_double), intent(out) :: primary, max_rmsd
    integer :: p1, p2, p3, p4, p5, p6
    character(len=64) :: tok
    p1 = index(line, char(9))
    rel = line(1:p1-1)
    p2 = p1 + index(line(p1+1:), char(9))
    pg = adjustl(line(p1+1:p2-1))
    p3 = p2 + index(line(p2+1:), char(9))
    tok = adjustl(line(p2+1:p3-1)); read(tok, *) order
    p4 = p3 + index(line(p3+1:), char(9))
    tok = adjustl(line(p3+1:p4-1)); read(tok, *) sigma
    p5 = p4 + index(line(p4+1:), char(9))
    tok = adjustl(line(p4+1:p5-1)); read(tok, *) primary
    p6 = p5 + index(line(p5+1:), char(9))
    tok = adjustl(line(p5+1:p6-1)); read(tok, *) max_rmsd
    source = adjustl(line(p6+1:))
  end subroutine

  subroutine load_xyz_file(path, z_all, xyz_all, nat, ios)
    character(len=*), intent(in) :: path
    integer(c_int), allocatable, intent(out) :: z_all(:)
    real(c_double), allocatable, intent(out) :: xyz_all(:,:)
    integer, intent(out) :: nat, ios
    integer :: uu, k
    character(len=256) :: buf
    open(newunit=uu, file=trim(path), status="old", action="read", iostat=ios)
    if (ios /= 0) return
    read(uu, *, iostat=ios) nat
    if (ios /= 0) then; close(uu); return; end if
    read(uu, "(A)", iostat=ios) buf
    allocate(z_all(nat), xyz_all(3, nat))
    do k = 1, nat
      read(uu, *, iostat=ios) sym, xx, yy, zz
      if (ios /= 0) then; close(uu); return; end if
      z_all(k) = int(atomic_number(trim(sym)), c_int)
      xyz_all(:, k) = [xx, yy, zz]
    end do
    close(uu)
  end subroutine

  integer function atomic_number(s)
    character(len=*), intent(in) :: s
    select case (s)
    case ("H"); atomic_number = 1
    case ("He"); atomic_number = 2
    case ("Li"); atomic_number = 3
    case ("Be"); atomic_number = 4
    case ("B"); atomic_number = 5
    case ("C"); atomic_number = 6
    case ("N"); atomic_number = 7
    case ("O"); atomic_number = 8
    case ("F"); atomic_number = 9
    case ("Ne"); atomic_number = 10
    case ("Na"); atomic_number = 11
    case ("Mg"); atomic_number = 12
    case ("Al"); atomic_number = 13
    case ("Si"); atomic_number = 14
    case ("P"); atomic_number = 15
    case ("S"); atomic_number = 16
    case ("Cl"); atomic_number = 17
    case ("Ar"); atomic_number = 18
    case ("K"); atomic_number = 19
    case ("Ca"); atomic_number = 20
    case ("Sc"); atomic_number = 21
    case ("Ti"); atomic_number = 22
    case ("V"); atomic_number = 23
    case ("Cr"); atomic_number = 24
    case ("Mn"); atomic_number = 25
    case ("Fe"); atomic_number = 26
    case ("Co"); atomic_number = 27
    case ("Ni"); atomic_number = 28
    case ("Cu"); atomic_number = 29
    case ("Zn"); atomic_number = 30
    case ("Ga"); atomic_number = 31
    case ("Ge"); atomic_number = 32
    case ("As"); atomic_number = 33
    case ("Se"); atomic_number = 34
    case ("Br"); atomic_number = 35
    case ("Kr"); atomic_number = 36
    case ("Rb"); atomic_number = 37
    case ("Sr"); atomic_number = 38
    case ("Y"); atomic_number = 39
    case ("Zr"); atomic_number = 40
    case ("Nb"); atomic_number = 41
    case ("Mo"); atomic_number = 42
    case ("Tc"); atomic_number = 43
    case ("Ru"); atomic_number = 44
    case ("Rh"); atomic_number = 45
    case ("Pd"); atomic_number = 46
    case ("Ag"); atomic_number = 47
    case ("Cd"); atomic_number = 48
    case ("In"); atomic_number = 49
    case ("Sn"); atomic_number = 50
    case ("Sb"); atomic_number = 51
    case ("Te"); atomic_number = 52
    case ("I"); atomic_number = 53
    case ("Xe"); atomic_number = 54
    case ("Cs"); atomic_number = 55
    case ("Ba"); atomic_number = 56
    case ("La"); atomic_number = 57
    case ("Hf"); atomic_number = 72
    case ("Ta"); atomic_number = 73
    case ("W"); atomic_number = 74
    case ("Re"); atomic_number = 75
    case ("Os"); atomic_number = 76
    case ("Ir"); atomic_number = 77
    case ("Pt"); atomic_number = 78
    case ("Au"); atomic_number = 79
    case ("Hg"); atomic_number = 80
    case default
      print *, "unknown element ", s
      stop 1
    end select
  end function

end program test_fortran
