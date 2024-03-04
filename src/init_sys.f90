module init_sys_mod

  implicit none
contains
  subroutine init_sys()
#include <petsc/finclude/petsc.h>
    use slepceps
    use petscmat
    use petsc_mod
    use petsc_wrapper
    use kinds
    use globals
    use readsys
    use lattice_mod
    use ft_mod
    use conquest_mod
    use kmat_mod
    use slepc_mod
    use init_ep_coupling
    use phonon_mod
    use loe_mod
    use misc
    use error_handler
    use k_on_demand_mod
    use dftu_mod
    use fix_mat_mod

    implicit none

    integer :: ierr, ikx, iky, imu1, imu2, imin, n, i1, i2, i3, iter, j1, j2
    integer, allocatable :: tomat_tmp(:, :), inao_tmp(:), neigh_tmp(:), ndim_tmp(:), atoms_tmp(:), nlist_tmp(:),&
    &nzrow_c(:), nzrow_l(:), nzrow_r(:)
    PetscScalar :: trace_c, trace_elec, fix_l, fix_r
    integer :: iunitscf, iunith1, iunith2, iunith3, ndfts, imu, jmu
    logical :: makeks
    character(strln) :: sdummy
    PetscErrorCode :: p_ierr

! petsc stuff
    real(dp) :: info(MAT_INFO_SIZE)
    PetscInt :: nl1, nl2, ii(1), jj(1), ione, nev, i, j, nlc1, nlc2, m, ic(2), ir(2), k
    PetscViewer :: viewer, xview
    PetscScalar :: zout(1), p_scal
    PetscReal :: norm, norm_b, norm_a, norm_m, d1, r0, rm, rp, deps,&
   &norm_r0, norm_rp, norm_rm, dd, d2, dm, ds1, ds2

    Mat :: p_tmp, p_helec, p_hecc, p_selec, p_tmp_read, p_mat_tmp
    Mat :: p_tmp1, p_tmp2, p_tmp3, p_tmp4, p_w, p_w2, p_ev
    Vec :: p_ew1, p_ew2, p_v1, p_v2
    integer :: irow(2), icol(1), ncols, idxm(1), idym(1), iat, ixyz, ik, nnz, diag_dim_save
    integer, allocatable :: cols(:), idx(:), idy(:)
    integer, pointer :: iidx(:), iidy(:)
    PetscScalar, allocatable :: p_vals(:), p_vals2(:)
    PetscScalar, allocatable, target :: t_x(:)
    PetscScalar, pointer :: p_x(:), p_y(:, :)
    PetscBool :: ldone
    character(5) :: titel
    character(strln) :: s1, s2, s3, what, tmpfile, atstr, xyzstr
    MatType ::  mattype
    MatFactorInfo :: matfacinfo(MAT_FACTORINFO_SIZE)
    PetscInt, pointer :: i_row_p(:), i_col_p(:)
    integer(8) :: counti, count_rate, countf
    logical :: l_init_ep_mat, l_read_dh 
    Mat, pointer :: p_tmp_k(:), p_tmp_cc(:, :)

! debug ep and loe
!~       l_ep_loe=.false.
!~       l_ep=.false.
!-------------------

! petsc stuff -----------------

    l_init_ep_mat = .true.

    write (pstr_out, fmt='(2A)') "C: ", trim(ecc_dir); call petsc_print_master()
    call read_sysinfo(xyz_ecc, species_ecc, imu_ecc, nat_ecc, &
   &   dlat_c, ncell_c, nmat_c, ef_c, n_electrons_ecc, nspecies_c, species_c, &
   &   trim(ecc_dir)//"/sys_info.dat")
!~ nmat_c read from sys_info.dat is obosolte. It should contain the
!~ total number of nonzeros over all the realspace unitcells of the realspace
!~ matrices (H, K, S). Calculate this later directly from K mat.

    if (l_use_sigma_l) then
      write (pstr_out, fmt='(2A)') "L electrode: ", trim(elec_l_dir)
      call petsc_print_master()
      call read_sysinfo(xyz_elec_l, species_elec_l, imu_elec_l, nat_elec_l, &
     &  dlat_l, ncell_l, nmat_l, ef_l, n_electrons_l, nspecies_l, species_l, &
     &  trim(elec_l_dir)//"/sys_info.dat")
    else
      l_no_noneq = .true.
      nat_elec_l = 0
      nmat_l = 0
      ncell_l = 0
      allocate(xyz_elec_l(1,1))
      write (pstr_out, fmt='(A)') "no left electrode -> don't add self energy"
      call petsc_print_master()
    end if
    
    if (l_use_sigma_r) then  
      write (pstr_out, fmt='(2A)') "R electrode: ", trim(elec_r_dir)
      call petsc_print_master()
      call read_sysinfo(xyz_elec_r, species_elec_r, imu_elec_r, nat_elec_r, &
     &  dlat_r, ncell_r, nmat_r, ef_r, n_electrons_r, nspecies_r, species_r, &
     &  trim(elec_r_dir)//"/sys_info.dat")
    else
      l_no_noneq = .true.
      nat_elec_r = 0
      nmat_r = 0
      ncell_r = 0
      allocate(xyz_elec_r(1,1))
      write (pstr_out, fmt='(A)') "no right electrode -> don't add self energy"
      call petsc_print_master()
    end if

    if (l_use_sigma_l.and.(.not.l_use_sigma_r)) then 
      ef_r = ef_l
    else if (l_use_sigma_r.and.(.not.l_use_sigma_l)) then 
      ef_l = ef_r
    else if ((.not.l_use_sigma_r).and.(.not.l_use_sigma_l)) then
      ef_l = ef_c_fix
      ef_r = ef_c_fix
      ef_c = ef_c_fix
!~       write (errormsg, fmt='(A,i8)') "no L and R electrode present ", 666
!~       call error()
    end if
  

    diag_dim_save = diag_dim
!~    diag_dim = 3
!~    if ((k_mat_mode .eq. 2)) then
!~      if (diag_dim .eq. 1) then
!~        ncell_c(2:3) = 0
!~        ncell_l(2:3) = 0
!~        ncell_r(2:3) = 0
!~      else if (diag_dim .eq. 2) then
!~        ncell_c(3) = 0
!~        ncell_l(3) = 0
!~        ncell_r(3) = 0
!~      end if
!~    end if

    allocate (lcr_info(nat_ecc), stat=ierr)
    if (ierr .ne. 0) then
      write (errormsg, fmt='(A,i8)') "allocation error init_sys ", ierr
      call error()
    end if

    if (.not.l_diag_skip_check) then
      call check_sys(xyz_ecc, xyz_elec_l, xyz_elec_r, nat_ecc, nat_elec_l, &
        nat_elec_r, lcr_info)
      write (pstr_out, fmt='(A)') "electrodes ok"; call petsc_print_master()
    end if

    nmu_c = sum(imu_ecc(nat_elec_l + 1:nat_ecc - nat_elec_r))
    nmu_l = sum(imu_ecc(1:nat_elec_l))
    nmu_r = sum(imu_ecc(nat_ecc - nat_elec_r + 1:nat_ecc))

    nmu_ecc = sum(imu_ecc)

    nmu_c = nmu_ecc

! mapping imu to atoms, we need this for ep coupling later
    allocate (imu_to_at(nmu_c), stat=ierr)
    if (ierr .ne. 0) then
      write (errormsg, fmt='(A,i8)') "allocation error imu_to_at ", ierr
      call error()
    end if
    k = 0
    do i = 1, nat_ecc
      do j = 1, imu_ecc(i)
        k = k + 1
        imu_to_at(k) = i
      end do
    end do

    makeks = .not. kpoints_from_file
    write (pstr_out, fmt='(A)') "ecc"; call petsc_print_master()
    call get_rlat(dlat_c, rlat_c, kp_c, wkp_c, nk_c, nkx_dez, nky_dez, nkz_diag, makeks)
    
    if (l_use_sigma_l) then
      write (pstr_out, fmt='(A)') "left electrode"; call petsc_print_master()
      call get_rlat(dlat_l, rlat_l, kp_l, wkp_l, nk_l, nkx_dez, nky_dez, nkz_diag, makeks)
    end if
    if (l_use_sigma_r) then
      write (pstr_out, fmt='(A)') "right electrode"; call petsc_print_master()
      call get_rlat(dlat_r, rlat_r, kp_r, wkp_r, nk_r, nkx_dez, nky_dez, nkz_diag, makeks)
    end if
    
    if (kpoints_from_file) then
      call read_kpoint_file(kp_l, wkp_l, nk_l)
      call read_kpoint_file(kp_r, wkp_r, nk_r)
    end if

! it seems the supercells is usually to small in x-y dimension so we have to sum over k-points.
! for the time being assume them to be equal in l, c, and r
    nk_l = nk_c
    nk_r = nk_c
    wkp_l = wkp_c
    wkp_r = wkp_c
    nktot = sum(wkp_c)

    write (pstr_out, fmt='(A,e24.12)') "ef_c=", ef_c; call petsc_print_master()
    write (pstr_out, fmt='(A,e24.12)') "ef_l=", ef_l; call petsc_print_master()
    write (pstr_out, fmt='(A,e24.12)') "ef_r=", ef_r; call petsc_print_master()
    write (pstr_out, fmt='(A,e24.12)') "ef_c-ef_l=", ef_c - ef_l; call petsc_print_master()
    write (pstr_out, fmt='(A,e24.12)') "ef_c-ef_r=", ef_c - ef_r; call petsc_print_master()
    write (pstr_out, fmt='(A,i8)') "nmu ecc ", nmu_ecc; call petsc_print_master()
    write (pstr_out, fmt='(A,i8)') "nmu c ", nmu_c; call petsc_print_master()
    write (pstr_out, fmt='(A,i8)') "nmu l ", nmu_l; call petsc_print_master()
    write (pstr_out, fmt='(A,i8)') "nmu r ", nmu_r; call petsc_print_master()
    write (pstr_out, fmt='(A,3i8)') "ncell_c ", ncell_c(1:3); call petsc_print_master()
    write (pstr_out, fmt='(A,2e24.12)') "number of electrons ECC ", n_electrons_ecc; call petsc_print_master()

! petsc: init arrays
    allocate (p_h00_cc(-ncell_c(1):ncell_c(1), -ncell_c(2):ncell_c(2)), &
              p_s00_cc(-ncell_c(1):ncell_c(1), -ncell_c(2):ncell_c(2)), &
              p_dmatxy(-ncell_c(1):ncell_c(1), -ncell_c(2):ncell_c(2), -ncell_c(3):ncell_c(3)), &
              p_dmatxy_l(-ncell_c(1):ncell_c(1), -ncell_c(2):ncell_c(2)), &
              p_dmatxy_r(-ncell_c(1):ncell_c(1), -ncell_c(2):ncell_c(2)), &
              p_h00_l(-ncell_c(1):ncell_c(1), -ncell_c(2):ncell_c(2)), &
              p_s00_l(-ncell_c(1):ncell_c(1), -ncell_c(2):ncell_c(2)), &
              p_k00_l(-ncell_c(1):ncell_c(1), -ncell_c(2):ncell_c(2)), &
              p_h01_l(-ncell_c(1):ncell_c(1), -ncell_c(2):ncell_c(2)), &
              p_s01_l(-ncell_c(1):ncell_c(1), -ncell_c(2):ncell_c(2)), &
              p_k01_l(-ncell_c(1):ncell_c(1), -ncell_c(2):ncell_c(2)), &
              p_h10_l(-ncell_c(1):ncell_c(1), -ncell_c(2):ncell_c(2)), &
              p_s10_l(-ncell_c(1):ncell_c(1), -ncell_c(2):ncell_c(2)), &
              p_k10_l(-ncell_c(1):ncell_c(1), -ncell_c(2):ncell_c(2)), &
              p_h00_r(-ncell_c(1):ncell_c(1), -ncell_c(2):ncell_c(2)), &
              p_s00_r(-ncell_c(1):ncell_c(1), -ncell_c(2):ncell_c(2)), &
              p_k00_r(-ncell_c(1):ncell_c(1), -ncell_c(2):ncell_c(2)), &
              p_h01_r(-ncell_c(1):ncell_c(1), -ncell_c(2):ncell_c(2)), &
              p_s01_r(-ncell_c(1):ncell_c(1), -ncell_c(2):ncell_c(2)), &
              p_k01_r(-ncell_c(1):ncell_c(1), -ncell_c(2):ncell_c(2)), &
              p_h10_r(-ncell_c(1):ncell_c(1), -ncell_c(2):ncell_c(2)), &
              p_s10_r(-ncell_c(1):ncell_c(1), -ncell_c(2):ncell_c(2)), &
              p_k10_r(-ncell_c(1):ncell_c(1), -ncell_c(2):ncell_c(2)), &              
              stat=ierr)
    if (ierr .ne. 0) then
      write (errormsg, fmt='(A,i8)') "allocation error init_sys ", ierr
      call error()
    end if

    p_dmatxy = PETSC_NULL_MAT


    if (calc_current_density) then
      allocate(p_vnl00_cc(-ncell_c(1):ncell_c(1), -ncell_c(2):ncell_c(2)), stat = ierr)
      if (ierr .ne. 0) then
        write (errormsg, fmt='(A,i8)') "allocation error init_sys ", ierr
        call error()
      end if
    end if

    if ((k_mat_mode .eq. 2) .or. (l_dftu)) then
      allocate(p_h01_cc(-ncell_c(1):ncell_c(1), -ncell_c(2):ncell_c(2)), &
        p_s01_cc(-ncell_c(1):ncell_c(1), -ncell_c(2):ncell_c(2)), &
        p_h10_cc(-ncell_c(1):ncell_c(1), -ncell_c(2):ncell_c(2)), &
        p_s10_cc(-ncell_c(1):ncell_c(1), -ncell_c(2):ncell_c(2)), &
        p_h00_diag(-ncell_c(1):ncell_c(1), -ncell_c(2):ncell_c(2),-ncell_c(3):ncell_c(3)), &
        p_s00_diag(-ncell_c(1):ncell_c(1), -ncell_c(2):ncell_c(2),-ncell_c(3):ncell_c(3)), &
        stat = ierr)
      if (ierr .ne. 0) then
        write (errormsg, fmt='(A,i8)') "allocation error init_sys ", ierr
        call error()
      end if
    end if

    if (l_dftu) then
      allocate(p_k00_diag(-ncell_c(1):ncell_c(1), -ncell_c(2):ncell_c(2),-ncell_c(3):ncell_c(3)), &
     &  stat = ierr)
      if (ierr .ne. 0) then
        write (errormsg, fmt='(A,i8)') "allocation error init_sys ", ierr
        call error()
      end if
    end if

! additional matrices for ep coupling
    if (l_ep) then
      n_ep_active_atoms = ep_active_atoms(2) - ep_active_atoms(1) + 1
      allocate (p_dh_cc(-ncell_c(1):ncell_c(1), -ncell_c(2):ncell_c(2), ep_active_atoms(1):ep_active_atoms(2), 3), &
                stat=ierr)
    end if
! petsc: init arrays -----------------

    call mpi_barrier(MPI_COMM_WORLD, ierr)

    if (l_ep) then
      nat3 = nat_ecc*3
      call init_phonons()
      allocate (p_dhk_cc(nk_c, ep_active_atoms(1):ep_active_atoms(2), 3), & ! i think would be better to exchange indices 2 and 3
                p_Kphonon00k_cc(nk_c), p_phonon_EV(nk_c), p_phonon_EW(nk_c), &
                n_ep_modes_k(nk_c), &
                stat=ierr)
      if (ierr .ne. 0) then
        write (errormsg, fmt='(A,i8)') "allocation failure init_sys ", ierr
        call error()
      end if
      call get_phonones()
    end if

!work around to extend system, remove this at some point
! i think this is depreciated and can be removed
!~       nmu_l=0
!~       nmu_r=0

! for now we assume nonzero structure following from Conquest partion, i.e.,
! the sparsity is given by the Conquest cutoff and explicit zeros are kept.
! While this increases the memory usage somewhat it allows us to use the same
! nonzero structure for all matrices which should speed up the matrix operations
! in PETSc and we can use MatDuplicate with MAT_SHARE_NONZERO_PATTERN

! Scattering region
    call system_clock(counti, count_rate)
    pstr_out = "read data of scattering region..."; call petsc_print_master()
    allocate (nrow_pproc(2, nprocs), stat=ierr)
    if (ierr .ne. 0) then
      write (pstr_out, fmt='(A,i8)') "allocation error  init_sys ", ierr; call petsc_print_master()
    end if

    ! H-matrix

    s1 = "0"
    s2 = "0"
    tmpfile = trim(ecc_dir)//"/K_"//trim(spin_str)//"_"//trim(s1)//"_"//trim(s2)//"_0"//"_petsc.dat"
    call petsc_mat_load(p_dmatxy(0, 0, 0), tmpfile, .false., .true., nrow_pproc, mattype_sparse, nzrow_out=nzrow_c)
    nzmax_row = maxval(nzrow_c)
    tmpfile = trim(elec_l_dir)//"/K_"//trim(spin_str)//"_"//trim(s1)//"_"//trim(s2)//"_0"//"_petsc.dat"
    if (l_use_sigma_l) call petsc_mat_load(p_dmatxy(0, 0, 0), tmpfile, .false., .true., nrow_pproc, mattype_sparse, nzrow_out=nzrow_l)
    tmpfile = trim(elec_r_dir)//"/K_"//trim(spin_str)//"_"//trim(s1)//"_"//trim(s2)//"_0"//"_petsc.dat"
    if (l_use_sigma_r) call petsc_mat_load(p_dmatxy(0, 0, 0), tmpfile, .false., .true., nrow_pproc, mattype_sparse, nzrow_out=nzrow_r)

!~     call nzs_to_procs(nrow_pproc, nzrow_c)
    call rows_to_procs(nrow_pproc, nmu_ecc) ! use rows to proc in accordance with pexsi

    do i1=1,nprocs
      write(pstr_out,fmt='(A,5i8)') "nz ",i1,nrow_pproc(1,i1),nrow_pproc(2,i1),&
      &nrow_pproc(2,i1)-nrow_pproc(1,i1)+1, nmu_ecc ; call petsc_print_master()
    end do

    nmat_c = 0

    do i1 = -ncell_c(1), ncell_c(1)
      do i2 = -ncell_c(2), ncell_c(2)

        s1 = int2str(i1)
        s2 = int2str(i2)

        ! K-matrix ! Loading the density matrix, which will be overwritten during SCF
        tmpfile = trim(ecc_dir)//"/K_"//trim(spin_str)//"_"//trim(s1)//"_"//trim(s2)//"_0"//"_petsc.dat"
        call petsc_mat_load(p_dmatxy(i1, i2, 0), tmpfile, .true., .false., nrow_pproc, &
          mattype_sparse, cols_loc, nzrow_loc, nnz_out = nnz)
        nmat_c = nmat_c + nnz

        ! H-matrix
        tmpfile = trim(ecc_dir)//"/H_"//trim(spin_str)//"_"//trim(s1)//"_"//trim(s2)//"_0"//"_petsc.dat"
        call MatDuplicate(p_dmatxy(i1, i2, 0), MAT_SHARE_NONZERO_PATTERN, p_h00_cc(i1, i2), ierr)
        call petsc_mat_load(p_h00_cc(i1, i2), tmpfile, .false., .false., nrow_pproc, mattype_sparse)
        
        if (ncell_c(3) .ne. 0) then
          ! K-matrix ! Loading the density matrix
          tmpfile = trim(ecc_dir)//"/K_"//trim(spin_str)//"_"//trim(s1)//"_"//trim(s2)//"_-1"//"_petsc.dat"
          call petsc_mat_load(p_dmatxy(i1, i2, -1), tmpfile, .true., .false., nrow_pproc, &
          mattype_sparse, cols_loc, nzrow_loc, nnz_out = nnz)
          nmat_c = nmat_c + nnz
  
  
          ! K-matrix ! Loading the density matrix
          tmpfile = trim(ecc_dir)//"/K_"//trim(spin_str)//"_"//trim(s1)//"_"//trim(s2)//"_1"//"_petsc.dat"
          call petsc_mat_load(p_dmatxy(i1, i2, 1), tmpfile, .true., .false., nrow_pproc, &
          mattype_sparse, cols_loc, nzrow_loc, nnz_out = nnz)
          nmat_c = nmat_c + nnz
        end if

        if ((k_mat_mode .eq. 2) .or. l_dftu) then
          if (diag_dim.eq.3) then
            tmpfile = trim(ecc_dir)//"/H_"//trim(spin_str)//"_"//trim(s1)//"_"//trim(s2)//"_-1"//"_petsc.dat"
            if (ncell_c(3) .ne. 0) then
              call MatDuplicate(p_dmatxy(i1, i2, -1), MAT_SHARE_NONZERO_PATTERN, p_h10_cc(i1, i2), ierr)
              call petsc_mat_load(p_h10_cc(i1, i2), tmpfile, .false., .false., nrow_pproc, mattype_sparse)
    
              tmpfile = trim(ecc_dir)//"/S_"//trim(s1)//"_"//trim(s2)//"_-1"//"_petsc.dat"
              call MatDuplicate(p_dmatxy(i1, i2, -1), MAT_SHARE_NONZERO_PATTERN, p_s10_cc(i1, i2), ierr)
              call petsc_mat_load(p_s10_cc(i1, i2), tmpfile, .false., .false., nrow_pproc, mattype_sparse)
    
              tmpfile = trim(ecc_dir)//"/H_"//trim(spin_str)//"_"//trim(s1)//"_"//trim(s2)//"_1"//"_petsc.dat"
              call MatDuplicate(p_dmatxy(i1, i2, 1), MAT_SHARE_NONZERO_PATTERN, p_h01_cc(i1, i2), ierr)
              call petsc_mat_load(p_h01_cc(i1, i2), tmpfile, .false., .false., nrow_pproc, mattype_sparse)
    
              tmpfile = trim(ecc_dir)//"/S_"//trim(s1)//"_"//trim(s2)//"_1"//"_petsc.dat"
              call MatDuplicate(p_dmatxy(i1, i2, 1), MAT_SHARE_NONZERO_PATTERN, p_s01_cc(i1, i2), ierr)
              call petsc_mat_load(p_s01_cc(i1, i2), tmpfile, .false., .false., nrow_pproc, mattype_sparse)
            end if
          end if

          p_h00_diag(i1,i2,0)%p=>p_h00_cc(i1,i2)
          p_s00_diag(i1,i2,0)%p=>p_s00_cc(i1,i2)
          
          if (ncell_c(3) .ne. 0) then
            p_h00_diag(i1,i2,-1)%p=>p_h10_cc(i1,i2)
            p_h00_diag(i1,i2,1)%p=>p_h01_cc(i1,i2)          
            p_s00_diag(i1,i2,-1)%p=>p_s10_cc(i1,i2)
            p_s00_diag(i1,i2,1)%p=>p_s01_cc(i1,i2)
          end if
          if (l_dftu) then
            do i3 = -ncell_c(3), ncell_c(3)
              p_k00_diag(i1,i2,i3)%p=>p_dmatxy(i1, i2, i3)  
            end do
          end if
        end if


        ! S-matrix
        tmpfile = trim(ecc_dir)//"/S_"//trim(s1)//"_"//trim(s2)//"_0"//"_petsc.dat"
        call MatDuplicate(p_dmatxy(i1, i2, 0), MAT_SHARE_NONZERO_PATTERN, p_s00_cc(i1, i2), ierr)
        call petsc_mat_load(p_s00_cc(i1, i2), tmpfile, .false., .false., nrow_pproc, mattype_sparse)

        if (l_ep) then

          call init_ep_active_mat(p_dmatxy(i1, i2, 0), p_tmp_ep_mat,&
          &ep_active_atoms(1), ep_active_atoms(2), imu_to_at)

          call MatDuplicate(p_dmatxy(i1, i2, 0), MAT_SHARE_NONZERO_PATTERN, p_tmp_read, ierr)

          if (l_ep_reuse_lambda) then
            call check_lambdas_on_disk(l_read_dh)
          else
            l_read_dh = .true.
          end if

          write (pstr_out, fmt='(A,l)') "read dh", l_read_dh; call petsc_print_master()
          if (l_read_dh) then
            do iat = ep_active_atoms(1), ep_active_atoms(2)
              atstr = int2str(iat)
              do ixyz = 1, 3
                xyzstr = int2str(ixyz)
                tmpfile = trim(ecc_dir)//"/dH_"//trim(atstr)//"_"//trim(xyzstr)//"_"//trim(s1)//"_"//trim(s2)//"_0"//"_petsc.dat"
                call petsc_mat_load(p_tmp_read, tmpfile, .false., .false., nrow_pproc, mattype_sparse)
                call MatDuplicate(p_tmp_ep_mat, MAT_SHARE_NONZERO_PATTERN, p_dh_cc(i1, i2, iat, ixyz), ierr)
                call petsc_aXpY(p_dh_cc(i1, i2, iat, ixyz), p_tmp_read, p_one, PETSC_FALSE, .false.)
              end do
            end do
!~             call petsc_mat_info(p_tmp_ep_mat,"p_tmp_ep_mat ",ierr)
            call MatDestroy(p_tmp_read, ierr)
          end if

        end if
        
        if (calc_current_density) then
          ! Vnl-matrix
          tmpfile = trim(ecc_dir)//"/Vnl_"//trim(s1)//"_"//trim(s2)//"_0"//"_petsc.dat"        
          call MatDuplicate(p_dmatxy(i1, i2, 0), MAT_SHARE_NONZERO_PATTERN, p_vnl00_cc(i1, i2), ierr)
          call petsc_mat_load(p_vnl00_cc(i1, i2), tmpfile, .false., .false., nrow_pproc, mattype_sparse)          
        end if

      end do
    end do


    write (pstr_out, fmt='(A,i0)') "total NNZ ", nmat_c; call petsc_print_master()
!work around to extend system, remove this at some point
    nmu_l = sum(imu_ecc(1:nat_elec_l))
    nmu_r = sum(imu_ecc(nat_ecc - nat_elec_r + 1:nat_ecc))

! -----------------------------------
! the nonzero structure could be matched to H_ecc_LL, ...
!
! -----------------------------------

    if (l_use_sigma_l) then
! Left electrode
      pstr_out = "read data of left electrode..."; call petsc_print_master()
      allocate (nrow_pproc_elec_l(2, nprocs), stat=ierr)
      if (ierr .ne. 0) then
        write (errormsg, fmt='(A,i8)') "allocation error init_sys ", ierr
        call error()
      end if
!~       s1="0"
!~       s2="0"
!~       tmpfile=trim(elec_l_dir)//"/K_"//trim(spin_str)//"_"//trim(s1)//"_"//trim(s2)//"_0"//"_petsc.dat"
!~       call petsc_mat_load(p_k00_l(0,0),tmpfile,.false.,.true.,nrow_pproc_elec_l,mattype_electrode)
      nzrow_l = nmu_l
      nmat_l = 0
      call nzs_to_procs(nrow_pproc_elec_l, nzrow_l)
      do i1 = -ncell_l(1), ncell_l(1)
        do i2 = -ncell_l(2), ncell_l(2)
  
          s1 = int2str(i1)
          s2 = int2str(i2)
  
          ! K_L_00
          tmpfile = trim(elec_l_dir)//"/K_"//trim(spin_str)//"_"//trim(s1)//"_"//trim(s2)//"_0"//"_petsc.dat"
          call petsc_mat_load(p_k00_l(i1, i2), tmpfile, .true., .false., nrow_pproc_elec_l, &
          &  mattype_electrode, nnz_out = nnz)
          nmat_l = nmat_l + nnz
          
          ! S_L_00
          tmpfile = trim(elec_l_dir)//"/S_"//trim(s1)//"_"//trim(s2)//"_0"//"_petsc.dat"
          call MatDuplicate(p_k00_l(i1, i2), MAT_SHARE_NONZERO_PATTERN, p_s00_l(i1, i2), ierr)
          call petsc_mat_load(p_s00_l(i1, i2), tmpfile, .false., .false., nrow_pproc_elec_l, mattype_electrode)
  
          ! H_L_00
          tmpfile = trim(elec_l_dir)//"/H_"//trim(spin_str)//"_"//trim(s1)//"_"//trim(s2)//"_0"//"_petsc.dat"
          call MatDuplicate(p_k00_l(i1, i2), MAT_SHARE_NONZERO_PATTERN, p_h00_l(i1, i2), ierr)
          call petsc_mat_load(p_h00_l(i1, i2), tmpfile, .false., .false., nrow_pproc_elec_l, mattype_electrode)
  
          ! K_L_01
          tmpfile = trim(elec_l_dir)//"/K_"//trim(spin_str)//"_"//trim(s1)//"_"//trim(s2)//"_1"//"_petsc.dat"
          call petsc_mat_load(p_k01_l(i1, i2), tmpfile, .true., .false., nrow_pproc_elec_l, &
        &  mattype_electrode, nnz_out = nnz)
          nmat_l = nmat_l + nnz
  
          ! S_L_01
          tmpfile = trim(elec_l_dir)//"/S_"//trim(s1)//"_"//trim(s2)//"_1"//"_petsc.dat"
          call MatDuplicate(p_k01_l(i1, i2), MAT_SHARE_NONZERO_PATTERN, p_s01_l(i1, i2), ierr)
          call petsc_mat_load(p_s01_l(i1, i2), tmpfile, .false., .false., nrow_pproc_elec_l, mattype_electrode)
  
          ! H_L_01 ( actually H_L-10)
          tmpfile = trim(elec_l_dir)//"/H_"//trim(spin_str)//"_"//trim(s1)//"_"//trim(s2)//"_1"//"_petsc.dat"
          call MatDuplicate(p_k01_l(i1, i2), MAT_SHARE_NONZERO_PATTERN, p_h01_l(i1, i2), ierr)
          call petsc_mat_load(p_h01_l(i1, i2), tmpfile, .false., .false., nrow_pproc_elec_l, &
        &  mattype_electrode, nnz_out = nnz)
          nmat_l = nmat_l + nnz
  
          ! k_L_10
          tmpfile = trim(elec_l_dir)//"/K_"//trim(spin_str)//"_"//trim(s1)//"_"//trim(s2)//"_-1"//"_petsc.dat"
          call petsc_mat_load(p_k10_l(i1, i2), tmpfile, .true., .false., nrow_pproc_elec_l, mattype_electrode)
  
          ! S_L_10
          tmpfile = trim(elec_l_dir)//"/S_"//trim(s1)//"_"//trim(s2)//"_-1"//"_petsc.dat"
          call MatDuplicate(p_k10_l(i1, i2), MAT_SHARE_NONZERO_PATTERN, p_s10_l(i1, i2), ierr)
          call petsc_mat_load(p_s10_l(i1, i2), tmpfile, .false., .false., nrow_pproc_elec_l, mattype_electrode)
  
          ! H_L_10
          tmpfile = trim(elec_l_dir)//"/H_"//trim(spin_str)//"_"//trim(s1)//"_"//trim(s2)//"_-1"//"_petsc.dat"
          call MatDuplicate(p_k10_l(i1, i2), MAT_SHARE_NONZERO_PATTERN, p_h10_l(i1, i2), ierr)
          call petsc_mat_load(p_h10_l(i1, i2), tmpfile, .false., .false., nrow_pproc_elec_l, mattype_electrode)
  
        end do
      end do
    

!~ in case the neighbor for the electrodes are smaller than for CC add dummy elements and zero them
!~ the opposite case should not happen as the electrodes are contained in CC
      do i1 = -ncell_c(1), ncell_c(1)
        do i2 = -ncell_c(2), ncell_c(2)
          if ((abs(i1).gt.ncell_l(1)).or.(abs(i2).gt.ncell_l(2))) then
            call MatDuplicate(p_k00_l(0, 0), MAT_SHARE_NONZERO_PATTERN, p_s00_l(i1, i2), ierr)
            call MatZeroEntries(p_s00_l(i1, i2), ierr)
            call MatDuplicate(p_k00_l(0, 0), MAT_SHARE_NONZERO_PATTERN, p_h00_l(i1, i2), ierr)
            call MatZeroEntries(p_h00_l(i1, i2), ierr)
            call MatDuplicate(p_k10_l(0, 0), MAT_SHARE_NONZERO_PATTERN, p_s10_l(i1, i2), ierr)
            call MatZeroEntries(p_s10_l(i1, i2), ierr)
            call MatDuplicate(p_k10_l(0, 0), MAT_SHARE_NONZERO_PATTERN, p_h10_l(i1, i2), ierr)
            call MatZeroEntries(p_h10_l(i1, i2), ierr)
            call MatDuplicate(p_k01_l(0, 0), MAT_SHARE_NONZERO_PATTERN, p_s01_l(i1, i2), ierr)
            call MatZeroEntries(p_s01_l(i1, i2), ierr)
            call MatDuplicate(p_k01_l(0, 0), MAT_SHARE_NONZERO_PATTERN, p_h01_l(i1, i2), ierr)
            call MatZeroEntries(p_h01_l(i1, i2), ierr)
          end if
        end do
      end do      
    end if

! Right electrode
    if (l_use_sigma_r) then
      pstr_out = "read data of right electrode "; call petsc_print_master()
      allocate (nrow_pproc_elec_r(2, nprocs), stat=ierr)
      if (ierr .ne. 0) then
        write (errormsg, *) "allocation error init_sys ", ierr
        call error()
      end if
!~       s1="0"
!~       s2="0"
!~       tmpfile=trim(elec_r_dir)//"/K_"//trim(spin_str)//"_"//trim(s1)//"_"//trim(s2)//"_0"//"_petsc.dat"
!~       call petsc_mat_load(p_k00_l(0,0),tmpfile,.false.,.true.,nrow_pproc_elec_r,mattype_electrode)
      nzrow_r = nmu_r
      nmat_r = 0
      call nzs_to_procs(nrow_pproc_elec_r, nzrow_r)
      do i1 = -ncell_r(1), ncell_r(1)
        do i2 = -ncell_r(2), ncell_r(2)
  
          s1 = int2str(i1)
          s2 = int2str(i2)
  
          ! K_r_00
          tmpfile = trim(elec_r_dir)//"/K_"//trim(spin_str)//"_"//trim(s1)//"_"//trim(s2)//"_0"//"_petsc.dat"
          call petsc_mat_load(p_k00_r(i1, i2), tmpfile, .true., .false., nrow_pproc_elec_r, &
        &  mattype_electrode, nnz_out = nnz)
          nmat_r = nmat_r + nnz
  
          ! S_r_00
          tmpfile = trim(elec_r_dir)//"/S_"//trim(s1)//"_"//trim(s2)//"_0"//"_petsc.dat"
          call MatDuplicate(p_k00_r(i1, i2), MAT_SHARE_NONZERO_PATTERN, p_s00_r(i1, i2), ierr)
          call petsc_mat_load(p_s00_r(i1, i2), tmpfile, .false., .false., nrow_pproc_elec_r, mattype_electrode)
  
          ! H_r_00
          tmpfile = trim(elec_r_dir)//"/H_"//trim(spin_str)//"_"//trim(s1)//"_"//trim(s2)//"_0"//"_petsc.dat"
          call MatDuplicate(p_k00_r(i1, i2), MAT_SHARE_NONZERO_PATTERN, p_h00_r(i1, i2), ierr)
          call petsc_mat_load(p_h00_r(i1, i2), tmpfile, .false., .false., nrow_pproc_elec_r, mattype_electrode)
  
          ! K_r_01
          tmpfile = trim(elec_r_dir)//"/K_"//trim(spin_str)//"_"//trim(s1)//"_"//trim(s2)//"_1"//"_petsc.dat"
          call petsc_mat_load(p_k01_r(i1, i2), tmpfile, .true., .false., nrow_pproc_elec_r, &
        &  mattype_electrode, nnz_out = nnz)
          nmat_r = nmat_r + nnz
  
          ! S_r_01
          tmpfile = trim(elec_r_dir)//"/S_"//trim(s1)//"_"//trim(s2)//"_1"//"_petsc.dat"
          call MatDuplicate(p_k01_r(i1, i2), MAT_SHARE_NONZERO_PATTERN, p_s01_r(i1, i2), ierr)
          call petsc_mat_load(p_s01_r(i1, i2), tmpfile, .false., .false., nrow_pproc_elec_r, mattype_electrode)
  
          ! H_r_01
          tmpfile = trim(elec_r_dir)//"/H_"//trim(spin_str)//"_"//trim(s1)//"_"//trim(s2)//"_1"//"_petsc.dat"
          call MatDuplicate(p_k01_r(i1, i2), MAT_SHARE_NONZERO_PATTERN, p_h01_r(i1, i2), ierr)
          call petsc_mat_load(p_h01_r(i1, i2), tmpfile, .false., .false., nrow_pproc_elec_r, mattype_electrode)
  
          ! K_r_10 ( actually K_r-10)
          tmpfile = trim(elec_r_dir)//"/K_"//trim(spin_str)//"_"//trim(s1)//"_"//trim(s2)//"_-1"//"_petsc.dat"
          call petsc_mat_load(p_k10_r(i1, i2), tmpfile, .true., .false., nrow_pproc_elec_r, &
        &  mattype_electrode, nnz_out = nnz)
          nmat_r = nmat_r + nnz        
  
          ! S_r_10
          tmpfile = trim(elec_r_dir)//"/S_"//trim(s1)//"_"//trim(s2)//"_-1"//"_petsc.dat"
          call MatDuplicate(p_k10_r(i1, i2), MAT_SHARE_NONZERO_PATTERN, p_s10_r(i1, i2), ierr)
          call petsc_mat_load(p_s10_r(i1, i2), tmpfile, .false., .false., nrow_pproc_elec_r, mattype_electrode)
  
          ! H_r_10
          tmpfile = trim(elec_r_dir)//"/H_"//trim(spin_str)//"_"//trim(s1)//"_"//trim(s2)//"_-1"//"_petsc.dat"
          call MatDuplicate(p_k10_r(i1, i2), MAT_SHARE_NONZERO_PATTERN, p_h10_r(i1, i2), ierr)
          call petsc_mat_load(p_h10_r(i1, i2), tmpfile, .false., .false., nrow_pproc_elec_r, mattype_electrode)
  
        end do
      end do

!~ in case the neighbor for the electrodes are smaller than for CC add dummy elements and zero them
!~ the opposite case should not happen as the electrodes are contained in CC
      do i1 = -ncell_c(1), ncell_c(1)
        do i2 = -ncell_c(2), ncell_c(2)
          if ((abs(i1).gt.ncell_r(1)).or.(abs(i2).gt.ncell_r(2))) then
            call MatDuplicate(p_k00_r(0, 0), MAT_SHARE_NONZERO_PATTERN, p_s00_r(i1, i2), ierr)
            call MatZeroEntries(p_s00_r(i1, i2), ierr)
            call MatDuplicate(p_k00_r(0, 0), MAT_SHARE_NONZERO_PATTERN, p_h00_r(i1, i2), ierr)
            call MatZeroEntries(p_h00_r(i1, i2), ierr)
            call MatDuplicate(p_k10_r(0, 0), MAT_SHARE_NONZERO_PATTERN, p_s10_r(i1, i2), ierr)
            call MatZeroEntries(p_s10_r(i1, i2), ierr)
            call MatDuplicate(p_k10_r(0, 0), MAT_SHARE_NONZERO_PATTERN, p_h10_r(i1, i2), ierr)
            call MatZeroEntries(p_h10_r(i1, i2), ierr)
            call MatDuplicate(p_k01_r(0, 0), MAT_SHARE_NONZERO_PATTERN, p_s01_r(i1, i2), ierr)
            call MatZeroEntries(p_s01_r(i1, i2), ierr)
            call MatDuplicate(p_k01_r(0, 0), MAT_SHARE_NONZERO_PATTERN, p_h01_r(i1, i2), ierr)
            call MatZeroEntries(p_h01_r(i1, i2), ierr)
          end if
        end do
      end do      
    end if



    call system_clock(countf)
    write (pstr_out, fmt='(A,e24.12)') "Matrices read. ", real(countf - counti, 8)/real(count_rate, 8); call petsc_print_master()

    
    
!~     write(0,*) "nmat_l ",nmat_l
!~     write(0,*) "nmat_r ",nmat_r
!~     write(0,*) "nmat_c ",nmat_c
    nmat_l3 = nmat_l * 3
    nmat_r3 = nmat_r * 3
    nat_elec_l3 = nat_elec_l * 3
    nat_elec_r3 = nat_elec_r * 3
    allocate (tomat(8, nmat_c), inao_c(nat_ecc), neigh_c(nat_ecc), ndim_c(nat_ecc), &
   &  tomatr(3, nmat_c), atoms_c(nat_ecc), stat=ierr)
    if (ierr .ne. 0) then
      write (errormsg, *) "allocation error init_sys ", ierr
      call error()
    end if
    allocate (tomat_l(8, nmat_l3), inao_l(nat_elec_l3), neigh_l(nat_elec_l3), &
   &  tomatr_l(3, nmat_l3), ndim_l(nat_elec_l3), atoms_l(nat_elec_l3), &
   &  xyz_elec_l3(3, nat_elec_l3),stat=ierr)
    if (ierr .ne. 0) then
      write (errormsg, *) "allocation error init_sys ", ierr
      call error()
    end if
    allocate (tomat_r(8, nmat_r3), inao_r(nat_elec_r3), neigh_r(nat_elec_r3), &
   &  tomatr_r(3, nmat_r3), ndim_r(nat_elec_r3), atoms_r(nat_elec_r3), &
   &  xyz_elec_r3(3, nat_elec_r3), stat=ierr)
    if (ierr .ne. 0) then
      write (errormsg, *) "allocation error init_sys ", ierr
      call error()
    end if

    ! work around for reduced ecc in DFT calculation
    nmu_l = 0
    nmu_r = 0
    pstr_out = "read_conquest_info "; call petsc_print_master()
    call read_conquest_info(ecc_dir, xyz_ecc, dlat_c, nat_ecc, tomat, tomatr, inao_c, neigh_c, nlist_c, ndim_c, atoms_c)
    nmu_l = sum(imu_ecc(1:nat_elec_l))
    nmu_r = sum(imu_ecc(nat_ecc - nat_elec_r + 1:nat_ecc))
    
    dlat_l3 = dlat_l
    dlat_l3(3,3) = dlat_l3(3,3) * 3d0
    dlat_r3 = dlat_r
    dlat_r3(3,3) = dlat_r3(3,3) * 3d0
    
    do i = 1, nat_elec_l
      xyz_elec_l3(1:3, i) = xyz_elec_l(1:3, i)
      xyz_elec_l3(1:3, nat_elec_l + i) = xyz_elec_l(1:3, i) + [ 0.0d0, 0.0d0, dlat_l(3,3) ]
      xyz_elec_l3(1:3, nat_elec_l * 2 + i) = xyz_elec_l(1:3, i) + [ 0.0d0, 0.0d0, dlat_l(3,3)*2d0 ]
    end do
    do i = 1, nat_elec_r
      xyz_elec_r3(1:3, i) = xyz_elec_r(1:3, i)
      xyz_elec_r3(1:3, nat_elec_r + i) = xyz_elec_r(1:3, i) + [ 0.0d0, 0.0d0, dlat_r(3,3) ]
      xyz_elec_r3(1:3, nat_elec_r * 2 + i) = xyz_elec_r(1:3, i) + [ 0.0d0, 0.0d0, dlat_r(3,3)*2d0 ]
    end do
    
    if (l_reaktor) then
      if (l_use_sigma_l) call read_conquest_info(trim(elec_l_dir)//"_reaktor", &
     & xyz_elec_l3, dlat_l3, nat_elec_l3, tomat_l, tomatr_l, inao_l, neigh_l, &
     & nlist_l, ndim_l, atoms_l)
      if (l_use_sigma_r )call read_conquest_info(trim(elec_r_dir)//"_reaktor", &
     & xyz_elec_r3, dlat_r3, nat_elec_r3, tomat_r, tomatr_r, inao_r, neigh_r, &
     & nlist_r, ndim_r, atoms_r)
    end if

    diag_dim = diag_dim_save 
    
    if ((k_mat_mode .eq. 2)) then
      if (diag_dim .eq. 1) then
        ncell_c(2:3) = 0
        ncell_l(2:3) = 0
        ncell_r(2:3) = 0
      else if (diag_dim .eq. 2) then
        ncell_c(3) = 0
        ncell_l(3) = 0
        ncell_r(3) = 0
      end if
    end if


    if ((l_dftu) .and. (k_mat_mode .eq. 2)) then
      p_scal = 0d0
      do i1 = -ncell_c(1), ncell_c(1)
        do i2 = -ncell_c(2), ncell_c(2)
          do i3 = -ncell_c(3), ncell_c(3)
            call MatMatMult(p_k00_diag(i1,i2,i3)%p, p_s00_diag(-i1,-i2,-i3)%p, &
           & MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, p_tmp, ierr)
            call MatGetTrace(p_tmp, trace_c, ierr)
            p_scal = p_scal + trace_c
            call MatDestroy(p_tmp, ierr)
          end do
        end do
      end do
      write (pstr_out, fmt='(A,2e24.12)') "N ", p_scal
      call petsc_print_master()
    end if
    


    if (l_ep) then
      write (errormsg, fmt='(A)') "EP does not work with k-on-demand at present"
      call error()
    end if

    if (k_mat_mode .eq. 1) then
      call MatDuplicate(p_h00_cc(0, 0), MAT_DO_NOT_COPY_VALUES, p_tmp, ierr)
      call fourier_trans_ik(p_h00_cc, p_tmp, kp_c, 1, dlat_c, ncell_c(1), ncell_c(2), 0)
      call MatDuplicate(p_tmp, MAT_COPY_VALUES, p_h00_ik_cc(1), ierr)
      call petsc_get_alloc_preG(p_tmp, p_preG, nmu_l, nmu_r)
      call petsc_get_a_with_b_c(p_invGr, p_tmp, p_tmp, mattype_sparse)
      call MatPreallocatorPreallocate(p_preG, PETSC_TRUE, p_invGr, ierr)
      call MatSetOption (p_invGr, MAT_NO_OFF_PROC_ENTRIES, PETSC_FALSE, ierr)        
      call petsc_mat_info(p_tmp, "Hk", ierr)
      call petsc_mat_info(p_invGr, "Gr", ierr)
      call MatZeroEntries(p_invGr, ierr)
      call MatDestroy(p_tmp, ierr)
      call MatDestroy(p_preG, ierr)

      call MatDuplicate(p_s00_cc(0, 0), MAT_DO_NOT_COPY_VALUES, p_tmp, ierr)
      call fourier_trans_ik(p_s00_cc, p_tmp, kp_c, 1, dlat_c, ncell_c(1), ncell_c(2), 0)
      call MatDuplicate(p_tmp, MAT_DO_NOT_COPY_VALUES, p_s00_ik_cc(1), ierr)
      call MatDestroy(p_tmp, ierr)
      
      if (calc_current_density) then
        call MatDuplicate(p_vnl00_cc(0, 0), MAT_DO_NOT_COPY_VALUES, p_tmp, ierr)
        call fourier_trans_ik(p_vnl00_cc, p_tmp, kp_c, 1, dlat_c, ncell_c(1), ncell_c(2), 0)
        call MatDuplicate(p_tmp, MAT_DO_NOT_COPY_VALUES, p_vnl00_ik_cc(1), ierr)
        call MatDestroy(p_tmp, ierr)
      end if

      if (l_use_sigma_l) then
        call MatDuplicate(p_h00_l(0, 0), MAT_DO_NOT_COPY_VALUES, p_tmp, ierr)
        call fourier_trans_ik(p_h00_l, p_tmp, kp_c, 1, dlat_c, ncell_c(1), ncell_c(2), 0)
        call MatDuplicate(p_tmp, MAT_DO_NOT_COPY_VALUES, p_h00_ik_l(1), ierr)
        call MatDestroy(p_tmp, ierr)

        call MatDuplicate(p_s00_l(0, 0), MAT_DO_NOT_COPY_VALUES, p_tmp, ierr)
        call fourier_trans_ik(p_s00_l, p_tmp, kp_c, 1, dlat_c, ncell_c(1), ncell_c(2), 0)
        call MatDuplicate(p_tmp, MAT_DO_NOT_COPY_VALUES, p_s00_ik_l(1), ierr)
        call MatDestroy(p_tmp, ierr)

        call MatDuplicate(p_h10_l(0, 0), MAT_DO_NOT_COPY_VALUES, p_tmp, ierr)
        call fourier_trans_ik(p_h10_l, p_tmp, kp_c, 1, dlat_c, ncell_c(1), ncell_c(2), 0)
        call MatDuplicate(p_tmp, MAT_DO_NOT_COPY_VALUES, p_h10_ik_l(1), ierr)
        call MatDestroy(p_tmp, ierr)

        call MatDuplicate(p_s10_l(0, 0), MAT_DO_NOT_COPY_VALUES, p_tmp, ierr)
        call fourier_trans_ik(p_s10_l, p_tmp, kp_c, 1, dlat_c, ncell_c(1), ncell_c(2), 0)
        call MatDuplicate(p_tmp, MAT_DO_NOT_COPY_VALUES, p_s10_ik_l(1), ierr)
        call MatDestroy(p_tmp, ierr)

        call MatDuplicate(p_h01_l(0, 0), MAT_DO_NOT_COPY_VALUES, p_tmp, ierr)
        call fourier_trans_ik(p_h01_l, p_tmp, kp_c, 1, dlat_c, ncell_c(1), ncell_c(2), 0)
        call MatDuplicate(p_tmp, MAT_DO_NOT_COPY_VALUES, p_h01_ik_l(1), ierr)
        call MatDestroy(p_tmp, ierr)

        call MatDuplicate(p_s01_l(0, 0), MAT_DO_NOT_COPY_VALUES, p_tmp, ierr)
        call fourier_trans_ik(p_s01_l, p_tmp, kp_c, 1, dlat_c, ncell_c(1), ncell_c(2), 0)
        call MatDuplicate(p_tmp, MAT_DO_NOT_COPY_VALUES, p_s01_ik_l(1), ierr)
        call MatDestroy(p_tmp, ierr)
      end if

      if (l_use_sigma_r) then
        call MatDuplicate(p_h00_r(0, 0), MAT_DO_NOT_COPY_VALUES, p_tmp, ierr)
        call fourier_trans_ik(p_h00_r, p_tmp, kp_c, 1, dlat_c, ncell_c(1), ncell_c(2), 0)
        call MatDuplicate(p_tmp, MAT_DO_NOT_COPY_VALUES, p_h00_ik_r(1), ierr)
        call MatDestroy(p_tmp, ierr)

        call MatDuplicate(p_s00_r(0, 0), MAT_DO_NOT_COPY_VALUES, p_tmp, ierr)
        call fourier_trans_ik(p_s00_r, p_tmp, kp_c, 1, dlat_c, ncell_c(1), ncell_c(2), 0)
        call MatDuplicate(p_tmp, MAT_DO_NOT_COPY_VALUES, p_s00_ik_r(1), ierr)
        call MatDestroy(p_tmp, ierr)

        call MatDuplicate(p_h10_r(0, 0), MAT_DO_NOT_COPY_VALUES, p_tmp, ierr)
        call fourier_trans_ik(p_h10_r, p_tmp, kp_c, 1, dlat_c, ncell_c(1), ncell_c(2), 0)
        call MatDuplicate(p_tmp, MAT_DO_NOT_COPY_VALUES, p_h10_ik_r(1), ierr)
        call MatDestroy(p_tmp, ierr)

        call MatDuplicate(p_s10_r(0, 0), MAT_DO_NOT_COPY_VALUES, p_tmp, ierr)
        call fourier_trans_ik(p_s10_r, p_tmp, kp_c, 1, dlat_c, ncell_c(1), ncell_c(2), 0)
        call MatDuplicate(p_tmp, MAT_DO_NOT_COPY_VALUES, p_s10_ik_r(1), ierr)
        call MatDestroy(p_tmp, ierr)

        call MatDuplicate(p_h01_r(0, 0), MAT_DO_NOT_COPY_VALUES, p_tmp, ierr)
        call fourier_trans_ik(p_h01_r, p_tmp, kp_c, 1, dlat_c, ncell_c(1), ncell_c(2), 0)
        call MatDuplicate(p_tmp, MAT_DO_NOT_COPY_VALUES, p_h01_ik_r(1), ierr)
        call MatDestroy(p_tmp, ierr)

        call MatDuplicate(p_s01_r(0, 0), MAT_DO_NOT_COPY_VALUES, p_tmp, ierr)
        call fourier_trans_ik(p_s01_r, p_tmp, kp_c, 1, dlat_c, ncell_c(1), ncell_c(2), 0)
        call MatDuplicate(p_tmp, MAT_DO_NOT_COPY_VALUES, p_s01_ik_r(1), ierr)
        call MatDestroy(p_tmp, ierr)
      end if
                
    else if (k_mat_mode .eq. 2) then
      call MatDuplicate(p_h00_diag(0, 0, 0)%p, MAT_DO_NOT_COPY_VALUES, p_tmp, ierr)
      call fourier3d_trans_ik(p_h00_diag, p_tmp, kp_c, 1, dlat_c, ncell_c(1), ncell_c(2), ncell_c(3), 0)
      call MatDuplicate(p_tmp, MAT_DO_NOT_COPY_VALUES, p_h00_ik_cc(1), ierr)
!~         call petsc_get_alloc_preG(p_tmp, p_preG, nmu_l, nmu_r)
      call MatDestroy(p_tmp, ierr)

      call MatDuplicate(p_s00_diag(0, 0, 0)%p, MAT_DO_NOT_COPY_VALUES, p_tmp, ierr)
      call fourier3d_trans_ik(p_s00_diag, p_tmp, kp_c, 1, dlat_c, ncell_c(1), ncell_c(2), ncell_c(3), 0)
      call MatDuplicate(p_tmp, MAT_DO_NOT_COPY_VALUES, p_s00_ik_cc(1), ierr)
      call MatDestroy(p_tmp, ierr)        
      
    else if (k_mat_mode .eq. 3) then
      call MatDuplicate(p_dmatxy(0, 0, 0), MAT_DO_NOT_COPY_VALUES, p_tmp, ierr)
      call fourier_trans_ik(p_dmatxy(:, :, 0), p_tmp, kp_c, 1, dlat_c, ncell_c(1), ncell_c(2), 0)
      call MatDuplicate(p_tmp, MAT_DO_NOT_COPY_VALUES, p_k00_ik_cc(1), ierr)
!~         call petsc_get_alloc_preG(p_tmp, p_preG, nmu_l, nmu_r)
      call MatDestroy(p_tmp, ierr)

      call MatDuplicate(p_s00_cc(0, 0), MAT_DO_NOT_COPY_VALUES, p_tmp, ierr)
      call fourier_trans_ik(p_s00_cc, p_tmp, kp_c, 1, dlat_c, ncell_c(1), ncell_c(2), 0)
      call MatDuplicate(p_tmp, MAT_DO_NOT_COPY_VALUES, p_s00_ik_cc(1), ierr)
      call MatDestroy(p_tmp, ierr)

    end if

    if (l_dftu) then
      if (.not.(k_mat_mode .eq. 3)) then
        call MatDuplicate(p_dmatxy(0, 0, 0), MAT_DO_NOT_COPY_VALUES, p_tmp, ierr)
        call fourier_trans_ik(p_dmatxy(:, :, 0), p_tmp, kp_c, 1, dlat_c, &
       &  ncell_c(1), ncell_c(2), 0)
        call MatDuplicate(p_tmp, MAT_DO_NOT_COPY_VALUES, p_k00_ik_cc(1), ierr)
!~         call petsc_get_alloc_preG(p_tmp, p_preG, nmu_l, nmu_r)
        call MatDestroy(p_tmp, ierr)
      end if
      call read_dftu(dftu_file_c, species_c)
!~       call read_dftu(dftu_file_l, species_l)
!~       call read_dftu(dftu_file_r, species_r)
    end if

!~       if ((k_mat_mode .eq. 1) .or. (l_diag_fixH)) then
    if (l_diag_fixH) then
      call fix_mat(p_h00_cc, p_h00_l, p_h00_r, p_h10_cc, p_h10_l, p_h01_cc,  p_h01_r)
      call fix_mat(p_s00_cc, p_s00_l, p_s00_r, p_s10_cc, p_s10_l, p_s01_cc,  p_s01_r)
    end if

    if ((k_mat_mode .eq. 2).and.(diag_dim .eq. 3).and.(1.eq.2)) then
      p_scal=cmplx(0.5d0,0d0)
      do i1 = -ncell_c(1), ncell_c(1)
        do i2 = -ncell_c(2), ncell_c(2)
          do i3 = -ncell_c(3), ncell_c(3)

            call MatTranspose(p_h00_diag(-i1,-i2,-i3)%p, MAT_INITIAL_MATRIX, &
           &  p_mat_tmp, ierr)

            call MatAXPY(p_h00_diag(i1, i2, i3)%p, p_one, p_mat_tmp , &
           &  DIFFERENT_NONZERO_PATTERN, ierr)
            call MatScale(p_h00_diag(i1,i2,i3)%p,p_scal,ierr)
            call MatDestroy(p_mat_tmp,ierr)

            if (.not.((i1.eq.0).and.(i2.eq.0).and.(i3.eq.0))) then
              call MatTranspose(p_h00_diag(i1,i2,i3)%p, MAT_INITIAL_MATRIX, &
             &  p_mat_tmp, ierr)
              call MatAYPX(p_h00_diag(-i1, -i2, -i3)%p, p_zero, p_mat_tmp, &
             &  DIFFERENT_NONZERO_PATTERN, ierr)
              call MatDestroy(p_mat_tmp,ierr)
            end if

            call MatTranspose(p_s00_diag(-i1,-i2,-i3)%p, MAT_INITIAL_MATRIX, &
           &  p_mat_tmp, ierr)

            call MatAXPY(p_s00_diag(i1, i2, i3)%p, p_one, p_mat_tmp, &
           &  DIFFERENT_NONZERO_PATTERN, ierr)
            call MatScale(p_s00_diag(i1,i2,i3)%p,p_scal,ierr)
            call MatDestroy(p_mat_tmp,ierr)

            if (.not.((i1.eq.0).and.(i2.eq.0).and.(i3.eq.0))) then
              call MatTranspose(p_s00_diag(i1,i2,i3)%p, MAT_INITIAL_MATRIX, &
             &  p_mat_tmp, ierr)
              call MatAYPX(p_s00_diag(-i1, -i2, -i3)%p, p_zero, p_mat_tmp, &
             &  DIFFERENT_NONZERO_PATTERN, ierr)
              call MatDestroy(p_mat_tmp,ierr)
            end if

          end do
        end do
      end do
    end if

    if (l_ep) then
      allocate (p_ep_lambda_k(n_ep_modes_max, nk_c), &
                stat=ierr)
      if (ierr .ne. 0) then
        write (errormsg, fmt='(A,i8)') "allocation failure init_sys ", ierr
        call error()
      end if

      do i1 = 1, nk_c
        do i2 = 1, n_ep_modes_max
          call MatDuplicate(p_tmp_ep_mat, MAT_SHARE_NONZERO_PATTERN, p_ep_lambda_k(i2, i1), ierr)
          call MatZeroEntries(p_ep_lambda_k(i2, i1), ierr)
        end do
        call get_lambda_ep(i1)
!~           call petsc_mat_info(p_dhk_cc(i1,ep_active_atoms(1),1),"p_dhk_cc ",ierr)
        if (l_read_dh) then
          do i2 = ep_active_atoms(1), ep_active_atoms(2)
            do i3 = 1, 3
              call MatDestroy(p_dhk_cc(i1, i2, i3), ierr)
            end do
          end do
        end if
        call MatDestroy(p_phonon_EV(i1), ierr)
      end do

      call MatDestroy(p_invsqrt_mass, ierr)
      call MatDestroy(p_mass_matrix, ierr)

    end if
    
    if (l_dump_nzs) then
      call dump_nonzero_structure(p_h00k_cc(1), "h_cc_nzs.dat", PETSC_COMM_WORLD, 0)
      if (l_use_sigma_l) call dump_nonzero_structure(p_h00k_l(1), "h_ll_nzs.dat", PETSC_COMM_WORLD, 0)
      if (l_use_sigma_r) call dump_nonzero_structure(p_h00k_r(1), "h_rr_nzs.dat", PETSC_COMM_WORLD, 0)
    end if

    call MPI_BARRIER(PETSC_COMM_WORLD, ierr)
    
    if (l_dftu) then
      call init_Uoperator()  
    end if

!~     if (l_dftu .and. (k_mat_mode .eq. 1)) then
    if (l_dftu) then
      call add_dftu_to_H(.false.)
      if (l_diag_fixH) then
        call fix_mat(p_h00_cc, p_h00_l, p_h00_r, p_h10_cc, p_h10_l, p_h01_cc,  p_h01_r)
      end if      
    end if

    
!~     call test_ne()
    
!~     call dump_kmat_lr("K_negf_matrix2.i00.p000000", trim(elec_l_dir)//"_reaktor", "l")
!~     call dump_kmat_lr("K_negf_matrix2.i00.p000000", trim(elec_r_dir)//"_reaktor", "r")  
!~     call dump_kmat(.false., "K_negf_matrix2.i00.p000000", ecc_dir)
!~     do i1 = -ncell_c(1), ncell_c(1)
!~       write(s1,*) i1
!~       s1=adjustl(s1)
!~       do i2 = -ncell_c(2), ncell_c(2)
!~         write(s2,*) i2
!~         s2=adjustl(s2)
!~         do i3 = -ncell_c(3), ncell_c(3)
!~           write(s3,*) i3
!~           s3=adjustl(s3)
!~           tmpfile="Ktest_"//trim(s1)//"_"//trim(s2)//"_"//trim(s3)//".dat"
!~           call dump_nonzero_structure(p_dmatxy(i1,i2,i3), tmpfile , PETSC_COMM_WORLD, 0, 1)
!~         end do
!~       end do
!~     end do


!~     call petsc_cleanup()
!~     call slepcfinalize(ierr)
!~     call MPI_FINALIZE(ierr)    
!~     stop

  end subroutine init_sys
end module init_sys_mod
