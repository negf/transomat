module phonon_mod
#include <petsc/finclude/petsc.h>
  use petscmat
  implicit none

  integer :: ii_mode
  PetscScalar :: mode_ii_energy

contains

! only to check if we can use the adaptive integrator which takes PETSc mats to use as integrator for scalar values
! by simply creating a matrix of size 1. maybe a bit stupid but we can reuse the code.
  subroutine integrator_test()
#include <petsc/finclude/petsc.h>
    use petscmat
    use petsc_mod
    use petsc_wrapper
    use kinds
    use globals
    use integrator_mod
    implicit none

    real(dp) :: d_low, d_high, eps, dd_eps, err_eq_int
    PetscScalar :: p_zint1, p_zint2, pp(1), dx, x1, x2, x
    Mat :: p_d_tmp1(1), p_d_tmp2(1)
    integer :: ierr, ii(1), jj(1), ik(2), i, n

    call MatCreate(PETSC_COMM_WORLD, p_tmpcc1, ierr)
    call MatSetType(p_tmpcc1, mattype_sparse, ierr)
    call MatSetSizes(p_tmpcc1, PETSC_DECIDE, PETSC_DECIDE, 1, 1, ierr)
    if (inode .eq. 0) then
      ik(1) = 0
      ik(2) = 1
      pp = 0d0
      call MatMPIAIJSetPreallocationCSR(p_tmpcc1, ik, ik, pp, ierr)
    else
      call MatMPIAIJSetPreallocationCSR(p_tmpcc1, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, pp, ierr)
    end if

    call MatAssemblyBegin(p_tmpcc1, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(p_tmpcc1, MAT_FINAL_ASSEMBLY, ierr)
    call MatDuplicate(p_tmpcc1, MAT_DO_NOT_COPY_VALUES, p_d_tmp1(1), ierr)
    call MatDuplicate(p_tmpcc1, MAT_DO_NOT_COPY_VALUES, p_d_tmp2(1), ierr)

    call MatAssemblyBegin(p_tmpcc1, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(p_tmpcc1, MAT_FINAL_ASSEMBLY, ierr)
    call petsc_mat_info(p_tmpcc1, "p_tmpcc1 ", ierr)

!~       call petsc_vec_getvalue(0,mode_ii_energy,p_phonon_EW(1))
    eps = 1d-12
    dd_eps = 1d-8
    d_high = zsqrt(2d0*mode_ii_energy/dd_eps + mode_ii_energy*mode_ii_energy + eta_ph_cc*eta_ph_cc*0.25d0)
    d_low = -d_high
    nint_order = 64
    l_output_progress = .false.
    call adaptive_int3(dr_alpha_petsc, d_low, d_high, p_d_tmp1(1), p_d_tmp2(1), maxsub*1000, eps, err_eq_int)
    call petsc_mat_getvalue(0, 0, p_zint1, p_d_tmp1(1), 1, PETSC_COMM_WORLD)

    write (6, *) "p_zint1 ", p_zint1
    call MatDestroy(p_tmpcc1, ierr)

  end subroutine integrator_test

  function rho_eq_alpha(x)
    use globals
    use kinds
    implicit none

    real(dp) :: x
    complex(dp) :: rho_eq_alpha

    rho_eq_alpha = 1d0/pi*(eta_ph_cc*0.5d0/((x - mode_ii_energy)**2 + eta_ph_cc**2/4d0) - eta_ph_cc*0.5d0/((x + mode_ii_energy)**2 + eta_ph_cc**2/4d0))

  end function rho_eq_alpha

  function dr_alpha(x)
    use globals
    use kinds
    implicit none

    real(dp) :: x
    complex(dp) :: dr_alpha

    dr_alpha = 2d0*mode_ii_energy/(x*x - mode_ii_energy*mode_ii_energy + zione*eta_ph_cc*x - eta_ph_cc*eta_ph_cc*0.25d0)

  end function dr_alpha

  function dr_alpha_petsc(x)
    use globals
    use kinds
    use petsc_wrapper, only: petsc_matassemble
    use petsc_mod
    implicit none

    real(dp) :: x

    real(dp) :: z
    integer :: ii(1), ierr, dr_alpha_petsc
    PetscScalar :: pp(1)

    dr_alpha_petsc = -1

    if (inode .eq. 0) then
      z = x
      pp(1) = dr_alpha(z)
      ii = 0
!~         write(6,fmt='(A,6e24.12)') "pp ",x,pp,mode_ii_energy
      call MatSetValues(p_tmpcc1, 1, ii, 1, ii, pp, INSERT_VALUES, ierr)
    end if

    call MatAssemblyBegin(p_tmpcc1, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(p_tmpcc1, MAT_FINAL_ASSEMBLY, ierr)

    call petsc_matassemble(p_tmpcc1)

    dr_alpha_petsc = 0

  end function dr_alpha_petsc

  subroutine init_phonons()
#include <petsc/finclude/petsc.h>
    use petscmat
    use petsc_mod
    use petsc_wrapper
    use kinds
    use globals
    use misc

    implicit none

    integer :: ierr, i1, i2, iunit, ii(1)
    PetscViewer :: v_infile
    character(strln) :: infile, str_i1, str_i2
    PetscScalar :: p_in(1)

    allocate (p_Kphonon00_cc(-ncell_c(1):ncell_c(1), -ncell_c(2):ncell_c(2)), &
              stat=ierr)

    call PetscViewerCreate(PETSC_COMM_WORLD, v_infile, ierr)
    call PetscViewerSetType(v_infile, PETSCVIEWERBINARY, ierr)
    call PetscViewerFileSetMode(v_infile, FILE_MODE_READ, ierr)

    write (pstr_out, fmt='(A)') "read dynamical matrix "; call petsc_print_master()
    do i1 = -ncell_c(1), ncell_c(1)
      do i2 = -ncell_c(2), ncell_c(2)
        str_i1 = int2str(i1)
        str_i2 = int2str(i2)
        infile = "FCmat_"//trim(str_i1)//"_"//trim(str_i2)//"_0_petsc.dat"
        call PetscViewerFileSetName(v_infile, trim(infile), ierr)
!~           write(pstr_out,fmt='(A)') "read "//trim(infile) ; call petsc_print_master()
        call MatCreate(PETSC_COMM_WORLD, p_Kphonon00_cc(i1, i2), ierr)
        call MatSetType(p_Kphonon00_cc(i1, i2), mattype_sparse, ierr)
        call MatLoad(p_Kphonon00_cc(i1, i2), v_infile, ierr)
      end do
    end do

    write (pstr_out, fmt='(A)') "read MassMatrix^(-1/2) "; call petsc_print_master()
    infile = "sqrtMinv_petsc.dat"
    call PetscViewerFileSetName(v_infile, trim(infile), ierr)
    call MatCreate(PETSC_COMM_WORLD, p_invsqrt_mass, ierr)
    call MatSetType(p_invsqrt_mass, mattype_sparse, ierr)
    call MatLoad(p_invsqrt_mass, v_infile, ierr)

    write (pstr_out, fmt='(A)') "read MassMatrix "; call petsc_print_master()
    infile = "Massmat_petsc.dat"
    call PetscViewerFileSetName(v_infile, trim(infile), ierr)
    call MatCreate(PETSC_COMM_WORLD, p_mass_matrix, ierr)
    call MatSetType(p_mass_matrix, mattype_sparse, ierr)
    call MatLoad(p_mass_matrix, v_infile, ierr)

    call PetscViewerDestroy(v_infile, ierr)

!~       call petsc_get_a_with_b_c(p_mass_matrix,p_Kphonon00_cc(0,0),p_Kphonon00_cc(0,0),mattype_cc)
!~       call MatShift(p_mass_matrix,p_one,ierr) ;  ! call MatZeroEntries(p_sqrt_mass,ierr)
!~       if (l_ionode) then
!~         open(newunit=iunit,file="mass_au.dat",action="read",status="old")
!~         ii=-1
!~         do i1=1,nat3
!~           ii=ii+1
!~           read(iunit,*) p_in(1)
!~           call MatSetValues(p_mass_matrix,1,ii,1,ii,p_in,INSERT_VALUES,ierr)
!~         end do
!~         close(iunit)
!~       end if
!~       call MatAssemblyBegin(p_mass_matrix,MAT_FINAL_ASSEMBLY,ierr)
!~       call MatAssemblyEnd(p_mass_matrix,MAT_FINAL_ASSEMBLY,ierr)
!~       if (inode.eq.0) then
!~         ii=24
!~         call MatGetValues(p_mass_matrix,1,ii,1,ii,p_in,ierr)
!~         write(6,*) "mass ",ii,p_in
!~       end if
!~       call petsc_mat_info(p_invsqrt_mass_matrix,"p_mass_matrix ",ierr)

  end subroutine init_phonons

  subroutine get_phonones()
#include <petsc/finclude/petsc.h>
    use petscmat
    use slepceps
    use kinds
    use petsc_mod
    use ft_mod
    use globals
    use slepc_mod

    implicit none

    Mat :: p_tmp, p_tmp1, p_tmp2, p_tmp3
    integer :: ia1, j, i, ii, ierr, i1, i2, jj, iroot
    PetscInt :: nl1, nl2
    PetscReal :: norm
    PetscScalar, pointer :: p_x(:)

    do i1 = 1, nk_c
      call MatDuplicate(p_Kphonon00_cc(0, 0), MAT_DO_NOT_COPY_VALUES, p_Kphonon00k_cc(i1), ierr)
      call MatCreateVecs(p_Kphonon00k_cc(i1), p_phonon_EW(i1), PETSC_NULL_vec, ierr)
      call fourier_trans(p_Kphonon00_cc, p_Kphonon00k_cc, kp_r, nk_r, dlat_r, ncell_c(1), ncell_c(2), 0)
    end do

    n_ep_modes_k = nat_ecc*3
    do i1 = 1, nk_c
      call petsc_mat_info(p_Kphonon00k_cc(i1), "p_Kphonon00k_cc(i1) ", ierr)
      call MatPtAP(p_Kphonon00k_cc(i1), p_invsqrt_mass, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, p_tmp, ierr)
      call petsc_get_densemat(p_Kphonon00k_cc(i1), p_tmp2, mattype_dense)
      call diag_mat2(p_tmp, PETSC_NULL_MAT, p_phonon_EW(i1), p_tmp2, nat3, -1d0, EPS_LARGEST_REAL)
      call MatGetOwnershipRange(p_tmp, nl1, nl2, ierr)
      call VecGetArrayF90(p_phonon_EW(i1), p_x, ierr)
      ! lets just take complex root, for negative frequencies (energies) this may fuck-up the transport but is good to check for instabilities
      ! this has to be handled somehow
      p_x = zsqrt(p_x) !sign(sqrt(abs((real(p_x,8)))),real(p_x,8))
      write (pstr_out, fmt='(A8,A8,A48,A48)') "k-point", "Mode", "Energy(H)", "Energy(eV)"; call petsc_print_master()

      do j = 0, nprocs - 1
        if (j .eq. inode) then
          i = 0
          do i2 = nl1, nl2 - 1
            i = i + 1
            ii = i2
            jj = 0
            iroot = -1
            if ((real(p_x(i)) .le. 1d-8) .and. (aimag(p_x(i)) .le. 1d-12) .and.&
            &(n_ep_modes_k(i1) .eq. nat3)) n_ep_modes_k(i1) = i2
            write (6, fmt='(2i8,4e24.12E4)') i1, i2, p_x(i), p_x(i)*27.2114d0
!~                 write(6,fmt='(A)') pstr_out
            call flush (6)
          end do
        end if
        call MPI_BARRIER(PETSC_COMM_WORLD, ierr)
      end do

      call MPI_BARRIER(PETSC_COMM_WORLD, ierr)
      call MPI_AllReduce(MPI_IN_PLACE, n_ep_modes_k(i1), 1, MPI_INTEGER, MPI_MIN, PETSC_COMM_WORLD, ierr)

      call MPI_BARRIER(PETSC_COMM_WORLD, ierr)
      write (pstr_out, fmt='(A,i8)') "number of nonzero modes=", n_ep_modes_k(i1); call petsc_print_master()
      call VecRestoreArrayF90(p_phonon_EW(i1), p_x, ierr)
      call MatDestroy(p_tmp, ierr)
      call MatMatMult(p_invsqrt_mass, p_tmp2, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, p_phonon_EV(i1), ierr) ! get mass weighted normal coordinates
!~         call MatConvert(p_phonon_EV(i1),mattype_sparse,MAT_INPLACE_MATRIX,p_phonon_EV(i1),ierr) ! this is not necessary anymore as MatMatMult Dense*AIJ inconcistency has been fixed in petsc 3.13
      call MatDestroy(p_tmp2, ierr)
      call MatHermitianTranspose(p_phonon_EV(i1), MAT_INITIAL_MATRIX, p_tmp2, ierr)
      call MatMatMatMult(p_tmp2, p_mass_matrix, p_phonon_EV(i1), MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, p_tmp3, ierr)
      call MatShift(p_tmp3, p_minus1, ierr)
      call MatNorm(p_tmp3, NORM_FROBENIUS, norm, ierr)
      write (pstr_out, fmt='(A,e24.12)') "norm(A'*M*A-1)=", norm; call petsc_print_master()

      call MatDestroy(p_tmp3, ierr)
      call MatDestroy(p_tmp2, ierr)
      call MatDestroy(p_Kphonon00k_cc(i1), ierr)

!~         call petsc_mat_info(p_ep_lambda_k(1,i1),"p_ep_lambda(i1) ",ierr)

    end do
    n_ep_modes_max = maxval(n_ep_modes_k)

    nullify (p_x)

  end subroutine get_phonones
end module phonon_mod
