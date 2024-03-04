module loe_mod
#include <petsc/finclude/petsc.h>
  use petscmat
  use kinds
  implicit none

! T0=Tr(Gr*GammaR*Ga*GammaL)
! Tec=2*Re[Tr(Gr*GammaR*Ga*GammaL*Gr*Sigma*Gr*Sigma)]
! TecL=Im[Tr(Gr*GammaR*Ga*GammaL*Gr*Sigma*Gr*GammaL*Ga*Sigma)]
! TecR=Im[Tr(Gr*GammaR*Ga*GammaL*Gr*Sigma*Gr*GammaR*Ga*Sigma)]
! Tin=Tr(Gr*GammaR*Ga*Sigma*Ga*GammaL*Gr*Sigma)
! TJL=Re[Tr(Gr*GammaR*Ga*GammaL*Gr*Sigma*Gr*GammaL*Ga*Sigma)]
! TJR=Re[Tr(Gr*GammaR*Ga*GammaL*Gr*Sigma*Gr*GammaR*Ga*Sigma)]
! TII=2*Re[Tr(Gr*GammaR*Ga*GammaL*Gr*Sigma)]

  real(dp) :: bias_ep

  complex(dp) :: T0
  complex(dp), allocatable :: Tec(:), TecL(:), TecR(:), Tin(:), TII(:), TJL(:), TJR(:), &
                              TpiLR(:), TpiLL(:), TpiRR(:), ImPi_r_fac(:)

  Mat :: p_gr_inv, p_gr, p_full_gammal, p_full_gammar

contains
  subroutine loe_run()
    use globals
    use petsc_mod
    use phonon_mod, only: ii_mode, mode_ii_energy
    use integrator_mod
    use misc
    use error_handler
    implicit none

    integer :: ik, imode, i, n, ib, iunit, ierr, junit
    real(dp) :: dd_eps, d_high, d_low, dx, x, bias_low, bias_high, dbias
    complex(dp) :: zz, current_out, neff_test
    complex(dp), allocatable :: current_loe(:), current_elastic_p_loe(:)
    PetscScalar :: dI0_ec_out, dI0_J_out, Iinel_out

    write (pstr_out, fmt='(A)') "init LOE matrices"; call petsc_print_master()
    call init_loe_matrices()

    bias_high = ep_bias(2)
    bias_low = ep_bias(1)
    n = ep_bias_n
    dbias = (bias_high - bias_low)/real(n)
    n = n + 1
    allocate (current_loe(n), current_elastic_p_loe(n), stat=ierr)
    if (ierr .ne. 0) then
      write (errormsg, *) "allocation error in loe_run ", ierr
      call error()
    end if
    current_loe = 0d0
    current_out = 0d0
    current_elastic_p_loe = 0d0

    do ik = nk_c, 1, -1

      write (pstr_out, fmt='(A,i8)') "loe current at k-point for all bias voltages ", ik; call petsc_print_master()

      call init_loe_coefficients(ik)

      write (pstr_out, fmt='(A,i8)') "get current for each mode", ik; call petsc_print_master()
      write (pstr_out, fmt='(A)') "mode "; call petsc_print_master(.false.)

      do imode = n_ep_modes_k(ik), 1, -1
        write (pstr_out, fmt='(i6)') imode; call petsc_print_master(.false.)
        call petsc_vec_getvalue(imode - 1, mode_ii_energy, p_phonon_EW(ik))

        ii_mode = imode
        do ib = 1 + inode, n, nprocs
          bias_ep = bias_low + dbias*real(ib - 1, 8)
!~ integrate  over rho_eq_phonon. does pretty much not differ from taking the limit to rho_eq_phonon->delta(omega) for eta->0
!~               dd_eps=1d-2
!~               d_high=sqrt(1d0/dd_eps)
!~               d_low=0d0

!~               call  loe_integrate(dI0_ec_petsc,d_low,d_high,current_out)
!~               current_loe(ib)=current_loe(ib)+current_out*wkp_l(ik)

!~               call  loe_integrate(Iinel_petsc,d_low,d_high,current_out)
!~               current_loe(ib)=current_loe(ib)+current_out*wkp_l(ik)

!~               dd_eps=1d-4
!~               d_high=zsqrt(2d0*mode_ii_energy/dd_eps+mode_ii_energy*mode_ii_energy+eta_ph_cc*eta_ph_cc*0.25d0)
!~               d_low=-d_high

!~               call  loe_integrate(dI0_J_petsc,d_low,d_high,current_out)
!~               current_loe(ib)=current_loe(ib)+current_out*wkp_l(ik)

          call loe_current(real(mode_ii_energy, 8), dI0_ec_out, dI0_J_out, Iinel_out)
          current_loe(ib) = current_loe(ib) + (dI0_ec_out + dI0_J_out + Iinel_out)*wkp_l(ik)
!~               write(pstr_out,fmt='(2i8,e24.12,i8,99e24.12)') ik,n-ib,bias_ep,imode,real(mode_ii_energy,8),real(current_loe(ib)),&
!~                 real(dI0_ec_out),real(dI0_J_out),real(Iinel_out) ; call petsc_print_master()

        end do ! bias

        call MPI_BARRIER(PETSC_COMM_WORLD, ierr)
      end do ! modes

      write (pstr_out, fmt='(A)') " finished"; call petsc_print_master()

      do ib = 1 + inode, n, nprocs
        bias_ep = bias_low + dbias*real(ib - 1, 8)
        current_elastic_p_loe(ib) = current_loe(ib) + T0*bias_ep*wkp_l(ik)
      end do ! bias

      call destroy_loe_coeffcients()

    end do !kpoit

    current_elastic_p_loe = current_elastic_p_loe/real(nktot, 8)
    current_loe = current_loe/real(nktot, 8)

    if (l_ionode) then
      call MPI_REDUCE(MPI_IN_PLACE, current_elastic_p_loe, size(current_elastic_p_loe), MPI_DOUBLE_COMPLEX, MPI_SUM, 0, PETSC_COMM_WORLD, ierr)
      call MPI_REDUCE(MPI_IN_PLACE, current_loe, size(current_loe), MPI_DOUBLE_COMPLEX, MPI_SUM, 0, PETSC_COMM_WORLD, ierr)
    else
      call MPI_REDUCE(current_elastic_p_loe, current_elastic_p_loe, size(current_elastic_p_loe), MPI_DOUBLE_COMPLEX, MPI_SUM, 0, PETSC_COMM_WORLD, ierr)
      call MPI_REDUCE(current_loe, current_loe, size(current_elastic_p_loe), MPI_DOUBLE_COMPLEX, MPI_SUM, 0, PETSC_COMM_WORLD, ierr)
    end if

    if (l_ionode) then
      open (newunit=iunit, file="current_loe.dat", action="write", status="replace")
      do ib = 1, n
        bias_ep = bias_low + dbias*real(ib - 1, 8)
        write (unit=iunit, fmt='(4es45.24e5)') bias_ep, real(current_elastic_p_loe(ib), 8), real(current_loe(ib), 8), real(current_elastic_p_loe(ib) - current_loe(ib), 8)
      end do
      close (iunit)

      do ik = 1, nk_c
        open (newunit=iunit, file="modes_"//trim(i2s(ik))//".dat", action="write", status="replace")
        do imode = 1, n_ep_modes_k(ik)
          call petsc_vec_getvalue(imode - 1, mode_ii_energy, p_phonon_EW(ik), .false.)
          write (iunit, fmt='(2e24.12)') real(mode_ii_energy), 0d0
          write (iunit, fmt='(2e24.12)') real(mode_ii_energy), -1d0
          write (iunit, fmt='(2e24.12)') real(mode_ii_energy), 1d0
          write (iunit, fmt='(2e24.12)') real(mode_ii_energy), 0d0
        end do
        close (iunit)
      end do

    end if

    call destroy_loe_matrices()

  end subroutine loe_run

  subroutine init_loe_matrices()
#include <petsc/finclude/petsc.h>
    use petscmat
    use petsc_mod
    use petsc_wrapper
    use kinds
    use globals
    use integrator_mod

    implicit none

    integer :: ierr

! initialize necessary matricies
    call petsc_get_densemat(p_h00k_r(1), p_grrr, mattype_surf)     ! elemental
    call petsc_get_densemat(p_h00k_l(1), p_gllr, mattype_surf)     ! elemental
    call petsc_get_densemat(p_h00k_r(1), p_sigmarr, mattype_dense) ! dense
    call petsc_get_densemat(p_h00k_r(1), p_gammar, mattype_dense)  ! dense
    call petsc_get_densemat(p_h00k_l(1), p_sigmalr, mattype_dense) ! dense
    call petsc_get_densemat(p_h00k_l(1), p_gammal, mattype_dense)  ! dense

    call petsc_get_a_with_b_c(p_full_gammal, p_h00k_cc(1), p_h00k_cc(1), mattype_sparse)
    call petsc_get_a_with_b_c(p_full_gammar, p_h00k_cc(1), p_h00k_cc(1), mattype_sparse)

    call petsc_alloc_mat_block(p_full_gammal, 0, nmu_l - 1, 0, nmu_l - 1)
    call petsc_alloc_mat_block(p_full_gammar, nmu_c - nmu_r, nmu_c - 1, nmu_c - nmu_r, nmu_c - 1)

  end subroutine init_loe_matrices

  subroutine destroy_loe_matrices()
#include <petsc/finclude/petsc.h>
    use petscmat
    use petsc_mod
    use petsc_wrapper
    use kinds
    use globals
    use integrator_mod

    implicit none

    integer :: ierr

! initialize necessary matricies
    call MatDestroy(p_grrr, ierr)
    call MatDestroy(p_gllr, ierr)
    call MatDestroy(p_sigmarr, ierr)
    call MatDestroy(p_gammar, ierr)
    call MatDestroy(p_sigmalr, ierr)
    call MatDestroy(p_gammal, ierr)

    call MatDestroy(p_full_gammal, ierr)
    call MatDestroy(p_full_gammar, ierr)

    call MatDestroy(p_full_gammal, ierr)
    call MatDestroy(p_full_gammar, ierr)

  end subroutine destroy_loe_matrices

  subroutine destroy_loe_coeffcients()
#include <petsc/finclude/petsc.h>
    use petscmat
    implicit none

    integer :: ierr

    deallocate (Tec)
    deallocate (TecL)
    deallocate (TecR)
    deallocate (Tin)
    deallocate (TII)
    deallocate (TJL)
    deallocate (TJR)
    deallocate (TpiLR)
    deallocate (TpiRR)
    deallocate (TpiLL)
    deallocate (ImPi_r_fac)

  end subroutine destroy_loe_coeffcients

  subroutine init_loe_coefficients(ik)
#include <petsc/finclude/petsc.h>
    use petscmat
    use petsc_mod
    use petsc_wrapper
    use kinds
    use globals
    use integrator_mod
    use error_handler

    implicit none

    integer :: ik

    Mat :: p_gr, p_GrGammaRGa, p_Ga, p_tmp1, p_tmp2, p_tmp3, p_GrGammaRGaGammaLGr, &
      p_GrGammaRGaGammaLGrSigmaGr, p_gaSigma

    integer :: ierr, ii(1), jj(1), imode
    complex(dp) :: ef, pp(1), pi_pm, pi_r
    PetscReal :: p_real
    PetscScalar :: p_tr

    write (pstr_out, fmt='(A,i8)') "init leo coefficients ", ik; call petsc_print_master()

    iik = ik

    allocate (Tec(n_ep_modes_k(ik)), TecL(n_ep_modes_k(ik)), TecR(n_ep_modes_k(ik)), &
              Tin(n_ep_modes_k(ik)), TII(n_ep_modes_k(ik)), TJL(n_ep_modes_k(ik)), &
              TJR(n_ep_modes_k(ik)), TpiLR(n_ep_modes_k(ik)), TpiRR(n_ep_modes_k(ik)), &
              TpiLL(n_ep_modes_k(ik)), ImPi_r_fac(n_ep_modes_k(ik)), stat=ierr)
    if (ierr .ne. 0) then
      write (errormsg, *) "allocation error in loe_run ", ierr
      call error()
    end if

    ef = 0.5d0*(ef_l + ef_r) ! probably ok for hetrojunctions

    mul = ef_l!+vl
    mur = ef_r!+vr

    nsim_inv = nmu_c

    ldouble_contour = .false.

    call init_gr(ef, p_gr_inv) ! get Gr^-1

    call petsc_add_sub_B_to_A(p_gammal, p_full_gammal, 0, 0, p_one, INSERT_VALUES, PETSC_TRUE)
    call petsc_add_sub_B_to_A(p_gammar, p_full_gammar, nmu_c - nmu_r, nmu_c - nmu_r, p_one, INSERT_VALUES, PETSC_TRUE)

!~ ! debug ----
!~       call petsc_mat_info(p_gammal,"p_gammal ",ierr)
!~       call petsc_mat_info(p_gammar,"p_gammar ",ierr)
!~       call petsc_mat_info(p_full_gammal,"p_full_gammal ",ierr)
!~       call petsc_mat_info(p_full_gammar,"p_full_gammar ",ierr)
!~ ! ----------

    call petsc_get_densemat(p_h00k_cc(1), p_gr, mattype_dense) ! dense
    call petsc_invert(p_gr_inv, p_gr, matsolvertype_cc, mattype_dense)
!~       call MatConvert(p_gr,mattype_sparse,MAT_INPLACE_MATRIX,p_gr,ierr) ! something has to be converted either GammaL,R to dense (elemental) or Gr to sparse
    call MatHermitianTranspose(p_gr, MAT_INITIAL_MATRIX, p_ga, ierr)
    call MatMatMatMult(p_gr, p_full_gammar, p_ga, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, p_GrGammaRGa, ierr)
    call MatMatMatMult(p_GrGammaRGa, p_full_gammal, p_gr, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, &
                       p_GrGammaRGaGammaLGr, ierr)

    call MatMatMult(p_GrGammaRGa, p_full_gammal, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, &
                    p_tmp3, ierr) !
    call MatGetTrace(p_tmp3, p_tr, ierr)
    T0 = p_tr
    write (pstr_out, fmt='(A,2e24.12)') "T0", T0; call petsc_print_master()
    call MatDestroy(p_tmp3, ierr)

! Tec=2*Re[Tr(Gr*GammaR*Ga*GammaL*Gr*Sigma*Gr*Sigma)]
! TecL= Im[Tr(Gr*GammaR*Ga*GammaL*Gr*Sigma*Gr*GammaL*Ga*Sigma)]
! TecR= Im[Tr(Gr*GammaR*Ga*GammaL*Gr*Sigma*Gr*GammaR*Ga*Sigma)]
! TJL=  Re[Tr(Gr*GammaR*Ga*GammaL*Gr*Sigma*Gr*GammaL*Ga*Sigma)]
! TJR=  Re[Tr(Gr*GammaR*Ga*GammaL*Gr*Sigma*Gr*GammaR*Ga*Sigma)]

! TpiLR=   Tr(Gr*GammaR*Ga*Sigma*Gr*GammaL*Ga*Sigma)
! TpiLL=   Tr(Gr*GammaL*Ga*Sigma*Gr*GammaL*Ga*Sigma)
! TpiRR=   Tr(Gr*GammaR*Ga*Sigma*Gr*GammaR*Ga*Sigma)
!
! TII=2*Re[Tr(Gr*GammaR*Ga*GammaL*Gr*Sigma)]
!
! Tin=     Tr(Gr*GammaR*Ga*Sigma*Ga*GammaL*Gr*Sigma)
!
    write (pstr_out, fmt='(A)') "mode "; call petsc_print_master(.false.)
    do imode = n_ep_modes_k(ik), 1, -1
      write (pstr_out, fmt='(i6)') imode; call petsc_print_master(.false.)
      call MatMatMatMult(p_GrGammaRGaGammaLGr, p_ep_lambda_k(imode, ik), p_gr, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, &
                         p_GrGammaRGaGammaLGrSigmaGr, ierr)

      call MatMatMult(p_GrGammaRGaGammaLGrSigmaGr, p_ep_lambda_k(imode, ik), MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, p_tmp1, ierr) !Tec
      call MatGetTrace(p_tmp1, p_tr, ierr)
      Tec(imode) = 2d0*real(p_tr, 8)
      call MatDestroy(p_tmp1, ierr)

      call MatMatMult(p_ga, p_ep_lambda_k(imode, ik), MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, p_gaSigma, ierr) !Ga*Sigma

      call MatMatMatMult(p_GrGammaRGaGammaLGrSigmaGr, p_full_gammal, p_gaSigma, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, &
                         p_tmp1, ierr) !TecL and TJL (in WBL)
      call MatGetTrace(p_tmp1, p_tr, ierr)
      TecL(imode) = aimag(p_tr)
      TJL(imode) = real(p_tr, 8)
      call MatDestroy(p_tmp1, ierr)

      call MatMatMatMult(p_GrGammaRGaGammaLGrSigmaGr, p_full_gammar, p_gaSigma, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, &
                         p_tmp1, ierr) !TecR and TJR (in WBL)
      call MatGetTrace(p_tmp1, p_tr, ierr)
      TecR(imode) = aimag(p_tr)
      TJR(imode) = real(p_tr, 8)
      call MatDestroy(p_tmp1, ierr)

      call MatMatMult(p_GrGammaRGaGammaLGr, p_ep_lambda_k(imode, ik), MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, &
                      p_tmp1, ierr) ! TII
      call MatGetTrace(p_tmp1, p_tr, ierr)
      TII(imode) = 2d0*real(p_tr, 8)
      call MatDestroy(p_tmp1, ierr)

      call MatMatMatMult(p_GrGammaRGa, p_ep_lambda_k(imode, ik), p_ga, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, &
                         p_tmp1, ierr) ! Gr*GammaR*Ga*Sigma*Ga
      call MatMatMatMult(p_tmp1, p_full_gammal, p_gr, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, &
                         p_tmp2, ierr) ! Gr*GammaR*Ga*Sigma*Ga*GammaL*Gr
      call MatDestroy(p_tmp1, ierr)
      call MatMatMult(p_tmp2, p_ep_lambda_k(imode, ik), MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, p_tmp1, ierr)
      call MatGetTrace(p_tmp1, p_tr, ierr)
      Tin(imode) = p_tr
      call MatDestroy(p_tmp1, ierr)
      call MatDestroy(p_tmp2, ierr)

! TpiLR=   Tr(Gr*GammaR*Ga*Sigma*Gr*GammaL*Ga*Sigma)
! TpiRR=   Tr(Gr*GammaR*Ga*Sigma*Gr*GammaR*Ga*Sigma)
! TpiLL=   Tr(Gr*GammaL*Ga*Sigma*Gr*GammaL*Ga*Sigma)

      call MatMatMatMult(p_Gr, p_full_gammar, p_gaSigma, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, &
                         p_tmp1, ierr) ! Gr*GammaR*Ga*Sigma
      call MatMatMatMult(p_tmp1, p_gr, p_full_gammal, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, &
                         p_tmp2, ierr) ! Gr*GammaR*Ga*Sigma*Gr*GammaL
      call MatMatMult(p_tmp2, p_gaSigma, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, &
                      p_tmp3, ierr) ! Gr*GammaR*Ga*Sigma*Gr*GammaL*Ga*Sigma
      call MatGetTrace(p_tmp3, p_tr, ierr)
      TpiLR(imode) = 0.5d0/pi*p_tr

      call MatDestroy(p_tmp3, ierr)
      call MatDestroy(p_tmp2, ierr)

      call MatMatMatMult(p_tmp1, p_gr, p_full_gammar, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, &
                         p_tmp2, ierr) ! Gr*GammaR*Ga*Sigma*Gr*GammaR
      call MatMatMult(p_tmp2, p_gaSigma, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, &
                      p_tmp3, ierr) ! Gr*GammaR*Ga*Sigma*Gr*GammaR*Ga*Sigma
      call MatGetTrace(p_tmp3, p_tr, ierr)
      TpiRR(imode) = 0.5d0/pi*p_tr

      call MatDestroy(p_tmp3, ierr)
      call MatDestroy(p_tmp2, ierr)
      call MatDestroy(p_tmp1, ierr)

      call MatMatMatMult(p_Gr, p_full_gammal, p_gaSigma, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, &
                         p_tmp1, ierr) ! Gr*GammaL*Ga*Sigma
      call MatMatMatMult(p_tmp1, p_gr, p_full_gammal, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, &
                         p_tmp2, ierr) ! Gr*GammaL*Ga*Sigma*Gr*GammaL
      call MatMatMult(p_tmp2, p_gaSigma, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, &
                      p_tmp3, ierr) ! Gr*GammaL*Ga*Sigma*Gr*GammaL*Ga*Sigma
      call MatGetTrace(p_tmp3, p_tr, ierr)
      TpiLL(imode) = 0.5d0/pi*p_tr

!~         write(pstr_out,fmt='(A,i8,2e24.12)') "TpiLL",imode,TpiLL(imode) ; call petsc_print_master()
!~         write(pstr_out,fmt='(A,i8,2e24.12)') "TpiRR",imode,TpiRR(imode) ; call petsc_print_master()
!~         write(pstr_out,fmt='(A,i8,2e24.12)') "TpiLR",imode,TpiLR(imode) ; call petsc_print_master()

      call MatDestroy(p_tmp3, ierr)
      call MatDestroy(p_tmp2, ierr)
      call MatDestroy(p_tmp1, ierr)
      call MatDestroy(p_gaSigma, ierr)
      call MatDestroy(p_GrGammaRGaGammaLGrSigmaGr, ierr)

      call MatDuplicate(p_gr, MAT_COPY_VALUES, p_tmp1, ierr)
      call MatImaginaryPart(p_tmp1, ierr)
      call MatMatMult(p_ep_lambda_k(imode, ik), p_tmp1, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, &
                      p_tmp2, ierr)
      call MatMatMult(p_tmp2, p_tmp2, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, &
                      p_tmp3, ierr)
      call MatGetTrace(p_tmp3, p_tr, ierr)
      ImPi_r_fac(imode) = -1d0/pi*p_tr

      call MatDestroy(p_tmp3, ierr)
      call MatDestroy(p_tmp2, ierr)
      call MatDestroy(p_tmp1, ierr)

!~         write(pstr_out,fmt='(A,i8,2e24.12)') "Tec(imode)",imode,Tec(imode) ; call petsc_print_master()
!~         write(pstr_out,fmt='(A,i8,2e24.12)') "TecL(imode)",imode,TecL(imode) ; call petsc_print_master()
!~         write(pstr_out,fmt='(A,i8,2e24.12)') "TecR(imode)",imode,TecR(imode) ; call petsc_print_master()
!~         write(pstr_out,fmt='(A,i8,2e24.12)') "Tin(imode)",imode,Tin(imode) ; call petsc_print_master()
!~         write(pstr_out,fmt='(A,i8,2e24.12)') "TII(imode)",imode,TII(imode) ; call petsc_print_master()
!~         write(pstr_out,fmt='(A,i8,2e24.12)') "TJL(imode)",imode,TJL(imode) ; call petsc_print_master()
!~         write(pstr_out,fmt='(A,i8,2e24.12)') "TJR(imode)",imode,TJR(imode) ; call petsc_print_master()
!~         write(pstr_out,fmt='(A,i8,2e24.12)') "TpiLR(imode)",imode,TpiLR(imode) ; call petsc_print_master()
!~         write(pstr_out,fmt='(A,i8,2e24.12)') "TpiRR(imode)",imode,TpiRR(imode) ; call petsc_print_master()
!~         write(pstr_out,fmt='(A,i8,2e24.12)') "TpiLL(imode)",imode,TpiLL(imode) ; call petsc_print_master()

    end do

    write (pstr_out, fmt='(A)') " finished"; call petsc_print_master()

    call MatDestroy(p_GrGammaRGaGammaLGr, ierr)
    call MatDestroy(p_ga, ierr)
    call MatDestroy(p_gr, ierr)
    call MatDestroy(p_GrGammaRGa, ierr)
    call MatDestroy(p_gr_inv, ierr)

  end subroutine init_loe_coefficients

  function ImPi_r_alpha(x)
    use globals
    use kinds
    use phonon_mod
    implicit none

    PetscScalar :: ImPi_r_alpha
    real(dp) :: x

    ImPi_r_alpha = ImPi_r_fac(ii_mode)*x

  end function ImPi_r_alpha

  function ImPi_pm_alpha(x)
    use globals
    use kinds
    use phonon_mod
    use integrator_mod, only: bose
    use petsc_mod
    implicit none

    PetscScalar :: ImPi_pm_alpha
    real(dp) :: x

    complex(dp) :: x_p_u, x_m_u, xx

    xx = x
    x_p_u = x + bias_ep
    x_m_u = x - bias_ep

    ImPi_pm_alpha = 0.5d0/pi*(TpiLR(ii_mode)*(x_p_u*bose(x_p_u, temperature_ph) + &
                                              (x_m_u)*bose(x_m_u, temperature_ph)) + (TpiLL(ii_mode) + TpiRR(ii_mode))*xx*bose(xx, temperature_ph))
!~         ImPi_pm_alpha=0.5d0/pi*(1d0*(x_p_u*bose(x_p_u,temperature_ph)+&
!~                       (x_m_u)*bose(x_m_u,temperature_ph))+(1d0+1d0)*xx*bose(xx,temperature_ph))

!~         write(pstr_out,fmt='(A,4e24.12)') "ImPi_pm_alpha, TpiLR(ii_mode) ",TpiLR(ii_mode) ; call petsc_print_master()
!~         write(pstr_out,fmt='(A,4e24.12)') "ImPi_pm_alpha, TpiLL(ii_mode) ",TpiLL(ii_mode) ; call petsc_print_master()
!~         write(pstr_out,fmt='(A,4e24.12)') "ImPi_pm_alpha, TpiRR(ii_mode) ",TpiRR(ii_mode) ; call petsc_print_master()
!~         write(pstr_out,fmt='(A,4e24.12)') "ImPi_pm_alpha, x_p_u ",x_p_u ; call petsc_print_master()
!~         write(pstr_out,fmt='(A,4e24.12)') "ImPi_pm_alpha, x_m_u ",x_m_u ; call petsc_print_master()
!~         write(pstr_out,fmt='(A,4e24.12)') "ImPi_pm_alpha, bose(x_p_u,temperature_ph) ",bose(x_p_u,temperature_ph) ; call petsc_print_master()
!~         write(pstr_out,fmt='(A,4e24.12)') "ImPi_pm_alpha, bose(x_m_u,temperature_ph) ",bose(x_m_u,temperature_ph) ; call petsc_print_master()

  end function ImPi_pm_alpha

  function Nneq_alpha(x)
    use globals
    use kinds
    use integrator_mod, only: bose
    use phonon_mod
    use petsc_mod
    implicit none

    PetscScalar :: Nneq_alpha
    real(dp) :: x
    complex(dp) :: x_p_u, x_m_u

    complex(dp) :: xx

    xx = x
    x_p_u = x + bias_ep
    Nneq_alpha = -0.5d0*(ImPi_pm_alpha(x) + bose(xx, temperature_ph)*eta_ph_cc*x/mode_ii_energy)/(ImPi_r_alpha(x) - eta_ph_cc*x/mode_ii_energy*0.5d0)

!~         write(pstr_out,fmt='(A,4e24.12)') "Nneq_alpha, x ",xx,x ; call petsc_print_master()
!~         write(pstr_out,fmt='(A,4e24.12)') "Nneq_alpha, mode_ii_energy ",mode_ii_energy ; call petsc_print_master()
!~         write(pstr_out,fmt='(A,4e24.12)') "Nneq_alpha, eta_ph_cc ",eta_ph_cc ; call petsc_print_master()
!~         write(pstr_out,fmt='(A,4e24.12)') "Nneq_alpha, ImPi_pm_alpha(x) ",ImPi_pm_alpha(x) ; call petsc_print_master()
!~         write(pstr_out,fmt='(A,4e24.12)') "Nneq_alpha, ImPi_r_alpha(x) ",ImPi_r_alpha(x) ; call petsc_print_master()
!~         write(pstr_out,fmt='(A,4e24.12)') "Nneq_alpha, bose(xx,temperature_ph)* ",bose(xx,temperature_ph) ; call petsc_print_master()

  end function Nneq_alpha

  function dI0_ec(x)
    use globals
    use kinds
    use petsc_mod
    use integrator_mod, only: bose
    use phonon_mod
    implicit none

    PetscScalar :: dI0_ec
    real(dp) :: x

    complex(dp) :: x_p_u, x_m_u, xx

    xx = x
    x_p_u = x + bias_ep
    x_m_u = x - bias_ep

!~       dI0_ec=rho_eq_alpha(x)*(Tec(ii_mode)*(2d0*Nneq_alpha(x)+1d0)*bias_ep+&
!~              (TecL(ii_mode)+TecR(ii_mode))*(x_m_u*bose(x_m_u,temperature_ph)-x_p_u*bose(x_p_u,temperature_ph)-bias_ep))
    dI0_ec = (Tec(ii_mode)*(2d0*Nneq_alpha(x) + 1d0)*bias_ep + &
              (TecL(ii_mode) + TecR(ii_mode))*(x_m_u*bose(x_m_u, temperature_ph) - x_p_u*bose(x_p_u, temperature_ph) - bias_ep))

  end function dI0_ec

  function dI0_ec_petsc(x)
    use globals
    use kinds
    use petsc_wrapper, only: petsc_matassemble
    use petsc_mod
    implicit none

    real(dp) :: x

    real(dp) :: z
    integer :: ii(1), ierr, dI0_ec_petsc
    PetscScalar :: pp(1)

    dI0_ec_petsc = -1

    if (inode .eq. 0) then
      z = x
      pp(1) = dI0_ec(z)
      ii = 0
      call MatSetValues(p_tmpcc1, 1, ii, 1, ii, pp, INSERT_VALUES, ierr)
    end if

    call MatAssemblyBegin(p_tmpcc1, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(p_tmpcc1, MAT_FINAL_ASSEMBLY, ierr)

    call petsc_matassemble(p_tmpcc1)

    dI0_ec_petsc = 0

  end function dI0_ec_petsc

  function dI0_J(x)
    use globals
    use kinds
    use petsc_mod
    use integrator_mod, only: bose
    use phonon_mod
    implicit none

    PetscScalar :: dI0_J
    real(dp) :: x

    complex(dp) :: x_p_u, x_m_u, xx

    xx = x
    x_p_u = x + bias_ep
    x_m_u = x - bias_ep

    dI0_J = real(dr_alpha(x))*(xx*bose(xx, temperature_ph) - (x_p_u)*bose(x_p_u, temperature_ph))

!~       write(pstr_out,fmt='(A,4e24.12)') "dI0_ec, x ",xx,x ; call petsc_print_master()
!~       write(pstr_out,fmt='(A,2e24.12)') "dI0_ec, temperature_ph ",temperature_ph ; call petsc_print_master()
!~       write(pstr_out,fmt='(A,2e24.12)') "dI0_ec, x_p_u ",x_p_u ; call petsc_print_master()
!~       write(pstr_out,fmt='(A,2e24.12)') "dI0_ec, x_m_u ",x_m_u ; call petsc_print_master()
!~       write(pstr_out,fmt='(A,2e24.12)') "dI0_ec, rho_eq_alpha(x) ",rho_eq_alpha(x) ; call petsc_print_master()
!~       write(pstr_out,fmt='(A,2e24.12)') "dI0_ec, Tec(ii_mode) ",Tec(ii_mode) ; call petsc_print_master()
!~       write(pstr_out,fmt='(A,2e24.12)') "dI0_ec, Nneq_alpha(x) ",Nneq_alpha(x) ; call petsc_print_master()
!~       write(pstr_out,fmt='(A,2e24.12)') "dI0_ec, TecL(ii_mode) ",TecL(ii_mode) ; call petsc_print_master()
!~       write(pstr_out,fmt='(A,2e24.12)') "dI0_ec, TecR(ii_mode) ",TecR(ii_mode) ; call petsc_print_master()
!~       write(pstr_out,fmt='(A,2e24.12)') "dI0_ec, bose(x_m_u,temperature_ph) ",bose(x_m_u,temperature_ph) ; call petsc_print_master()
    dI0_J = -1d0/pi*(TJR(ii_mode) - TJL(ii_mode))*dI0_J

  end function dI0_J

  function dI0_J_petsc(x)
    use globals
    use kinds
    use petsc_wrapper, only: petsc_matassemble
    use petsc_mod
    implicit none

    real(dp) :: x

    real(dp) :: z
    integer :: ii(1), ierr, dI0_J_petsc
    PetscScalar :: pp(1)

    dI0_J_petsc = -1

    if (inode .eq. 0) then
      z = x
      pp(1) = dI0_J(z)
      ii = 0
      call MatSetValues(p_tmpcc1, 1, ii, 1, ii, pp, INSERT_VALUES, ierr)
    end if

    call MatAssemblyBegin(p_tmpcc1, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(p_tmpcc1, MAT_FINAL_ASSEMBLY, ierr)

    call petsc_matassemble(p_tmpcc1)

    dI0_J_petsc = 0

  end function dI0_J_petsc

  function dI0_J_scalar(x)
    use globals
    use kinds
    use petsc_wrapper, only: petsc_matassemble
    use petsc_mod
    use integrator_scalar
    implicit none

    real(dp) :: x

    real(dp) :: z
    integer :: ii(1), ierr, dI0_J_scalar
    PetscScalar :: pp(1)

    zint_scalar_output = dI0_J(x)
    dI0_J_scalar = 1

  end function dI0_J_scalar

  function Iinel(x)
    use globals
    use kinds
    use petsc_mod
    use integrator_mod, only: bose
    use phonon_mod
    implicit none

    PetscScalar :: Iinel
    real(dp) :: x

    complex(dp) :: x_p_u, x_m_u, xx

    xx = x
    x_p_u = x + bias_ep
    x_m_u = x - bias_ep

!~       Iinel=rho_eq_alpha(x)*(2d0*Nneq_alpha(x)*bias_ep+x_m_u*bose(x_m_u,temperature_ph)-x_p_u*bose(x_p_u,temperature_ph))
    Iinel = (2d0*Nneq_alpha(x)*bias_ep + x_m_u*bose(x_m_u, temperature_ph) - x_p_u*bose(x_p_u, temperature_ph))

    Iinel = Tin(ii_mode)*Iinel

  end function Iinel

  function Iinel_petsc(x)
    use globals
    use kinds
    use petsc_wrapper, only: petsc_matassemble
    use petsc_mod
    implicit none

    real(dp) :: x

    real(dp) :: z
    integer :: ii(1), ierr, Iinel_petsc
    PetscScalar :: pp(1)

    Iinel_petsc = -1

    if (inode .eq. 0) then
      z = x
      pp(1) = Iinel(z)
      ii = 0
      call MatSetValues(p_tmpcc1, 1, ii, 1, ii, pp, INSERT_VALUES, ierr)
    end if

    call MatAssemblyBegin(p_tmpcc1, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(p_tmpcc1, MAT_FINAL_ASSEMBLY, ierr)

    call petsc_matassemble(p_tmpcc1)

    Iinel_petsc = 0

  end function Iinel_petsc

  subroutine loe_current(energy, dI0_ec_out, dI0_J_out, Iinel_out)
    use globals
    use kinds
    use phonon_mod
    implicit none

    PetscScalar :: dI0_ec_out, dI0_J_out, Iinel_out
    real(dp) :: energy
    real(dp) :: dd_eps, d_high, d_low

    dI0_ec_out = 0d0
    dI0_J_out = 0d0
    Iinel_out = 0d0

    dI0_ec_out = dI0_ec(energy)

    dd_eps = 1d-4
    d_high = zsqrt(2d0*energy/dd_eps + mode_ii_energy*energy + eta_ph_cc*eta_ph_cc*0.25d0)
    d_low = -d_high
!~       call loe_integrate(dI0_J_petsc,d_low,d_high,dI0_J_out)
    call loe_integrate_scalar(dI0_J_scalar, d_low, d_high, dI0_J_out)
    Iinel_out = Iinel(energy)

  end subroutine loe_current

!~ this uses the integrator for PETSC matrices with a 1x1 matrix to
!~ emulate a scalar quantity, this has of course a huge overhead and
!~ cannot simply parallelized (in bias points or modes for example)
  subroutine loe_integrate(f, d_low, d_high, zout, what)
#include <petsc/finclude/petsc.h>
    use petscmat
    use petsc_mod
    use petsc_wrapper
    use kinds
    use globals
    use integrator_mod
    use phonon_mod, only: mode_ii_energy
    implicit none

    integer, external :: f
    real(dp) :: d_low, d_high
    complex(dp) :: zout
    character(*), optional :: what

    real(dp) :: err_eq_int, eps
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

!~       call petsc_mat_info(p_tmpcc1,"p_tmpcc1 ",ierr)

    eps = 1d-5
    nint_order = 32
    l_output_progress = .false.
    call adaptive_int3(f, d_low, d_high, p_d_tmp1(1), p_d_tmp2(1), maxsub*10, eps, err_eq_int)
    call petsc_mat_getvalue(0, 0, zout, p_d_tmp1(1), 1, PETSC_COMM_WORLD)

    if (present(what)) then
      write (pstr_out, fmt='(A,X,4e24.12)') trim(what), zout, d_low, d_high; call petsc_print_master()
    end if
    call MatDestroy(p_tmpcc1, ierr)
    call MatDestroy(p_d_tmp1(1), ierr)
    call MatDestroy(p_d_tmp2(1), ierr)

  end subroutine loe_integrate

  subroutine loe_integrate_scalar(f, d_low, d_high, zout, what)
    use petsc_mod
    use kinds
    use globals
    use integrator_scalar
    use phonon_mod, only: mode_ii_energy
    implicit none

    integer, external :: f
    real(dp) :: d_low, d_high
    complex(dp) :: zout
    character(*), optional :: what

    real(dp) :: err_eq_int, eps
    PetscScalar :: p_zint1, p_zint2, pp, dx, x1, x2, x
    complex(dp) :: p_d_tmp1, p_d_tmp2
    integer :: ierr, ii(1), jj(1), ik(2), i, n

    eps = 1d-5
    nint_order = 32
    l_output_progress = .false.
    call adaptive_int3_scalar(f, d_low, d_high, p_d_tmp1, p_d_tmp2, maxsub*10, eps, err_eq_int)
    zout = p_d_tmp1

    if (present(what)) then
      write (pstr_out, fmt='(A,X,4e24.12)') trim(what), zout, d_low, d_high; call petsc_print_master()
    end if

  end subroutine loe_integrate_scalar

end module loe_mod
