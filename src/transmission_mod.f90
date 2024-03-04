module transmission_mod
#include <petsc/finclude/petsc.h>
  use petscmat
  implicit none
  
  Mat :: G13, G31
  

contains

!   subroutine get_transmission(x,trans,dosc,dosl,dosr)
! We want Tr[Gr*GammaR*Ga*GammaL] -> Tr[Gr13*GammaR33*Gr13'*GammaL11]=
! Tr[Gr13*GammaR33*(GammaL11*G13)']=Tr[A*B'].
! 1. get G=(G13 G23 G33)
! 2. A=G13*GammaR33
! 3. B=GammaL11*G13
! 4. use petsc_matmatconj_restricted to calculate A*B' add to p_trace (only diagonal allocated)
! 5. MatTrace(p_trace)
! ------
! get transmission eigenvalues from  diagonalizing GammaR^1/2*Ga*GammaL*Gr*GammaR^1/2 ->
! GammaR33^1/2*Gr13'*GammaL11*Gr13*GammaR33^1/2
! 1.  get GammaR33^1/2
! 2.  A=G13*GammaR33^1/2
! 3.  C=(G13*GammaR33^1/2)'=GammaR33^1/2*Gr13'=A'
! 4.  T=C*GammaL11*A=A'*GammaL11*A
! 5.  Diag(T)
! !!!! Gamma matricies are covariant thus writing Gamma=Gamma^(1/2)*Gamma^(1/2) should be Gamma^(1/2)*S^-1*Gamma^(1/2)
! can be easily checked by writing down as operators instead of matricies.
  subroutine get_transmission(x, trans, t_i, dosc, dosl, dosr)
#include <slepc/finclude/slepceps.h>
    use slepceps
    use petscmat
    use petsc_mod
    use petsc_wrapper
    use globals
    use integrator_mod
    use slepc_mod
    use eigenchannel_mod

    real(dp), intent(in) :: x
    complex(dp), intent(out) :: trans
    Vec, optional :: t_i(:)
    real(dp), intent(out), optional :: dosc, dosl, dosr

    logical :: ldosl, ldosr, ldosc

    Mat :: f, p_gammal_block, p_gammar_block, p_x, b, p_g_block, p_rhs, &
      p_g13, p_g13gammar, p_g31conjgammal, p_gammalg13, p_tmp, p_sqrt_gammar, &
      p_tt, p_tmp2, p_tmp3, p_mone

    Mat, pointer :: p_gr_inv

    complex(dp) :: z

    integer :: ierr, i
    real(dp) :: vbb
    PetscScalar :: p_scal
    PetscScalar, pointer :: xx_v(:)
    PetscReal :: norm

    integer(8) :: counti, count_rate, countf

    ldosc = present(dosc)
    ldosl = present(dosl)
    ldosr = present(dosr)

    z = x

    p_gr_inv => p_invGr ! see eval_gr

    call init_gr(z, p_gr_inv)

    call petsc_get_a_with_b_c(p_g_block, p_gr_inv, p_gammar, mattype_dense)

    solver_mode = 1 ! for now assume solver mode 1
    if (solver_mode .eq. 2) then
      call petsc_get_a_with_b_c(p_rhs, p_gr_inv, p_gammar, mattype_dense)
      call petsc_one_mat(p_rhs, nmu_c - nmu_r + 1 - 1, nmu_c - 1) ! prepare RHS nmu_c,nmu_r ( .., 0 0 .. 1 0, 0 0 .. 0 1 )   C style 1=0
      call petsc_call_solver(p_gr_inv, p_rhs, p_g_block, matsolvertype_cc, solver_mode)
    else if ((solver_mode .eq. 1).or.(solver_mode .eq. 3).or.(solver_mode .eq. 4)) then
      call petsc_get_a_with_b_c(p_rhs, p_gr_inv, p_gammar, mattype_sparse)
      call petsc_one_mat(p_rhs, nmu_c - nmu_r + 1 - 1, nmu_c - 1,2) ! prepare RHS nmu_c,nmu_r ( .., 0 0 .. 1 0, 0 0 .. 0 1 )   C style 1=0
      call petsc_call_solver(p_gr_inv, p_rhs, p_g_block, matsolvertype_cc, solver_mode)
    end if

    call MatDestroy(p_rhs, ierr)

!~     call MatDestroy(p_gr_inv, ierr)

    call petsc_split_matrix(p_g_block, p_tmp, 1, nmu_l, 1, nmu_r, MAT_INITIAL_MATRIX)

    call MatDestroy(p_g_block, ierr)
    call MatDuplicate(p_gammar, MAT_SHARE_NONZERO_PATTERN, p_g13, ierr)
    if (lget_ti) then
      call MatDuplicate(p_gammar, MAT_SHARE_NONZERO_PATTERN, p_sqrt_gammar, ierr)
    end if

    call petsc_aXpY(p_g13, p_tmp, p_one, petsc_false, .false.)

!~     call MatCopy(p_tmp,p_g13,DIFFERENT_NONZERO_PATTERN,ierr)  ! workaround: MatMatMult needs p_g13 local structure to be compatible with p_gammaL/R

    call MatDestroy(p_tmp, ierr)

!~     if (lget_ti.and.1.eq.2) then
    if (lget_ti) then
      call calc_eigenchannel_wf(t_i(iik)) 
    end if
    if (lget_ti.and.1.eq.2) then
!~   if (lget_ti) then
! transmission channels t_i (no wavefunctions at the moment)

      call petsc_get_sqrt_mat(p_gammar, PETSC_NULL_MAT, p_sqrt_gammar)
      call MatMatMult(p_g13, p_sqrt_gammar, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, p_tmp, ierr)
      call MatDuplicate(p_gammar, MAT_SHARE_NONZERO_PATTERN, p_tmp2, ierr)

      call petsc_aXpY(p_tmp2, p_tmp, p_one, petsc_false, .false.)
!~       call MatCopy(p_tmp,p_tmp2,DIFFERENT_NONZERO_PATTERN,ierr)
      call MatTranspose(p_tmp, MAT_INPLACE_MATRIX, p_tmp, ierr); CHKERRQ(ierr)
      call MatConjugate(p_tmp, ierr)

      call MatMatMult(p_gammal, p_tmp2, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, p_tmp3, ierr)
      call MatMatMult(p_tmp, p_tmp3, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, p_tt, ierr)

      call MatGetTrace(p_tt, trans, ierr)
!~       write(pstr_out,fmt='(A,2e24.12)') "trans !!!!!! ",trans 
!~       call petsc_print_master()
      
      call diag_mat2(p_tt, PETSC_NULL_MAT, t_i(iik), PETSC_NULL_MAT, nmu_r, eps_which_in = EPS_LARGEST_REAL)

      call MatDestroy(p_tt, ierr)
      call MatDestroy(p_tmp2, ierr)
      call MatDestroy(p_tmp3, ierr)
      call MatDestroy(p_sqrt_gammar, ierr)
    end if
! ----


    call MatMatMult(p_g13, p_gammar, MAT_INITIAL_MATRIX, 1d0, p_g13gammar, ierr)

!~     call MatHermitianTranspose(p_g13,MAT_INPLACE_MATRIX,p_g13,ierr)
!~     call MatMatMult(p_g13,p_gammal,MAT_INITIAL_MATRIX,1d0,p_g31conjgammal,ierr)

    call MatMatMult(p_gammal, p_g13, MAT_INITIAL_MATRIX, 1d0, p_gammalg13, ierr)

    call petsc_matmatconj_restricted(p_g13gammar, p_gammalg13, p_trace)


    call MatGetTrace(p_trace, trans, ierr)
!~    write(0,*) "trace ",iik, trans
    call MatDestroy(p_g13, ierr)
    call MatDestroy(p_g13gammar, ierr)
    call MatDestroy(p_gammalg13, ierr)
!~     call MatDestroy(p_trace,ierr)

  end subroutine get_transmission

  subroutine get_transmission_v2(x, trans, t_i, dosc, dosl, dosr)
#include <slepc/finclude/slepceps.h>
    use petsc
    use petsc_mod
    use petsc_wrapper
    use globals
    use integrator_mod
    use slepc_mod
    use eigenchannel_mod

    real(dp), intent(in) :: x
    complex(dp), intent(out) :: trans
    Vec, optional :: t_i(:)
    real(dp), intent(out), optional :: dosc, dosl, dosr

    logical :: ldosl, ldosr, ldosc

    Mat :: f, p_gammal_block, p_gammar_block, p_x, b, p_g_block, p_rhs, &
      p_g13, p_g13gammar, p_g31conjgammal, p_gammalg13, p_tmp, p_sqrt_gammar, &
      p_tt, p_tmp2, p_tmp3, p_mone, p_g31gammal

    Mat, pointer :: p_gr_inv

    complex(dp) :: z

    integer :: ierr, i
    real(dp) :: vbb
    PetscScalar :: p_scal
    PetscScalar, pointer :: xx_v(:)
    PetscReal :: norm

    integer(8) :: counti, count_rate, countf

    ldosc = present(dosc)
    ldosl = present(dosl)
    ldosr = present(dosr)

    z = x

    p_gr_inv => p_invGr ! see eval_gr

    call init_gr(z, p_gr_inv)

    call petsc_get_a_with_b_c(p_g_block, p_gr_inv, p_gammar, mattype_dense)

    solver_mode = 4 ! for now assume solver mode 1
    if (solver_mode .eq. 2) then
      call petsc_get_a_with_b_c(p_rhs, p_gr_inv, p_gammar, mattype_dense)
      call petsc_one_mat(p_rhs, nmu_c - nmu_r + 1 - 1, nmu_c - 1) ! prepare RHS nmu_c,nmu_r ( .., 0 0 .. 1 0, 0 0 .. 0 1 )   C style 1=0
      call petsc_call_solver(p_gr_inv, p_rhs, p_g_block, matsolvertype_cc, solver_mode)
    else if ((solver_mode .eq. 1).or.(solver_mode .eq. 3).or.(solver_mode .eq. 4)) then
      call petsc_get_a_with_b_c(p_rhs, p_gr_inv, p_gammar, mattype_sparse)
      call petsc_one_mat(p_rhs, nmu_c - nmu_r + 1 - 1, nmu_c - 1,2) ! prepare RHS nmu_c,nmu_r ( .., 0 0 .. 1 0, 0 0 .. 0 1 )   C style 1=0
      call petsc_call_solver(p_gr_inv, p_rhs, p_g_block, matsolvertype_cc, solver_mode)
    end if

    call MatDestroy(p_rhs, ierr)

!~     call petsc_mat_info(p_gammar,    "p_gammar", ierr)
!~     call petsc_mat_info(p_gammal,    "p_gammal", ierr)
!~     call petsc_mat_info(p_g_block,   "p_gblock", ierr)

    call petsc_get_a_with_b_c(p_g13, p_gammal, p_gammar, mattype_dense)

    call petsc_split_matrix(p_g_block, p_tmp, 1, nmu_l, 1, nmu_r, MAT_INITIAL_MATRIX)
    
!~     call petsc_mat_info(p_tmp,   "p_tmp   ", ierr)
    call petsc_add_sub_B_to_A(p_tmp, p_g13, 0, 0, p_one, INSERT_VALUES, PETSC_FALSE)
!~     call MatView(p_g13, PETSC_VIEWER_STDOUT_SELF, ierr)
!~     call petsc_mat_info(p_g13,  "p_g13   ", ierr)
    call MatDestroy(p_tmp, ierr)
    call MatDestroy(p_g_block, ierr)
    
    call MatMatMult(p_g13, p_gammar, MAT_INITIAL_MATRIX, 1d0, p_g13gammar, ierr)
!~     call petsc_mat_info(p_g13gammar,  "p_g13gr ", ierr)
    call MatConjugate(p_g13, ierr)
    call MatTransposeMatMult(p_g13, p_gammal, MAT_INITIAL_MATRIX, 1d0, p_g31gammal, ierr)
!~     call petsc_mat_info(p_g31gammal,  "p_g31gr ", ierr)
    call MatMatMult(p_g13gammar, p_g31gammal, MAT_INITIAL_MATRIX, 1d0, p_tmp, ierr)
!~     call petsc_mat_info(p_tmp,  "p_trace ", ierr)
    
    
    
    
!~     call slepcfinalize(ierr)
!~     stop    

    call MatGetTrace(p_tmp, trans, ierr)
!~    write(0,*) "trace ",iik, trans
    call MatDestroy(p_g13, ierr)
    call MatDestroy(p_g13gammar, ierr)
    call MatDestroy(p_g31gammal, ierr)
    call MatDestroy(p_tmp, ierr)
!~     call MatDestroy(p_trace,ierr)

  end subroutine get_transmission_v2

  subroutine transmission()
#include <petsc/finclude/petsc.h>
    use petscmat
    use petsc_mod
    use globals

    integer(8) :: counti, count_rate, countf
    integer :: ie, ierr, i1, ik, nl1, nl2
    complex(dp) :: z, trans, zdoscc, transk, tisum
    real(dp) :: energy, de, energy_restart
    integer :: iunit_dosl, iunit_dosr, iunit_dosc, iunit_trans
    Vec, allocatable :: p_tik(:)
    PetscScalar, pointer :: xx_v(:)
    PetscScalar, allocatable :: ti(:)
    logical :: l_trans_restart

    call petsc_get_densemat(p_h00k_r(1), p_grrr, mattype_surf)     ! elemental
    call petsc_get_densemat(p_h00k_l(1), p_gllr, mattype_surf)     ! elemental
    call petsc_get_densemat(p_h00k_r(1), p_sigmarr, mattype_dense) ! dense
    call petsc_get_densemat(p_h00k_r(1), p_gammar, mattype_dense)  ! dense
    call petsc_get_densemat(p_h00k_l(1), p_sigmalr, mattype_dense) ! dense
    call petsc_get_densemat(p_h00k_l(1), p_gammal, mattype_dense)  ! dense

    call petsc_get_a_with_b_c(p_trace, p_gammal, p_gammal, mattype_cc) ! setup AIJ matrix dim(nmu_l,nmu_l) which is going to hold only the trace, no allocation is made here

    call MatShift(p_trace, p_one, ierr); call MatZeroEntries(p_trace, ierr) ! use MatShift to allocate the trace and set it to 0

    if (l_ionode) then
      inquire (file=trim(trans_file), exist=l_trans_restart)
      write (pstr_out, fmt='(A,l)') "transmission file exist ", l_trans_restart
      call petsc_print_master()
      if (l_trans_restart) then
        energy_restart = 0d0
        l_trans_restart = .false.
        open (newunit=iunit_dosr, file="dos_r.dat", action="write", status="old", position="append")
        open (newunit=iunit_dosl, file="dos_l.dat", action="write", status="old", position="append")
        open (newunit=iunit_dosc, file="dos_cc.dat", action="write", status="old", position="append")
        open (newunit=iunit_trans, file=trim(trans_file), action="readwrite", status="old", position="rewind")
        do
          read (unit=iunit_trans, fmt=*, iostat=ierr) energy
          if (ierr .ne. 0) exit
          l_trans_restart = .true.
          energy_restart = energy
        end do
        write (pstr_out, fmt='(A,e24.12)') "last energy found ", energy_restart
        call petsc_print_master()
      else
        open (newunit=iunit_dosr, file="dos_r.dat", action="write", status="replace", position="append")
        open (newunit=iunit_dosl, file="dos_l.dat", action="write", status="replace", position="append")
        open (newunit=iunit_dosc, file="dos_cc.dat", action="write", status="replace", position="append")
        open (newunit=iunit_trans, file=trim(trans_file), action="write", status="replace", position="append")
      end if
    end if

    call MPI_bcast(energy_restart, 1, MPI_DOUBLE_PRECISION, 0, PETSC_COMM_WORLD, ierr)
    call MPI_bcast(l_trans_restart, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr)

    if (lget_ti) then
      allocate (p_tik(nk_c), ti(nmu_r))
      do ik = 1, nk_c
        call MatCreateVecs(p_gammar, p_tik(ik), PETSC_NULL_VEC, ierr)
      end do
    end if

    de = (eend - estart)/real(n_energy_steps, 8)
    energy = estart

    if (l_trans_restart) then
      i_energy_start = nint((energy_restart - (estart))/de) + 1
      write (pstr_out, fmt='(A,i8,e24.12)') "restart transmission calculation ", i_energy_start, energy_restart
      call petsc_print_master()
    end if

    do ie = i_energy_start, i_energy_end
      energy = estart + de*real(ie, 8)
      z = energy
      trans = 0d0
      ti = 0d0
      zdosl = 0d0
      zdosr = 0d0
      zdoscc = 0d0
      call system_clock(counti, count_rate)

      do ik = 1, nk_c

        iik = ik
        call get_transmission(energy, transk, p_tik)

        if (lget_ti) then
          call VecSum(p_tik(ik), tisum, ierr)
          call VecGetOwnershipRange(p_tik(ik), nl1, nl2, ierr)
          call VecGetArrayReadF90(p_tik(ik), xx_v, ierr)
          ti(nl1 + 1:nl2) = ti(nl1 + 1:nl2) + xx_v(:)*real(wkp_r(ik), dp)
          pstr_out = ""
          if (l_ionode) then
            write (pstr_out, fmt='(i8,A,2e16.8)') ie, " ti(E,k_n) ", real(z), real(tisum)
          end if
          do i1 = 1, size(xx_v)
            if ((abs(real(xx_v(i1))/real(transk)) .lt. 0.01) .or. (len(trim(pstr_out)) + 16 .ge. 256)) exit
            write (pstr_out, fmt='(A,e16.8)') trim(pstr_out), real(xx_v(i1))
          end do
          call PetscSynchronizedPrintf(PETSC_COMM_WORLD, trim(pstr_out), ierr)
          call PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT, ierr)

          write (pstr_out, fmt=*) ""; call petsc_print_master(.true.)
          call VecRestoreArrayReadF90(p_tik(iik), xx_v, ierr)

        end if
        write (pstr_out, fmt='(i8,A,e16.8,e16.8,A,i6,2e16.8,f8.4)') ie, " T(E,k_n)  ", real(z), real(transk), " kxy_n, w_n ", ik, kp_l(1:2, ik), real(wkp_r(ik), dp); call petsc_print_master()
        trans = trans + transk*real(wkp_r(ik), dp)
        if (lget_ti) then
          if (inode .eq. 0) then
            call MPI_REDUCE(MPI_IN_PLACE, ti, size(ti), MPI_DOUBLE_COMPLEX, MPI_SUM, 0, PETSC_COMM_WORLD, ierr)
          else
            call MPI_REDUCE(ti, ti, size(ti), MPI_DOUBLE_COMPLEX, MPI_SUM, 0, PETSC_COMM_WORLD, ierr)
          end if
        end if
      end do
      call system_clock(countf)
      trans = trans/real(nkx_dez*nky_dez, dp)
      if (lget_ti) ti = ti/real(nkx_dez*nky_dez, dp)

      if (l_ionode) then
        if (abs(aimag(trans)) .ge. 1d-12) then
          write (pstr_out, fmt='(i8,A,2e16.8)') ie, "warning large imaginary part for transmission ", trans
          call petsc_print_master()
        end if
        if (lget_ti) then
          write (pstr_out, fmt='(i8,A,2e16.8)') ie, " ti(E)     ", real(z), real(sum(ti))
          do i1 = 1, size(ti)
            if ((abs(real(ti(i1))/real(trans)) .lt. 0.01) .or. (len(trim(pstr_out)) + 16 .ge. 256)) exit
            write (pstr_out, fmt='(A,e16.8)') trim(pstr_out), real(ti(i1))
          end do
        end if
        write (pstr_out, fmt='(i8,A,2e16.8,A,e16.8)') ie, " T(E)      ", real(z), real(trans), " timing: ", real(countf - counti, dp)/real(count_rate, dp); call petsc_print_master()
        pstr_out = ""; call petsc_print_master(.true.)
        if (lget_ti) then
          write (unit=iunit_trans, fmt='(22es45.24e5)') real(z, 8), real(trans), real(ti(1:min(20, size(ti))))
        else
          write (unit=iunit_trans, fmt='(22es45.24e5)') real(z, 8), real(trans)
        end if
!~         zdosl=zdosl/real(nkx_dez*nky_dez,dp)
!~         zdosr=zdosr/real(nkx_dez*nky_dez,dp)
!~         zdoscc=zdoscc/real(nkx_dez*nky_dez,dp)
        write (unit=iunit_dosr, fmt='(3es45.24e5)') energy, -aimag(zdosr(1))/pi, -aimag(zdosr(2))/pi !,zdosl!-aimag(zdosat_r)/pi,-aimag(zdosat_l)/pi
        write (unit=iunit_dosl, fmt='(3es45.24e5)') energy, -aimag(zdosl(1))/pi, -aimag(zdosl(2))/pi !,zdosl!-aimag(zdosat_r)/pi,-aimag(zdosat_l)/pi
        write (unit=iunit_dosc, fmt='(2es45.24e5)') energy, -aimag(zdoscc)/pi

        call flush (iunit_dosr)
        call flush (iunit_dosl)
        call flush (iunit_dosc)
        call flush (iunit_trans)
      end if

    end do

    if (inode .eq. 0) then
      close (iunit_dosr)
      close (iunit_dosl)
      close (iunit_dosc)
      close (iunit_trans)
    end if

    if (lget_ti) then
      do ik = 1, nk_c
        call VecDestroy(p_tik(ik), ierr)
      end do
    end if

    call MatDestroy(p_grrr, ierr)
    call MatDestroy(p_gllr, ierr)
    call MatDestroy(p_sigmarr, ierr)
    call MatDestroy(p_gammar, ierr)
    call MatDestroy(p_sigmalr, ierr)
    call MatDestroy(p_gammal, ierr)

    call MatDestroy(p_trace, ierr)

  end subroutine transmission

  subroutine transmission_k_on_demand()
#include <petsc/finclude/petsc.h>
    use petscmat
    use petsc_mod
    use k_on_demand_mod
    use globals
    use error_handler

    integer(8) :: counti, count_rate, countf
    integer :: ie, ierr, i1, ik, nl1, nl2, i_energy_start_new, ik_start, ik_done
    integer :: ne, ne_in, n_done,n_total
    complex(dp) :: z, trans, zdoscc, transk, tisum
    real(dp) :: energy, de, energy_restart, energy_in, e1, e2, timing, timing_sum, time_remain
    integer :: iunit_dosl, iunit_dosr, iunit_dosc, iunit_trans, iunit_k_trans, iunit_trans_in
    Vec, allocatable :: p_tik(:)
    PetscScalar, pointer :: xx_v(:)
    PetscScalar, allocatable :: ti(:)
    logical :: l_trans_restart, l_exist
    real(dp), allocatable :: trans_in1(:), trans_in2(:)
    character(256) :: instr

    inquire (file=trim(trans_file), exist=l_trans_restart)
    if (l_trans_restart) then
      write (pstr_out, fmt='(A)') "transmission file exists no need to recalculate"
      call petsc_print_master()
      return
    end if

    call MPI_BARRIER(PETSC_COMM_WORLD, ierr)

    if (l_ionode) then
      inquire (file="trans.k.tmp", exist=l_trans_restart)
      write (pstr_out, fmt='(A,l)') "transmission file exist ", l_trans_restart
      call petsc_print_master()
      if (l_trans_restart) then
        energy_restart = 0d0
        l_trans_restart = .false.
!~         open (newunit=iunit_dosr, file="dos_r.k.dat", action="write", status="old", position="append")
!~         open (newunit=iunit_dosl, file="dos_l.k.dat", action="write", status="old", position="append")
!~         open (newunit=iunit_dosc, file="dos_cc.k.dat", action="write", status="old", position="append")
        open (newunit=iunit_k_trans, file="trans.k.tmp", action="read", status="old", position="rewind")
        read (unit=iunit_k_trans, fmt=*, iostat=ierr) instr, ik_start
        if (ierr .ne. 0) ik_start = 1
        do
          read (unit=iunit_k_trans, fmt=*, iostat=ierr) energy
          if (ierr .ne. 0) exit
          l_trans_restart = .true.
          energy_restart = energy
        end do
        close (iunit_k_trans)
        write (pstr_out, fmt='(A,i8,e24.12)') "at k-point last energy found ", ik_start, energy_restart
        call petsc_print_master()
      else
!~         open (newunit=iunit_dosr, file="dos_r.k.dat", action="write", status="replace", position="append")
!~         open (newunit=iunit_dosl, file="dos_l.k.dat", action="write", status="replace", position="append")
!~         open (newunit=iunit_dosc, file="dos_cc.k.dat", action="write", status="replace", position="append")
        open (newunit=iunit_k_trans, file="trans.k.tmp", action="write", status="replace", position="append")
        ik_start = 1
      end if
      ie = 0
      if (lget_ti) ie = min(20, nmu_l)
      allocate (trans_in1(1 + ie), trans_in2(1 + ie))
      trans_in1 = 0d0
      inquire (file="trans.tmp", exist=l_exist)
      if ((.not. l_exist) .and. (ik_start .eq. 1)) then
!~ set up a dummy trans.tmp file for step 0.
        open (newunit=iunit_trans_in, file="trans.tmp", action="write", status="new")
        de = (eend - estart)/real(n_energy_steps, 8)
        do ie = i_energy_start, i_energy_end
          energy = estart + de*real(ie, 8)
          write (unit=iunit_trans_in, fmt='(22es45.24e5)') energy, trans_in1
        end do
        close (iunit_trans_in)
      else if ((.not. l_exist) .and. (ik_start .ne. 1)) then
        write (pstr_out, fmt='(A)') "missing transmission from previous k-points delete trans.k.tmp and restart"
        call petsc_print_master()
        errormsg = pstr_out
        call error()
      end if

    end if

    call MPI_bcast(ik_start, 1, MPI_DOUBLE_PRECISION, 0, PETSC_COMM_WORLD, ierr)
    call MPI_bcast(energy_restart, 1, MPI_DOUBLE_PRECISION, 0, PETSC_COMM_WORLD, ierr)
    call MPI_bcast(l_trans_restart, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr)

    de = (eend - estart)/real(n_energy_steps, 8)
    energy = estart
    i_energy_start_new = i_energy_start

    if (l_trans_restart) then
      i_energy_start_new = nint((energy_restart - (estart))/de) + 1
      i_energy_start_new = min(i_energy_start_new, i_energy_end)
      write (pstr_out, fmt='(A,2i8,2e24.12)') "restart transmission calculation ", ik_start, i_energy_start, energy_restart, (energy_restart - (estart))/de
      call petsc_print_master()
    end if

    call petsc_get_densemat(p_h00_ik_r(1), p_grrr, mattype_surf)     ! elemental
    call petsc_get_densemat(p_h00_ik_l(1), p_gllr, mattype_surf)     ! elemental
    call petsc_get_densemat(p_h00_ik_r(1), p_sigmarr, mattype_dense) ! dense
    call petsc_get_densemat(p_h00_ik_r(1), p_gammar, mattype_dense)  ! dense
    call petsc_get_densemat(p_h00_ik_l(1), p_sigmalr, mattype_dense) ! dense
    call petsc_get_densemat(p_h00_ik_l(1), p_gammal, mattype_dense)  ! dense


    if (lget_ti) then
      allocate (p_tik(nk_c), ti(nmu_r))
      do ik = 1, nk_c
        call MatCreateVecs(p_gammar, p_tik(ik), PETSC_NULL_VEC, ierr)
      end do
    end if

    call petsc_get_a_with_b_c(p_trace, p_gammal, p_gammal, mattype_cc) ! setup AIJ matrix dim(nmu_l,nmu_l) which is going to hold only the trace, no allocation is made here

    call MatShift(p_trace, p_one, ierr); call MatZeroEntries(p_trace, ierr) ! use MatShift to allocate the trace and set it to 0

    timing_sum=0d0
    n_done=0
    n_total=(nk_c-ik_start)*(i_energy_end-i_energy_start+1)+(i_energy_end-i_energy_start_new+1)

    write (pstr_out, fmt='(A,i16,4i8)') "n_total ",n_total,nk_c,ik_start,i_energy_end,i_energy_start_new
    call petsc_print_master()


    do ik = ik_start, nk_c
      iik = ik
      if (l_k_on_demand) call k_on_demand(ik, 1)
      if (l_ionode) then
        if (l_trans_restart) then
          open (newunit=iunit_k_trans, file="trans.k.tmp", action="readwrite", status="old", position="append")
        else
          open (newunit=iunit_k_trans, file="trans.k.tmp", action="readwrite", status="replace")
          write (unit=iunit_k_trans, fmt='(A,i8)') "#", ik
        end if
      end if

      do ie = i_energy_start_new, i_energy_end
        n_done=n_done+1
        energy = estart + de*real(ie, 8)
        z = energy
        trans = 0d0
        ti = 0d0
        call system_clock(counti, count_rate)

        call get_transmission_v2(energy, transk, p_tik)

        if (lget_ti) then
          call VecSum(p_tik(ik), tisum, ierr)
          call VecGetOwnershipRange(p_tik(ik), nl1, nl2, ierr)
          call VecGetArrayReadF90(p_tik(ik), xx_v, ierr)
          write(0,*) inode, xx_v
          ti(nl1 + 1:nl2) = ti(nl1 + 1:nl2) + xx_v(:)
!~           pstr_out = " "
!~           if (l_ionode) then
!~             write (pstr_out, fmt='(i8,A,2e16.8)') ie, " ti(E,k_n) ", real(z), real(tisum)
!~             call petsc_print_master(.false.)
!~           end if
          
!~           pstr_out = " fuck "
!~           do i1 = 1, size(xx_v)
!~             if ((abs(real(xx_v(i1))/real(transk)) .lt. real(tisum)*0.01d0) .or. (len(trim(pstr_out)) + 16 .ge. 256)) exit
!~             write (pstr_out, fmt='(A,e16.8)') trim(pstr_out), real(xx_v(i1))
!~           end do
!~           pstr_out = adjustl(pstr_out)
!~           if (trim(pstr_out).ne."") call PetscSynchronizedPrintf(PETSC_COMM_WORLD, pstr_out, ierr)
!~           call PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT, ierr)

!~           call VecRestoreArrayReadF90(p_tik(iik), xx_v, ierr)
!~           pstr_out=""
!~           call petsc_print_master()

        end if
        call system_clock(countf)
        timing = real(countf - counti, dp)/real(count_rate, dp)
        timing_sum = timing_sum + timing
        time_remain = timing_sum/real(n_done, dp)*real(n_total-n_done,dp)/60d0/60d0
        write (pstr_out, fmt='(i8,A,e16.8,e16.8,A,i6,2e16.8,f8.4,2e14.6,2i8)') &
          ie, " T(E,k_n)  ", real(z), real(transk), " kxy_n, w_n ", ik, kp_l(1:2, ik), &
          real(wkp_r(ik), dp), timing, time_remain, n_total, n_done
        call petsc_print_master()
        
        if (lget_ti) then
          if (inode .eq. 0) then
            call MPI_REDUCE(MPI_IN_PLACE, ti, size(ti), MPI_DOUBLE_COMPLEX, MPI_SUM, 0, PETSC_COMM_WORLD, ierr)
          else
            call MPI_REDUCE(ti, ti, size(ti), MPI_DOUBLE_COMPLEX, MPI_SUM, 0, PETSC_COMM_WORLD, ierr)
          end if
        end if

        if (l_ionode) then
          if (abs(aimag(trans)) .ge. 1d-12) then
            write (pstr_out, fmt='(i8,A,2e16.8)') ie, "warning large imaginary part for transmission ", trans
            call petsc_print_master()
          end if
          if (lget_ti) then
! this is stupid should be change          
            write (unit=iunit_k_trans, fmt='(22es45.24e5)') real(z, 8), real(transk), &
              real(ti(1:min(20, size(ti))))
              
            write (pstr_out, fmt='(i8,A,2e16.8)') ie, " ti(E,k_n) ", real(z), real(tisum)
            do i1 = 1, nmu_c
              if ((abs(real(ti(i1))/real(transk)) .lt. 0.001d0) .or. (len(trim(pstr_out)) + 16 .ge. 256)) exit
              write (pstr_out, fmt='(A, e16.8)') trim(pstr_out), real(ti(i1))
            end do
            call petsc_print_master()
            write(0, *) real(ti)
          else
            write (unit=iunit_k_trans, fmt='(22es45.24e5)') real(z, 8), real(transk)
          end if
          call flush (iunit_k_trans)
        end if

      end do ! energy_loop

      if (l_trans_restart) then
        l_trans_restart = .false.
        i_energy_start_new = i_energy_start
      end if

      if (l_ionode) then
!~ accumulate date for current k-point with previous ones
        call rename("trans.tmp", "trans_old.tmp")
        open (newunit=iunit_trans, file="trans.tmp", action="write", status="replace")
        open (newunit=iunit_trans_in, file="trans_old.tmp", action="read", status="old")
        rewind (iunit_k_trans)
        read (iunit_k_trans, *) instr
        do ie = i_energy_start, i_energy_end
          read (iunit_trans_in, *) e1, trans_in1
          read (iunit_k_trans, *) e2, trans_in2
          trans_in2 = trans_in2*real(wkp_r(ik), dp)/real(nkx_dez*nky_dez, dp)
          if (abs(e2 - e1) .gt. eps_real) then
            write (pstr_out, fmt='(A,2e24.12)') "energy mismatch in trans_old.tmp and trans.k.tmp", e1, e2
            call petsc_print_master()
            errormsg = pstr_out
            call error()
          end if
          trans_in1 = trans_in1 + trans_in2
          write (unit=iunit_trans, fmt='(22es45.24e5)') e1, trans_in1
        end do
        close (iunit_trans)
        close (iunit_k_trans)
        close (iunit_trans_in)
      end if

      call MPI_BARRIER(PETSC_COMM_WORLD, ierr)

    end do ! k_loop

    if (l_ionode) then
      call rename("trans.tmp", trim(trans_file))
      open (newunit=iunit_trans, file="trans.k.tmp", action="write", status="replace")
      close (iunit_trans, status="delete")
      open (newunit=iunit_trans, file="trans_out.tmp", action="write", status="replace")
      close (iunit_trans, status="delete")
      open (newunit=iunit_trans, file="trans_old.tmp", action="write", status="replace")
      close (iunit_trans, status="delete")
    end if

    if (lget_ti) then
      do ik = 1, nk_c
        call VecDestroy(p_tik(ik), ierr)
      end do
    end if

    call MatDestroy(p_grrr, ierr)
    call MatDestroy(p_gllr, ierr)
    call MatDestroy(p_sigmarr, ierr)
    call MatDestroy(p_gammar, ierr)
    call MatDestroy(p_sigmalr, ierr)
    call MatDestroy(p_gammal, ierr)

    call MatDestroy(p_trace, ierr)

  end subroutine transmission_k_on_demand

end module transmission_mod

