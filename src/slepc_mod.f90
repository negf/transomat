module slepc_mod

  implicit none

contains
  subroutine diag_mat(a, b, ew, ev, nev, eps_which_in)
#include <slepc/finclude/slepceps.h>
    use slepceps
    use petsc_mod
    use globals, only: inode, mattype_dense
    use kinds
    use error_handler

    implicit none

    Mat :: a, b, ev, aa, bb, cc, dd
    Vec :: ew
    EPSWhich, optional :: eps_which_in
    PetscInt :: nev

    PetscInt :: nconv, ncv, mpd, ni(1), mi(1), i, nl1, nl2, ml1, ml2, j, maxits
    integer :: ierr

    PetscInt, allocatable :: idxr(:), idxc(:), idxvec(:)
    EPSProblemType :: epsproblem
    EPS :: eps
    EPSType :: tname
    PetscScalar kr(1), ki(1)
    PetscReal :: tol, error_ev, norm
    Vec :: xr, xi, xd, xdd
    PetscScalar, pointer :: xx_v(:)
    EPSWhich :: eps_which
    character(256) :: sdummy
    real(dp) :: info(MAT_INFO_SIZE)

    eps_which = EPS_LARGEST_MAGNITUDE
    if (present(eps_which_in)) eps_which = eps_which_in

    call VecSetOption(ew, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE, ierr)

    if (ev .ne. PETSC_NULL_MAT) then
      call MatCreateVecs(ev, xr, PETSC_NULL_VEC, ierr)
      call MatCreateVecs(ev, xi, PETSC_NULL_VEC, ierr)
    end if

    if (b .eq. PETSC_NULL_MAT) then
      epsproblem = EPS_HEP
    else
      epsproblem = EPS_GHEP
    end if

    !     ** Create eigensolver context
    call EPSCreate(PETSC_COMM_WORLD, eps, ierr)

!~     call EPSSetConvergenceTest(eps, EPS_CONV_ABS, ierr)

    call EPSSetDimensions(eps, nev, PETSC_DEFAULT_INTEGER,&
   &  PETSC_DEFAULT_INTEGER, ierr)

!     ** Set operators. In this case, it is a standard eigenvalue problem
    call EPSSetOperators(eps, a, b, ierr)
    call EPSSetProblemType(eps, epsproblem, ierr)
    call EPSSetWhichEigenpairs(eps, eps_which, ierr)

!     ** Set solver parameters at runtime
    call EPSSetFromOptions(eps, ierr)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Solve the eigensystem
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     ** Optional: Get some information from the solver and display it
    call EPSGetType(eps, tname, ierr)

    call EPSSolve(eps, ierr)
    call EPSGetConverged(eps, nconv, ierr)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Display solution and clean up
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    allocate (idxc(1), idxvec(1), stat=ierr)
    if (ierr .ne. 0) then
      write (errormsg, *) "allocation error ", ierr
      call error()
    end if

    do i = 0, nev - 1
      if (ev .ne. PETSC_NULL_MAT) then
        call EPSGetEigenpair(eps, i, kr, PETSC_NULL_SCALAR, xr, PETSC_NULL_VEC, ierr)
      else
        call EPSGetEigenpair(eps, i, kr, PETSC_NULL_SCALAR, PETSC_NULL_VEC, PETSC_NULL_VEC, ierr)
      end if
      call EPSComputeError(eps, i, EPS_ERROR_RELATIVE, error_ev, ierr)

      call VecGetOwnershipRange(ew, ml1, ml2, ierr)
!~         if ((i.ge.ml1).and.(i.le.ml2-1)) then
!~           write(6,*) "ev ",i,"=",kr
!~         end if

      if (ev .ne. PETSC_NULL_MAT) then
        call VecGetOwnershipRange(xr, nl1, nl2, ierr)
        call VecGetArrayReadF90(xr, xx_v, ierr)
        allocate (idxr(0:nl2 - nl1 - 1), stat=ierr)
        if (ierr .ne. 0) then
          write (errormsg, *) "allocation error ", ierr
          call error()
        end if
        do j = 0, nl2 - nl1 - 1
          idxr(j) = j + nl1
        end do
        idxc(1) = i

        call MatSetValuesBlocked(ev, nl2 - nl1, idxr, 1, idxc, xx_v, INSERT_VALUES, ierr)

        deallocate (idxr)

      end if

      idxvec(1) = i*sign(1, ml2 - i)*sign(1, i - ml1 - 1)

      call VecSetValues(ew, 1, idxvec, kr, INSERT_VALUES, ierr)

    end do
    call VecAssemblyBegin(ew, ierr)
    call VecAssemblyEnd(ew, ierr)
    if (ev .ne. PETSC_NULL_MAT) then
      call MatAssemblyBegin(ev, MAT_FINAL_ASSEMBLY, ierr)
      call MatAssemblyEnd(ev, MAT_FINAL_ASSEMBLY, ierr)
      call VecRestoreArrayReadF90(xr, xx_v, ierr)
    end if

    if (ev .ne. PETSC_NULL_MAT) then
      call VecDestroy(xr, ierr)
      call VecDestroy(xi, ierr)
    end if

    call EPSDestroy(eps, ierr)

  end subroutine diag_mat

  subroutine diag_mat2(a_in, b_in, ew, ev, nev, eps_calc, eps_which_in, nev_calc)
#include <slepc/finclude/slepceps.h>
    use slepceps
    use petsc_mod
    use globals, only: inode, mattype_dense
    use kinds
    use error_handler

    implicit none

    Mat :: a_in, b_in, ev, aa, bb, cc, dd
    Mat :: a, b
    Vec :: ew
    EPSWhich, optional :: eps_which_in
    integer, optional :: nev_calc
    PetscInt :: nev
    real(dp), optional :: eps_calc

    PetscInt :: nconv, ncv, mpd, ni(1), mi(1), i, nl1, nl2, ml1, ml2, j, maxits, nsim, iter
    integer :: ierr, nrow, ncol, ii

    PetscInt, allocatable :: idxr(:), idxc(:), idxvec(:)
    EPSProblemType :: epsproblem
    EPS :: eps
    EPSType :: tname
    PetscScalar kr(1), ki(1)
    PetscReal :: tol, error_ev, norm
    Vec :: xr, xi, xd, xdd
    Vec, allocatable :: def_space(:)
    PetscScalar, pointer :: xx_v(:), ev_pointer(:, :), ev_pointer2(:, :)
    PetscScalar, target :: t_dummy(1, 1)
    logical :: lfinished_evs, l_eps_calc
    EPSWhich :: eps_which
    EPSConvergedReason :: reason
    character(256) :: sdummy
    real(dp) :: info(MAT_INFO_SIZE)
    MatType :: a_mattype, b_mattype 


    eps_which = EPS_LARGEST_MAGNITUDE
    if (present(eps_which_in)) eps_which = eps_which_in

    l_eps_calc = .false.
    if (present(eps_calc)) l_eps_calc = .true.

    allocate (def_space(nev), idxc(1), idxvec(1), stat=ierr)
    if (ierr .ne. 0) then
      write (errormsg, *) "allocation error ", ierr
      call error()
    end if

!~ SLEPc crashes for MATDENSE. So if MATDENSE convert to MATSCALAPACK here and back to MATDENSE at the end.
!~ MATSCALAPACK changes the distribution of local rows per process which causes problems with EW and EV not 
!~ not matching to A. So convert MATDENSE to MATAIJ, while stupid there will probably be anyway some
!~ conversion within SLEPc happening.
    call MatGetType(a_in, a_mattype, ierr)
    if (b_in .ne. PETSC_NULL_MAT) then
      call MatGetType(b_in, b_mattype, ierr)
    else
      b_mattype = MATMPIAIJ
    end if
    if (a_mattype .eq. MATMPIDENSE) then
      call MatConvert(a_in, MATMPIAIJ, MAT_INITIAL_MATRIX, a, ierr)
    else
      a = a_in
    end if
    if (b_mattype .eq. MATMPIDENSE) then
      call MatConvert(b_in, MATMPIAIJ, MAT_INITIAL_MATRIX, b, ierr)
    else
      b = b_in
    end if

    call VecSetOption(ew, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE, ierr)

    if (b .eq. PETSC_NULL_MAT) then
      epsproblem = EPS_HEP
    else
      epsproblem = EPS_GHEP !EPS_GHEP
    end if

    !     ** Create eigensolver context
    call EPSCreate(PETSC_COMM_WORLD, eps, ierr)

!~     nsim = min(100,nev)

!~     call EPSSetConvergenceTest(eps, EPS_CONV_ABS, ierr)

    call EPSSetDimensions(eps, nev, PETSC_DEFAULT_INTEGER , PETSC_DEFAULT_INTEGER, ierr)
!~     call EPSSetDimensions(eps, nsim, ncv , mpd, ierr)

!     ** Set operators. In this case, it is a standard eigenvalue problem
    call EPSSetOperators(eps, a, b, ierr)
    call EPSSetProblemType(eps, epsproblem, ierr)
    call EPSSetWhichEigenpairs(eps, eps_which, ierr)

!     ** Set solver parameters at runtime
    call EPSSetFromOptions(eps, ierr)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Solve the eigensystem
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     ** Optional: Get some information from the solver and display it
    call EPSGetType(eps, tname, ierr)

    call VecGetOwnershipRange(ew, ml1, ml2, ierr)

! complicated fortran vs C index, probably chaning i=1 to i=0 would be better, but well
    i = 1
    if (ev .ne. PETSC_NULL_MAT) then
      call MatDenseGetArrayF90(ev, ev_pointer, ierr)
      if (size(ev_pointer) .eq. 0) then
        t_dummy = 0d0
        ev_pointer => t_dummy
      end if
    end if

    lfinished_evs = .false.

    call MPI_BARRIER(PETSC_COMM_WORLD, ierr)
    iter = 0
    call EPSGetDimensions(eps, nsim, ncv, mpd, ierr)
!~     call petsc_print_master()
!~     write(pstr_out,fmt='(A,5i8)') "Dims ", nev, nsim, ncv, mpd, nev
!~     call petsc_print_master()
    do
!~       write(pstr_out,fmt='(A)') "solve "
!~       call petsc_print_master()
      iter = iter + 1
!~       write(pstr_out,fmt='(i8)') iter
!~       call petsc_print_master()
      if (iter.ge.10000) exit
!~       call EPSGetDimensions(eps, nsim, ncv, mpd, ierr)
      call EPSSolve(eps, ierr)
      if (ierr .ne. 0) then
        write (errormsg, fmt='(A,i8)') "EPSSoler error diag_mat2 ", ierr
        call error()
      end if
      call EPSGetConverged(eps, nconv, ierr)
      call EPSGetConvergedReason(eps, reason, ierr)
      call EPSGetIterationNumber (eps, ii, ierr)
      if (ierr .ne. 0) then
        write (errormsg, fmt='(A,i8)') "EPSGetConverged error diag_mat2 ", ierr
        call error()
      end if

!~       write(pstr_out,fmt='(A,8i8)') "i, nconv, nsim ", i, nconv, nsim, ncv, mpd, ii, reason
!~       call petsc_print_master()

      nconv = min(nconv, nev)
!~       ncv = nconv * 2 + 1
!~       call EPSSetDimensions(eps, nconv, ncv , mpd, ierr)
      do j = i, i + nconv - 1

        if (ev .ne. PETSC_NULL_MAT) then
          call VecCreateMPIWithArray(PETSC_COMM_WORLD, 1, ml2 - ml1, PETSC_DECIDE, ev_pointer(:, j), def_space(j), ierr)
        else
          call MatCreateVecs(a, def_space(j), PETSC_NULL_VEC, ierr)
        end if
!~         write(pstr_out,fmt='(A,2i8)') "j- i",j,j - i
!~         call petsc_print_master()
        call EPSGetEigenpair(eps, j - i, kr, PETSC_NULL_SCALAR, def_space(j), PETSC_NULL_VEC, ierr)

        if (ierr .ne. 0) then
          write (errormsg, fmt='(A,i8)') "EPSGetEigenpair error diag_mat2 ", ierr
          call error()
        end if

        idxvec(1) = (j - 1)*sign(1, ml2 - (j - 1))*sign(1, (j - 1) - ml1 - 1)

        call VecSetValues(ew, 1, idxvec, kr, INSERT_VALUES, ierr)
!~         write(pstr_out,fmt='(A,2i8,2e24.12)') "ew ",j,kr
!~         call petsc_print_master()
        call VecAssemblyBegin(ew, ierr)
        call VecAssemblyEnd(ew, ierr)


        if (l_eps_calc) then
          if (abs(kr(1)) .le. eps_calc) then
!~             write(pstr_out,fmt='(A,2i8,2e24.12)') "ew i",i,j,kr ; call petsc_print_master()
            nconv = j - i + 1
            lfinished_evs = .true.
            exit
          end if
        end if
        if (j .eq. nev) then
          nconv = j - i + 1
          lfinished_evs = .true.
          exit
        end if
      end do

      i = i + nconv
!~       write(pstr_out,fmt='(A,2i8)') "Def space j,i",j, i
!~       call petsc_print_master()
      call EPSSetDeflationSpace(eps, i - 1, def_space(1:i-1), ierr);

      if ((i .ge. nev) .or. (lfinished_evs)) exit
    end do


    if (ev .ne. PETSC_NULL_MAT) call MatDenseRestoreArrayF90(ev, ev_pointer, ierr)

    call EPSDestroy(eps, ierr)
    do j = 1, i - 1
      call VecDestroy(def_space(j), ierr)
    end do
!
    if (a_mattype .eq. MATMPIDENSE) then
      call MatDestroy(a, ierr)
    end if
    if (b_mattype .eq. MATMPIDENSE) then
      call MatDestroy(b, ierr)
    end if

    if (present(nev_calc)) nev_calc = i

  end subroutine diag_mat2

end module slepc_mod
