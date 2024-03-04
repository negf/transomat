module petsc_wrapper

  implicit none

contains

  subroutine petsc_matassemble(a)
#include <petsc/finclude/petscmat.h>
    use petscmat
    implicit none

    Mat :: a
    integer :: ierr

    call MatAssemblyBegin(a, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(a, MAT_FINAL_ASSEMBLY, ierr)

  end subroutine petsc_matassemble

!~   subroutine petsc_matmatconj_restricted(A,B,C)
!~ C=A*B' that is we multiply row_A*conjg(row_B)
!~ C can probably take the petsc options offproc values disable
  subroutine petsc_matmatconj_restricted(A, B, C)
#include <petsc/finclude/petscmat.h>
    use petscmat
!~    use mpi
    use petsc_mod
    use globals, only: inode, nprocs
    use kinds
    use error_handler
    implicit none

    Mat :: A, B, C, D

    Mat :: Arow(0:nprocs - 1), Brow(0:nprocs - 1)
    PetscScalar, pointer :: p_array_B(:, :), p_array_A(:, :)
    PetscScalar, allocatable :: x_row_B(:)
    PetscScalar :: x(1)
    double complex :: zz(10)
    integer :: jnode, ierr, ncols, nrows, irow, icol, ir(1), ic(1), nrow_A, ncol_A, nrow_B, ncol_B, nrow_C, ncol_C, jcol
    integer :: nrow1_local_A, nrow2_local_A, nrow1_local_B, nrow2_local_B, nrow1_local_C, nrow2_local_C
    integer :: nzcols_C, nrow_max, iroot
    integer, allocatable ::C_cols(:)
    PetscErrorCode :: p_ierr

    call MatGetSize(A, nrow_A, ncol_A, ierr)
    call MatGetSize(B, nrow_B, ncol_B, ierr)
    call MatGetSize(C, nrow_C, ncol_C, ierr)

    call MatGetOwnershipRange(B, nrow1_local_B, nrow2_local_B, ierr)
    call MatGetOwnershipRange(A, nrow1_local_A, nrow2_local_A, ierr)
    call MatGetOwnershipRange(C, nrow1_local_C, nrow2_local_C, ierr)

    call MatSetOption(C, MAT_IGNORE_OFF_PROC_ENTRIES, PETSC_FALSE, ierr)

    allocate (x_row_B(ncol_B), C_cols(ncol_C), stat=ierr)
    if (ierr .ne. 0) then
      write (errormsg, fmt='(A,i8)') "allocation error ", ierr
      call error()
    end if

    call MatDenseGetArrayF90(B, p_array_B, p_ierr)
    call MatDenseGetArrayF90(A, p_array_A, p_ierr)

    do irow = 0, nrow_B - 1
      iroot = -1

      if ((irow .ge. nrow1_local_B) .and. (irow .lt. nrow2_local_B)) then
        x_row_B(1:ncol_B) = p_array_B(irow - nrow1_local_B + 1, 1:ncol_B)
        iroot = inode
      end if
      call MPI_AllReduce(MPI_IN_PLACE, iroot, 1, MPI_INTEGER, MPI_MAX, PETSC_COMM_WORLD, ierr)  ! is this actually a clever way to distribute iroot? looks stupid.
      call MPI_Bcast(x_row_B, ncol_B, MPI_DOUBLE_COMPLEX, iroot, PETSC_COMM_WORLD, ierr)

      iroot = -1
      nzcols_C = -1
      if ((irow .ge. nrow1_local_C) .and. (irow .lt. nrow2_local_C)) then
        call MatGetRow(C, irow, nzcols_C, C_cols, PETSC_NULL_SCALAR, p_ierr)
        iroot = inode
      end if
      call MPI_AllReduce(MPI_IN_PLACE, iroot, 1, MPI_INTEGER, MPI_MAX, PETSC_COMM_WORLD, ierr)  ! is this actually a clever way to distribute iroot? looks stupid.
      call MPI_Bcast(C_cols, ncol_C, MPI_INTEGER, iroot, PETSC_COMM_WORLD, ierr)
      call MPI_Bcast(nzcols_C, 1, MPI_INTEGER, iroot, PETSC_COMM_WORLD, ierr)

      ir = irow

      do icol = 1, nzcols_C
        jcol = C_cols(icol)
        if ((jcol .ge. nrow1_local_A) .and. (jcol .lt. nrow2_local_A)) then
          ic = jcol
!~           x(1)=dotc(x_row_B,p_array_A(jcol-nrow1_local_A+1,1:ncol_B))
          x(1) = dot_product(x_row_B, p_array_A(jcol - nrow1_local_A + 1, 1:ncol_B))
          call MatSetValues(C, 1, ic, 1, ir, x, INSERT_VALUES, ierr)
        end if
      end do

      call MatAssemblyBegin(C, MAT_FINAL_ASSEMBLY, ierr)
      call MatAssemblyEnd(C, MAT_FINAL_ASSEMBLY, ierr)

      if ((irow .ge. nrow1_local_C) .and. (irow .lt. nrow2_local_C)) &
     &call MatRestoreRow(C, irow, nzcols_C, C_cols, PETSC_NULL_SCALAR, p_ierr)

    end do

    call MatDenseRestoreArrayF90(B, p_array_B, p_ierr)
    call MatDenseRestoreArrayF90(A, p_array_A, p_ierr)

    deallocate (x_row_B, C_cols, stat=ierr)
    if (ierr .ne. 0) then
      write (errormsg, fmt='(A,i8)') "deallocation error ", ierr
      call error()
    end if

    nullify (p_array_B)
    nullify (p_array_A)
    call MatSetOption(C, MAT_IGNORE_OFF_PROC_ENTRIES, PETSC_FALSE, ierr)
  end subroutine petsc_matmatconj_restricted

  subroutine petsc_lu_fac(A, F, matsolvtype)
#include <petsc/finclude/petscmat.h>
    use petscmat
    implicit none

    Mat :: A, F
    MatSolverType :: matsolvtype

    IS :: isrow, iscol
    integer :: ierr
    PetscErrorCode :: p_ierr
    MatFactorInfo :: matfacinfo(MAT_FACTORINFO_SIZE)

!~     call MatGetOrdering(A,MATORDERINGNATURAL,isrow,iscol,ierr)
    call MatGetFactor(A, matsolvtype, MAT_FACTOR_LU, f, ierr)
    call MatFactorInfoInitialize(matfacinfo, p_ierr)

    call MatLUFactorSymbolic(F, A, PETSC_NULL_IS, PETSC_NULL_IS, matfacinfo, p_ierr)
!~     call MatLUFactorSymbolic(F,A,isrow,iscol,matfacinfo,p_ierr)
    call MatLUFactorNumeric(F, A, matfacinfo, p_ierr)

!~     call ISDestroy(isrow,ierr)
!~     call ISDestroy(iscol,ierr)

  end subroutine petsc_lu_fac

  subroutine petsc_ksp_lu_fac(ksp, A, F, comm)
#include <petsc/finclude/petsc.h>
    use petsc
    implicit none

    TYPE(MPI_Comm) :: comm
    Mat :: A, F
    KSP :: ksp

    PC :: pc
    integer :: ierr



    call KSPCreate(comm, ksp, ierr)
    call KSPSetOperators(ksp, A, A, ierr)
    call KSPSetType(ksp, KSPPREONLY, ierr)
    call KSPGetPC(ksp, pc, ierr)

    call PCSetType(pc, PCLU, ierr)
    call PCSetFromOptions(pc, ierr)
    call PCSetUp(pc, ierr)

    call KSPSetFromOptions(ksp, ierr)
    call KSPSetUp(ksp, ierr)


    call PCFactorSetUpMatSolverType(pc, ierr)
    call PCFactorGetMatrix(pc, F, ierr)

  end subroutine petsc_ksp_lu_fac


!~  reduce the number of simultanious RHS to nsim to avoid problems with SuperLU_dist
!~  and in the future for possibly very large matrices.
!~  basically using the pointer to the array (part) of the full matrix (MatDenseGetArrayF90 )
!~  and use the pointer as storage for the reduced matrix MatDensePlaceArray).
!~  This avoids creating and allocating memory as well as copying between matrices, this
!~  should only introduce negligible overhead.
  subroutine petsc_solve_cols(A, B, X, nsim_in, matsolvertype)
#include <petsc/finclude/petsc.h>
    use petsc
    use petsc_mod
    use globals, only: inode, nprocs, l_ionode, mattype_dense,&
    &inode_sub, inode_group, psubcomm, ngroups, nodes_group, mattype_surf
    use timing
    use misc
    use kinds
    use error_handler
    implicit none

    Mat :: A, B, X
    integer :: nsim_in
    MatSolverType :: matsolvertype

    integer :: ierr, nrow, ncol, nruns, irun, nrow_loc, nrow_loc1, nrow_loc2, &
      ioff, i, icol(1), nrest, nsim
    integer, allocatable :: irows_loc(:)
    Mat :: F, BB, XX
    KSP :: ksp
    PC :: pc
    PetscScalar, pointer :: p_bb(:,:),p_xx(:,:)

    integer(8) :: counti, count_rate, countf


    call MatSetOption(X, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE, ierr)
!~ I think this option is save to use as the temporary XX, BB have the same
!~ row distribution as X, B
    call MatSetOption(X, MAT_IGNORE_OFF_PROC_ENTRIES, PETSC_FALSE, ierr)

    call KSPCreate(PETSC_COMM_WORLD, ksp, ierr)
    call KSPSetOperators(ksp, A, A, ierr)
    call KSPSetType(ksp, KSPPREONLY, ierr)
    call KSPGetPC(ksp, pc, ierr)

    call PCSetType(pc, PCLU, ierr)
    call PCSetFromOptions(pc, ierr)
    call PCSetUp(pc, ierr)

    call KSPSetFromOptions(ksp, ierr)
    call KSPSetUp(ksp, ierr)

    call PCFactorSetUpMatSolverType(pc, ierr)
    call PCFactorGetMatrix(pc, F, ierr)


    call MatGetSize(B, nrow, ncol, ierr)
    nsim = min(ncol, nsim_in)


    call MatGetOwnershipRange(B, nrow_loc1, nrow_loc2, ierr)
    nrow_loc = nrow_loc2 - nrow_loc1

    allocate(irows_loc(nrow_loc), stat = ierr)
    if (ierr .ne. 0) then
      write (errormsg, fmt='(A,i8)') "allocation error ", ierr
      call error()
    end if

!~ I wonder if there is an easier way than to use this with MatSetValues
    do i = 1,nrow_loc
      irows_loc(i) = nrow_loc1 + i - 1
    end do


    call MatCreateDense(PETSC_COMM_WORLD, nrow_loc, PETSC_DECIDE, PETSC_DECIDE, nsim, &
      PETSC_NULL_SCALAR, XX, ierr)
    call MatDenseGetArrayF90(XX, p_xx, ierr)

    call MatDenseGetArrayF90(B, p_bb, ierr)
    call MatCreateDense(PETSC_COMM_WORLD, nrow_loc, PETSC_DECIDE, PETSC_DECIDE, nsim, &
      p_bb(1:nrow_loc, 1:nsim), BB, ierr)

    nruns = ncol / nsim
    nrest = ncol - nruns * nsim


    call system_clock(counti, count_rate)
    do irun = 0, nruns - 1
!~       call system_clock(counti, count_rate)
      ioff = 1 + nsim * irun
      call MatDensePlaceArray(BB, p_bb(1, ioff), ierr)
      call MatMatSolve(F, BB, XX, ierr)
      if (ierr.ne.0) then
        write (errormsg, fmt='(A,i8)') "solver rest error ", ierr
        call error()
      end if

      call MatDenseResetArray(BB, ierr)

      do i = 1, nsim
        icol = ioff + i - 2  ! -1 as ioff starts at 1 and -1 for 1 index to 0 index
        call MatSetValues(X, nrow_loc, irows_loc, 1, icol , p_xx(1:nrow_loc, i), &
          INSERT_VALUES, ierr)
      end do
      call MatAssemblyBegin(X, MAT_FINAL_ASSEMBLY, ierr) ! I think the values are cached.
      call MatAssemblyEnd(X, MAT_FINAL_ASSEMBLY, ierr)

!~       call system_clock(countf)
!~       write (pstr_out, fmt='(A,i8,e24.12)') "test ",nruns,real(countf - counti, 8)/real(count_rate, 8)
!~       call petsc_print_master()
    end do


    call MatDenseRestoreArrayF90(XX, p_xx, ierr)
    call MatDenseRestoreArrayF90(B, p_bb, ierr)




    if (nrest .gt. 0) then

      ioff = 1 + nsim * nruns

      if ( (ioff + nrest - 1) .ne. (ncol) ) then
        write (errormsg, fmt='(A,3i12)') "ioff + nrest - 1 .ne. ncol  ", ioff, nrest, ncol
        call error()
      end if

      call KSPDestroy(ksp,ierr)

      call KSPCreate(PETSC_COMM_WORLD, ksp, ierr)
      call KSPSetOperators(ksp, A, A, ierr)
      call KSPSetType(ksp, KSPPREONLY, ierr)
      call KSPGetPC(ksp, pc, ierr)

      call PCSetType(pc, PCLU, ierr)
      call PCSetFromOptions(pc, ierr)
      call PCSetUp(pc, ierr)

      call KSPSetFromOptions(ksp, ierr)
      call KSPSetUp(ksp, ierr)

      call PCFactorSetUpMatSolverType(pc, ierr)
      call PCFactorGetMatrix(pc, F, ierr)

      call MatDestroy(BB, ierr)
      call MatDestroy(XX, ierr)


      call MatCreateDense(PETSC_COMM_WORLD, nrow_loc, PETSC_DECIDE, PETSC_DECIDE, nrest, &
        PETSC_NULL_SCALAR, XX, ierr)
      call MatDenseGetArrayF90(XX, p_xx, ierr)

      call MatDenseGetArrayF90(B, p_bb, ierr)
      call MatCreateDense(PETSC_COMM_WORLD, nrow_loc, PETSC_DECIDE, PETSC_DECIDE, nrest, &
        p_bb(1:nrow_loc, ioff : ioff + nrest - 1), BB, ierr)

      call MatMatSolve(F, BB, XX, ierr)
      if (ierr.ne.0)  then
        write (errormsg, fmt='(A,i8)') "solver rest error ", ierr
        call error()
      end if
      do i = 1, nrest
        icol = ioff + i - 2  ! -1 as ioff starts at 1 and -1 for 1 index to 0 index
        call MatSetValues(X, nrow_loc, irows_loc, 1, icol , p_xx(1:nrow_loc, i), &
        INSERT_VALUES, ierr)
      end do
      call MatAssemblyBegin(X, MAT_FINAL_ASSEMBLY, ierr)
      call MatAssemblyEnd(X, MAT_FINAL_ASSEMBLY, ierr)

      call MatDenseRestoreArrayF90(XX, p_xx, ierr)
      call MatDenseRestoreArrayF90(B, p_bb, ierr)

    end if




    call system_clock(countf)
    timing_matmat_solve(1) = timing_matmat_solve(1) + real(countf - counti, 8)/real(count_rate, 8)

    call MatDestroy(BB, ierr)
    call MatDestroy(XX, ierr)

    call KSPDestroy(ksp,ierr)

  end subroutine petsc_solve_cols

  subroutine petsc_solve_cols_sparse_input(A, B, X, nsim_in, matsolvertype)
#include <petsc/finclude/petsc.h>
    use petsc
    use petsc_mod
    use globals, only: inode, nprocs, l_ionode, mattype_dense,&
    &inode_sub, inode_group, psubcomm, ngroups, nodes_group, mattype_surf
    use timing
    use misc
    use kinds
    use error_handler
    implicit none

    Mat :: A, B, X
    integer :: nsim_in
    MatSolverType :: matsolvertype

    integer :: ierr, nrow, ncol, nruns, irun, nrow_loc, nrow_loc1, nrow_loc2, &
      ioff, i, icol(1), nrest, nsim
    integer, allocatable :: irows_loc(:), icols_sim(:), icols_off(:)
    Mat :: F, BB, XX
    KSP :: ksp
    PC :: pc
    PetscScalar, pointer :: p_bb(:,:),p_xx(:,:)

    integer(8) :: counti, count_rate, countf


    call MatSetOption(X, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE, ierr)
!~ I think this option is save to use as the temporary XX, BB have the same
!~ row distribution as X, B
    call MatSetOption(X, MAT_IGNORE_OFF_PROC_ENTRIES, PETSC_TRUE, ierr)

    call KSPCreate(PETSC_COMM_WORLD, ksp, ierr)
    call KSPSetOperators(ksp, A, A, ierr)
    call KSPSetType(ksp, KSPPREONLY, ierr)
    call KSPGetPC(ksp, pc, ierr)

    call PCSetType(pc, PCLU, ierr)
    call PCSetFromOptions(pc, ierr)
    call PCSetUp(pc, ierr)

    call KSPSetFromOptions(ksp, ierr)
    call KSPSetUp(ksp, ierr)


    call PCFactorSetUpMatSolverType(pc, ierr)
    call PCFactorGetMatrix(pc, F, ierr)


    call MatGetSize(B, nrow, ncol, ierr)
    nsim = min(ncol, nsim_in)


    call MatGetOwnershipRange(B, nrow_loc1, nrow_loc2, ierr)
    nrow_loc = nrow_loc2 - nrow_loc1

    allocate(irows_loc(nrow_loc), icols_sim(nsim), icols_off(nsim), stat = ierr)
    if (ierr .ne. 0) then
      write (errormsg, fmt='(A,i8)') "allocation error ", ierr
      call error()
    end if

!~ I wonder if there is an easier way than to use this with MatSetValues
    do i = 1,nrow_loc
      irows_loc(i) = nrow_loc1 + i - 1
    end do
    do i = 1, nsim
      icols_sim(i) = i - 1
    end do


    call MatCreateDense(PETSC_COMM_WORLD, nrow_loc, PETSC_DECIDE, PETSC_DECIDE, nsim, &
      PETSC_NULL_SCALAR, XX, ierr)
    call MatDenseGetArrayF90(XX, p_xx, ierr)

    call MatCreateDense(PETSC_COMM_WORLD, nrow_loc, PETSC_DECIDE, PETSC_DECIDE, nsim, &
      PETSC_NULL_SCALAR, BB, ierr)
    call MatDenseGetArrayF90(BB, p_bb, ierr)


    call MatAssemblyBegin(BB, MAT_FINAL_ASSEMBLY, ierr) ! I think the values are cached.
    call MatAssemblyEnd(BB, MAT_FINAL_ASSEMBLY, ierr)


    nruns = ncol / nsim
    nrest = ncol - nruns * nsim


    call system_clock(counti, count_rate)
    do irun = 0, nruns - 1
!~       call system_clock(counti, count_rate)
      ioff = nsim * irun

      do i = 1, nsim
        icol = ioff + i - 1
        call MatGetValues(B, nrow_loc, irows_loc, 1, icol, p_bb(1:nrow_loc,i), ierr)
      end do
      call MatMatSolve(F, BB, XX, ierr)
      if (ierr.ne.0) then
        write (errormsg, fmt='(A,i8)') "solver rest error ", ierr
        call error()
      end if

      do i = 1, nsim
        icol = ioff + i - 1
        call MatSetValues(X, nrow_loc, irows_loc, 1, icol , p_xx(1:nrow_loc,i), &
          INSERT_VALUES, ierr)
      end do

      call MatAssemblyBegin(X, MAT_FLUSH_ASSEMBLY, ierr) ! I think the values are cached.
      call MatAssemblyEnd(X, MAT_FLUSH_ASSEMBLY, ierr)

!~       call system_clock(countf)
!~       write (pstr_out, fmt='(A,2i8,e24.12)') "test ",nruns,irun,real(countf - counti, 8)/real(count_rate, 8)
!~       call petsc_print_master()
    end do

    call MatDenseRestoreArrayF90(XX, p_xx, ierr)
    call MatDenseRestoreArrayF90(BB, p_bb, ierr)




    if (nrest .gt. 0) then

      ioff = nsim * nruns

      if ( (ioff + nrest ) .ne. (ncol) ) then
        write (errormsg, fmt='(A,3i12)') "ioff + nrest - 1 .ne. ncol  ", ioff, nrest, ncol
        call error()
      end if

      call KSPDestroy(ksp,ierr)

      call KSPCreate(PETSC_COMM_WORLD, ksp, ierr)
      call KSPSetOperators(ksp, A, A, ierr)
      call KSPSetType(ksp, KSPPREONLY, ierr)
      call KSPGetPC(ksp, pc, ierr)

      call PCSetType(pc, PCLU, ierr)
      call PCSetFromOptions(pc, ierr)
      call PCSetUp(pc, ierr)

      call KSPSetFromOptions(ksp, ierr)
      call KSPSetUp(ksp, ierr)

      call PCFactorSetUpMatSolverType(pc, ierr)
      call PCFactorGetMatrix(pc, F, ierr)

      call MatDestroy(BB, ierr)
      call MatDestroy(XX, ierr)


      call MatCreateDense(PETSC_COMM_WORLD, nrow_loc, PETSC_DECIDE, PETSC_DECIDE, nrest, &
        PETSC_NULL_SCALAR, XX, ierr)
      call MatDenseGetArrayF90(XX, p_xx, ierr)

      call MatCreateDense(PETSC_COMM_WORLD, nrow_loc, PETSC_DECIDE, PETSC_DECIDE, nrest, &
        PETSC_NULL_SCALAR, BB, ierr)
      call MatDenseGetArrayF90(BB, p_bb, ierr)
      call MatAssemblyBegin(BB, MAT_FINAL_ASSEMBLY, ierr)
      call MatAssemblyEnd(BB, MAT_FINAL_ASSEMBLY, ierr)


      do i = 1, nrest
        icol = ioff + i - 1
        call MatGetValues(B, nrow_loc, irows_loc, 1, icol, p_bb(1:nrow_loc,i), ierr)
      end do



      call MatMatSolve(F, BB, XX, ierr)
      if (ierr.ne.0)  then
        write (errormsg, fmt='(A,i8)') "solver rest error ", ierr
        call error()
      end if
      do i = 1, nrest
        icol = ioff + i - 1  ! -1 for 1 index to 0 index
        call MatSetValues(X, nrow_loc, irows_loc, 1, icol , p_xx(1:nrow_loc, i), &
        INSERT_VALUES, ierr)
      end do
      call MatAssemblyBegin(X, MAT_FINAL_ASSEMBLY, ierr)
      call MatAssemblyEnd(X, MAT_FINAL_ASSEMBLY, ierr)

      call MatDenseRestoreArrayF90(XX, p_xx, ierr)
      call MatDenseRestoreArrayF90(BB, p_bb, ierr)

    end if




    call system_clock(countf)
    timing_matmat_solve(1) = timing_matmat_solve(1) + real(countf - counti, 8)/real(count_rate, 8)

    call MatDestroy(BB, ierr)
    call MatDestroy(XX, ierr)

    call KSPDestroy(ksp,ierr)

!~     call dump_nonzero_structure(X, "X.dat",PETSC_COMM_WORLD, 0)
!~     write(0,*) "------finished------"
!~     call MPI_BARRIER(PETSC_COMM_WORLD,ierr)

!~     call slepcfinalize(ierr)
!~     stop
!~     call MatSetOption(X, MAT_IGNORE_OFF_PROC_ENTRIES, PETSC_FALSE, ierr)

  end subroutine petsc_solve_cols_sparse_input

  subroutine petsc_solve_sparse_input_onsubcoms(A_in, B_in, X, nsim_in, matsolvertype)
#include <petsc/finclude/petsc.h>
    use petsc
    use MPI
    use petsc_mod
    use globals, only: inode, nprocs, l_ionode, mattype_dense, group_range, &
      inode_sub, inode_group, psubcomm, ngroups, nodes_group, mattype_surf, &
      mattype_sparse
    use timing
    use misc
    use kinds
    use error_handler
    implicit none

    Mat :: A_in, B_in, X
    integer :: nsim_in
    MatSolverType :: matsolvertype

    integer :: ierr, nrow, ncol, nruns, irun, nrow_loc, nrow_loc1, nrow_loc2, j, nz, nsim_min, &
      ioff, i, icol(1), nrest, nsim, ncols_subcomm(2), nrow_loc_subcom, ncols_on_subcomm, &
      ii(1), jj(2), nrow_local_range(2), nrow_subcomm_loc1, nrow_subcomm_loc2, nruns_subcom
    integer, allocatable :: irows_loc(:), icols_sim(:), icols_off(:), cols_offset(:), &
    cols_count(:), irows_loc_subcom(:), irows(:), jrows(:)
    Mat :: F, BB, XX, A, B, Bt
    KSP :: ksp
    PC :: pc
    PetscScalar, pointer :: p_xx(:,:)
    PetscScalar, allocatable :: p_batch(:), p_batch_subcom(:)
    PetscScalar :: pp
    PetscReal :: norm
    integer(8) :: counti, count_rate, countf


    call MatSetOption(X, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE, ierr)
!~     call MatSetOption(B, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE, ierr)
!~ we need this as we are flipping cols and rows
!~     call MatSetOption(X, MAT_IGNORE_OFF_PROC_ENTRIES, PETSC_FALSE, ierr)
!~     call MatSetOption(B, MAT_IGNORE_OFF_PROC_ENTRIES, PETSC_FALSE, ierr)

!~ Transpose the matrix to have the columns as rows on the correct subcommunicator
!~ this could probably done more efficient, but this is super easy. so for now
!~ let it be.

!~ this should rebalance the columns for the processors, in particular if B and X are not
!~ square matrices. Maybe allocates a few to many NZ elements but should be fne.
!~ one could set up the correct preallocator matrix.
!~     call MatGetSize(B_in, nrow, ncol, ierr)
!~     call MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,ncol,nrow,0,PETSC_NULL_INTEGER, &
!~       0,PETSC_NULL_INTEGER,B,ierr)
!~     call MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY, ierr) ! I think the values are cached.
!~     call MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY, ierr)
!~     call petsc_aXpY(B, B_in, p_one, PETSC_TRUE, .true.)

    call petsc_get_alloc_preG(B_in, Bt, 0, 0, .true.)
    call petsc_get_a_with_b_c(B, Bt, Bt, mattype_sparse)
    call MatPreallocatorPreallocate(Bt, PETSC_TRUE, B, ierr)
    call MatSetOption (B, MAT_NO_OFF_PROC_ENTRIES, PETSC_FALSE, ierr)
    call MatDestroy(Bt, ierr)
    call petsc_aXpY(B, B_in, p_one, PETSC_FALSE, .true.)


!~     call MatTranspose(B_in,  MAT_INITIAL_MATRIX,B,ierr)
!~     call petsc_mat_info(B,"B ",ierr)
    call MatGetSize(B, nrow, ncol, ierr)
!~     write(0,*) "size ",nrow,ncol

    call MatGetOwnershipRange(B, nrow_loc1, nrow_loc2, ierr)
    nrow_loc = nrow_loc2 - nrow_loc1
    nrow_local_range(1) = nrow_loc1
    nrow_local_range(2) = nrow_loc2 - 1
!~     write(0,fmt='(A,5i8)') "B ",inode, inode_group,nrow_loc,nrow_loc1,nrow_loc2

    allocate(irows_loc(nrow_loc), irows(ncol), jrows(ncol), stat = ierr)
    if (ierr .ne. 0) then
      write (errormsg, fmt='(A,i8)') "allocation error ", ierr
      call error()
    end if

!~ I wonder if there is an easier way than to use this with MatSetValues
    do i = 1,nrow_loc
      irows_loc(i) = nrow_loc1 + i - 1
    end do

    do i = 1, ncol
      irows(i) = i - 1
    end do


    allocate (cols_offset(nprocs), cols_count(nprocs), stat=ierr)
    if (ierr .ne. 0) then
      write (errormsg, fmt='(A,i8)') "allocation error ", ierr
      call error()
    end if
!~     call petsc_cols_to_subcomms(ncol, ncols_subcomm, cols_offset, cols_count)
    call flush(0)
    call MPI_BARRIER(PETSC_COMM_WORLD,ierr)
    ncols_subcomm(1) = nrow_local_range(1)
    ncols_subcomm(2) = nrow_local_range(2)
    call MPI_bcast(ncols_subcomm(1), 1, MPI_INTEGER, 0, psubcomm, ierr)
    call MPI_bcast(ncols_subcomm(2), 1, MPI_INTEGER, nodes_group(inode_group) - 1, psubcomm, ierr)
    ncols_on_subcomm = ncols_subcomm(2) - ncols_subcomm(1) + 1
!~     write(0,fmt='(6i8)') inode, inode_group, inode_sub, ncols_on_subcomm, ncols_subcomm

    call MatCreateRedundantMatrix(A_in, 0, psubcomm, MAT_INITIAL_MATRIX, A, ierr)


    call MatGetOwnershipRange(A, nrow_subcomm_loc1, nrow_subcomm_loc2, ierr)
    nrow_loc_subcom = nrow_subcomm_loc2 - nrow_subcomm_loc1


    allocate(p_batch(ncol),stat = ierr)
    if (ierr .ne. 0) then
      write (errormsg, fmt='(A,i8)') "allocation error ", ierr
      call error()
    end if

    allocate(irows_loc_subcom(nrow_loc_subcom), stat = ierr)
    if (ierr .ne. 0) then
      write (errormsg, fmt='(A,i8)') "allocation error ", ierr
      call error()
    end if

!~ I wonder if there is an easier way than to use this with MatSetValues
    do i = 1,nrow_loc_subcom
      irows_loc_subcom(i) = nrow_subcomm_loc1 + i - 1
    end do


    nsim = min(ncols_on_subcomm, nsim_in)
    call MPI_Allreduce(nsim, nsim_min, 1, MPI_INTEGER, MPI_MIN, PETSC_COMM_WORLD, ierr)
    if (nsim_min .le. 0) then
      write (errormsg, fmt='(A,3i8)') "min(nsim)=0 on group, ngroups to large ",ncols_on_subcomm, nsim, nsim_min
      call error()
    end if

    call MatCreateDense(psubcomm, nrow_loc_subcom, PETSC_DECIDE, PETSC_DECIDE, nsim, &
      PETSC_NULL_SCALAR, XX, ierr)

    call MatCreateDense(psubcomm, nrow_loc_subcom, PETSC_DECIDE, PETSC_DECIDE, nsim, &
      PETSC_NULL_SCALAR, BB, ierr)

    call MatAssemblyBegin(BB, MAT_FINAL_ASSEMBLY, ierr) ! I think the values are cached.
    call MatAssemblyEnd(BB, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyBegin(XX, MAT_FINAL_ASSEMBLY, ierr) ! I think the values are cached.
    call MatAssemblyEnd(XX, MAT_FINAL_ASSEMBLY, ierr)
!~ we need this as we are flipping cols and rows
    call MatSetOption(BB, MAT_IGNORE_OFF_PROC_ENTRIES, PETSC_FALSE, ierr)


    call KSPCreate(psubcomm, ksp, ierr)
    call KSPSetOperators(ksp, A, A, ierr)
    call KSPSetType(ksp, KSPPREONLY, ierr)
    call KSPGetPC(ksp, pc, ierr)

    call PCSetType(pc, PCLU, ierr)
    call PCSetFromOptions(pc, ierr)
    call PCSetUp(pc, ierr)

    call KSPSetFromOptions(ksp, ierr)
    call KSPSetUp(ksp, ierr)


    call PCFactorSetUpMatSolverType(pc, ierr)
    call PCFactorGetMatrix(pc, F, ierr)


!~     write(0,*) inode,inode_group,ncols_on_subcomm
    nruns_subcom = ncols_on_subcomm / nsim
    nrest = ncols_on_subcomm - nruns_subcom * nsim

    call MPI_AllReduce(nruns_subcom, nruns, 1, MPI_INTEGER, MPI_MAX, PETSC_COMM_WORLD, ierr)
!~     write(0,fmt='(6i8)') inode, inode_group, ncols_on_subcomm, nruns_subcom, nruns, nrest

    call system_clock(counti, count_rate)

    do irun = 0, nruns - 1
!~       call system_clock(counti, count_rate)
      if (irun .le. nruns_subcom - 1) then
        ioff = nsim * irun + ncols_subcomm(1)

        call MatZeroEntries(BB,ierr)
        do i = 1, nsim

          icol = ioff + i - 1
          if ((icol(1) .ge. nrow_loc1) .and. (icol(1) .lt. nrow_loc2)) then
            j = icol(1)
            call MatGetRow(B,j,nz,jrows,p_batch,ierr)
            icol = i - 1
            call MatSetValues(BB, nz, jrows, 1, icol, p_batch(1:nz), INSERT_VALUES, ierr)
            if (ierr .ne. 0) then
              write (errormsg, fmt='(A,i8)') "MatSetValues BB ", ierr
              call error()
            end if
            call MatRestoreRow(B,j,nz,jrows,p_batch,ierr)
          end if

        end do
        call MatAssemblyBegin(BB, MAT_FINAL_ASSEMBLY, ierr) ! I think the values are cached.
        call MatAssemblyEnd(BB, MAT_FINAL_ASSEMBLY, ierr)

        call MatMatSolve(F, BB, XX, ierr)
        if (ierr .ne. 0) then
          write (errormsg, fmt='(A,i8)') "solver rest error ", ierr
          call error()
        end if

        call MatDenseGetArrayF90(XX,p_xx,ierr)
        do i = 1, nsim
          icol = ioff + i - 1
          call MatSetValues(X, nrow_loc_subcom, irows_loc_subcom, 1, icol ,&
            p_xx(1:nrow_loc_subcom,i), INSERT_VALUES, ierr)
          if (ierr .ne. 0) then
            write (errormsg, fmt='(A,i8)') "MatSetValues X ", ierr
            call error()
          end if
        end do
        call MatDenseRestoreArrayF90(XX,p_xx,ierr)
      end if


      call MatAssemblyBegin(X, MAT_FLUSH_ASSEMBLY, ierr) ! I think the values are cached.
      call MatAssemblyEnd(X, MAT_FLUSH_ASSEMBLY, ierr)

!~       call system_clock(countf)
!~       write (pstr_out, fmt='(A,2i8,e24.12)') "test ",nruns,irun,real(countf - counti, 8)/real(count_rate, 8)
!~       call petsc_print_master()
    end do




    if (nrest .gt. 0) then


      if ( (nsim * nruns_subcom  + nrest) .ne. (ncols_on_subcomm) ) then
        write (errormsg, fmt='(A,3i12)') "ioff + nrest - 1 .ne. ncol  ", ioff, nrest, ncols_on_subcomm
        call error()
      end if

      ioff = nsim * nruns_subcom + ncols_subcomm(1)

      call KSPDestroy(ksp,ierr)

      call KSPCreate(psubcomm, ksp, ierr)
      call KSPSetOperators(ksp, A, A, ierr)
      call KSPSetType(ksp, KSPPREONLY, ierr)
      call KSPGetPC(ksp, pc, ierr)

      call PCSetType(pc, PCLU, ierr)
      call PCSetFromOptions(pc, ierr)
      call PCSetUp(pc, ierr)

      call KSPSetFromOptions(ksp, ierr)
      call KSPSetUp(ksp, ierr)

      call PCFactorSetUpMatSolverType(pc, ierr)
      call PCFactorGetMatrix(pc, F, ierr)

      call MatDestroy(BB, ierr)
      call MatDestroy(XX, ierr)


      call MatCreateDense(psubcomm, nrow_loc_subcom, PETSC_DECIDE, PETSC_DECIDE, nrest, &
        PETSC_NULL_SCALAR, XX, ierr)

      call MatCreateDense(psubcomm, nrow_loc_subcom, PETSC_DECIDE, PETSC_DECIDE, nrest, &
        PETSC_NULL_SCALAR, BB, ierr)

      call MatAssemblyBegin(BB, MAT_FINAL_ASSEMBLY, ierr)
      call MatAssemblyEnd(BB, MAT_FINAL_ASSEMBLY, ierr)
      call MatAssemblyBegin(XX, MAT_FINAL_ASSEMBLY, ierr)
      call MatAssemblyEnd(XX, MAT_FINAL_ASSEMBLY, ierr)


      do i = 1, nrest
        icol = ioff + i - 1
        if ((icol(1) .ge. nrow_loc1) .and. (icol(1) .lt. nrow_loc2)) then
          j = icol(1)
          call MatGetRow(B,j,nz,jrows,p_batch,ierr)
          icol = i - 1
          call MatSetValues(BB, nz, jrows, 1, icol, p_batch(1:nz), INSERT_VALUES, ierr)
          call MatRestoreRow(B,j,nz,jrows,p_batch,ierr)
        end if
      end do
      call MatAssemblyBegin(BB, MAT_FINAL_ASSEMBLY, ierr) ! I think the values are cached.
      call MatAssemblyEnd(BB, MAT_FINAL_ASSEMBLY, ierr)

      call MatMatSolve(F, BB, XX, ierr)
      if (ierr.ne.0) then
        write (errormsg, fmt='(A,i8)') "solver rest error ", ierr
        call error()
      end if

      call MatDenseGetArrayF90(XX,p_xx,ierr)
      do i = 1, nrest
        icol = ioff + i - 1
        call MatSetValues(X, nrow_loc_subcom, irows_loc_subcom, 1, icol ,&
          p_xx(1:nrow_loc_subcom,i), INSERT_VALUES, ierr)
      end do
      call MatDenseRestoreArrayF90(XX,p_xx,ierr)

    end if

    call MatAssemblyBegin(X, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(X, MAT_FINAL_ASSEMBLY, ierr)




    call system_clock(countf)
    timing_matmat_solve(1) = timing_matmat_solve(1) + real(countf - counti, 8)/real(count_rate, 8)

    call MatDestroy(B, ierr)
    call MatDestroy(BB, ierr)
    call MatDestroy(XX, ierr)
    call MatDestroy(A, ierr)

    call KSPDestroy(ksp,ierr)

!~     call dump_nonzero_structure(X, "X_subcom.dat",PETSC_COMM_WORLD, 0)
!~     write(0,*) "------finished------"
!~     call MPI_BARRIER(PETSC_COMM_WORLD,ierr)

!~     call slepcfinalize(ierr)
!~     stop
  end subroutine petsc_solve_sparse_input_onsubcoms

  subroutine petsc_solve_direct(f, b, x)
#include <petsc/finclude/petsc.h>
    use petscmat
    use timing
    implicit none

    Mat :: f, b, x

    integer :: ierr
    integer(8) :: counti, count_rate, countf
    character(256) :: stype

    call system_clock(counti, count_rate)
    call MatMatSolve(f, b, x, ierr)
    call system_clock(countf)
    call MatGetType(x, stype, ierr)
    if (trim(adjustl(stype)) .eq. "elemental") then
      timing_matmat_solve(2) = timing_matmat_solve(2) + real(countf - counti, 8)/real(count_rate, 8)
    else
      timing_matmat_solve(1) = timing_matmat_solve(1) + real(countf - counti, 8)/real(count_rate, 8)
    end if

  end subroutine petsc_solve_direct

  subroutine petsc_solve_by_col_iter(a, b_in, x, matsolvertype)
#include <petsc/finclude/petsc.h>
    use petsc
    use petsc_mod
    use globals, only: inode, nprocs, inode_sub
    implicit none

    Mat :: a, b_in, x
    MatSolverType :: matsolvertype

    Vec :: b, xx
    Mat :: f
    integer :: ierr, irow, ngr, ngc, nl1, nl2, ii(1), nloc, ncols, j, icol
    integer, allocatable :: cols(:), idx(:), idy(:), rows(:)
    PetscScalar, pointer :: xx_v(:), b_mat_p(:, :), b_vec_p(:)
    PetscScalar, allocatable :: vals(:)
    PetscReal :: norm
    PetscScalar :: p_scale
    PetscErrorCode :: p_ierr
    PC :: pc
    KSP :: ksp
    PetscViewerAndFormat vf
    IS :: is1
    integer(8) :: counti, count_rate, countf

    call MatSetOption(x, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE, ierr)
    call MatSetOption(x, MAT_IGNORE_OFF_PROC_ENTRIES, PETSC_TRUE, ierr)

    call MatCreateVecs(a, b, PETSC_NULL_VEC, ierr)
    call MatCreateVecs(a, xx, PETSC_NULL_VEC, ierr)

    call MatGetSize(b_in, ngr, ngc, ierr)

    call MatGetOwnershipRange(b_in, nl1, nl2, ierr)

    nloc = nl2 - nl1

    allocate (rows(nloc), vals(nloc))

    j = 0
    do icol = nl1, nl2 - 1
      j = j + 1
      rows(j) = icol
    end do

!~     do icol=1,ngr
!~       rows(icol)=icol
!~     end do

    call VeczeroEntries(b, ierr)

    call KSPCreate(PETSC_COMM_WORLD, ksp, ierr)
    call KSPSetOperators(ksp, a, a, ierr)
    call KSPSetType(ksp, KSPPREONLY, ierr)

    call KSPGetPC(ksp, pc, ierr)
    call PCSetType(pc, PCLU, ierr)
    call PCFactorSetMatSolverType(pc, matsolvertype, ierr)
!~     call PCFactorSetUpMatSolverType(pc,ierr)
!~     call PCFactorGetMatrix(pc,f,ierr)

    call KSPSetFromOptions(ksp, ierr)

    call KSPSetUp(ksp, ierr)
!~     call KSPSetInitialGuessNonzero(ksp,PETSC_TRUE,ierr)

    call VecGetArrayF90(xx, xx_v, p_ierr)
    call VecGetArrayF90(b, b_vec_p, p_ierr)

    do icol = 0, ngc - 1

      ii = icol

      call MatGetValues(b_in, nloc, rows, 1, ii, b_vec_p(1:nloc), ierr)
!~       call system_clock(counti, count_rate)
      call KSPSolve(ksp, b, xx, ierr)
!~       call system_clock(countf)
!~       call KSPGetResidualNorm(ksp,norm,ierr)
!~       if (inode.eq.0) write(0,fmt='(3i8,2e24.12)') inode,icol,ngc,norm,real(countf - counti, 8)/real(count_rate, 8)

      call MatSetValues(x, nloc, rows, 1, ii, xx_v, INSERT_VALUES, ierr)

    end do

    call VecRestoreArrayF90(b, b_vec_p, p_ierr)
    call VecRestoreArrayF90(xx, xx_v, p_ierr)
    call MatAssemblyBegin(x, MAT_FINAL_ASSEMBLY, ierr) !is this call necessary ?
    call MatAssemblyEnd(x, MAT_FINAL_ASSEMBLY, ierr)

    call VecDestroy(b, ierr)
    call VecDestroy(xx, ierr)
    call KSPDestroy(ksp, ierr)
    call MatSetOption(x, MAT_IGNORE_OFF_PROC_ENTRIES, PETSC_FALSE, ierr)

  end subroutine petsc_solve_by_col_iter

  subroutine petsc_invert(a, x, matsolvertype, mattype)

    use petsc_mod
    implicit none

    Mat :: a, x
    MatSolverType :: matsolvertype
    MatType :: mattype

    Mat :: b, f
    integer :: ierr

    call petsc_get_densemat(a, b, mattype)
    call MatScale(b, p_zero, ierr)
    call MatShift(b, p_one, ierr)
    call petsc_lu_fac(a, f, matsolvertype)
    call petsc_solve_direct(f, b, x)
    call MatDestroy(b, ierr)
    call MatDestroy(f, ierr)

  end subroutine petsc_invert

  subroutine petsc_get_sqrt_mat(p_A, p_B, p_sqrtA)
#include <petsc/finclude/petscmat.h>
    use petscmat
    use petsc_mod
    use slepc_mod
    use globals, only: inode, nprocs, mattype_dense, mattype_surf, mattype_cc, mattype_sparse
    use error_handler
    implicit none

    Mat :: p_A, p_B, p_sqrtA

    Mat :: p_EV, p_tmp1, p_tmp2, p_tmp3, p_diag, p_tmp4
    Vec :: p_ew
    PetscScalar, pointer :: pp_ew(:)
    integer :: nev, ierr, nl1, nl2, i, j
    PetscReal :: norm

    call petsc_get_densemat(p_A, p_ev, mattype_dense)
    call MatCreateVecs(p_ev, p_ew, PETSC_NULL_vec, ierr)
    call MatDuplicate(p_A, MAT_COPY_VALUES, p_tmp4, ierr)
    call MatConvert(p_tmp4,mattype_sparse, MAT_INPLACE_MATRIX,p_tmp4,ierr)
    call VecGetSize(p_ew, nev, ierr)

    call diag_mat2(p_tmp4, p_B, p_ew, p_EV, nev, 1d-12)
!~    call diag_mat(p_A,p_B,p_ew,p_EV,nev)

    call VecGetOwnershipRange(p_ew, nl1, nl2, ierr)
!~     write(pstr_out,fmt='(A,3i8)') "ew: ",inode,nl1,nl2 ;  call petsc_print(pstr_out)

    j = 0
    call VecGetArrayReadF90(p_ew, pp_ew, ierr)

    do i = nl1, nl2 - 1

      j = j + 1
      if ((real(pp_ew(j)) .lt. 0d0) .and. 1 .eq. 2) then !.or.(abs(aimag(pp_ew(j))).gt.1d-6)) then
        if (real(pp_ew(j)) + 1d-6 .lt. 0d0) then
          write (errormsg, fmt='(A,2i8,2e24.12)') "EW error for SQRT ", i, j, pp_ew(j)
          call error()
        else
          pp_ew(j) = 0d0
        end if
      end if
      pp_ew(j) = zsqrt(pp_ew(j))
    end do

    call VecRestoreArrayReadF90(p_ew, pp_ew, ierr)

    call petsc_get_a_with_b_c(p_diag, p_A, p_A, mattype_cc)
    call MatSetOption(p_diag, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE, ierr)

    call MatDiagonalSet(p_diag, p_ew, INSERT_VALUES, ierr)

!~     call MatMatMult(p_B,p_ev,MAT_INITIAL_MATRIX,PETSC_DEFAULT_REAL,p_tmp1,ierr)

!~     call MatMatMult(p_tmp1,p_diag,MAT_INITIAL_MATRIX,PETSC_DEFAULT_REAL,p_tmp2,ierr)

!~     call MatTranspose(p_tmp1, MAT_INPLACE_MATRIX,p_tmp1,ierr) ; CHKERRQ(ierr)
!~     call MatConjugate(p_tmp1,ierr)

!~     call MatMatMult(p_tmp2,p_tmp1,MAT_INITIAL_MATRIX,PETSC_DEFAULT_REAL,p_tmp3,ierr)
!~     call MatMatMult(p_B,p_ev,MAT_INITIAL_MATRIX,PETSC_DEFAULT_REAL,p_tmp1,ierr)

    call MatMatMult(p_ev, p_diag, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, p_tmp2, ierr)

    call MatTranspose(p_ev, MAT_INPLACE_MATRIX, p_ev, ierr); CHKERRQ(ierr)
    call MatConjugate(p_ev, ierr)

    call MatMatMult(p_tmp2, p_ev, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, p_tmp1, ierr)

!~     call petsc_mat_info(p_tmp1,"p_tmp1     ",ierr)
!~     call petsc_mat_info(p_sqrtA,"p_sqrtA  ",ierr)
    call petsc_aXpY(p_sqrtA, p_tmp1, p_one, petsc_false, .false.)
!~     call MatCopy(p_tmp1,p_sqrtA,DIFFERENT_NONZERO_PATTERN,ierr)

!~     call petsc_mat_info(p_sqrtA,"p_sqrtA  ",ierr)

!~     call MatNorm(p_sqrtA,NORM_FROBENIUS,norm,ierr)
!~     write(pstr_out,fmt='(A,e24.12)') "p_sqrtA: ",norm ; call petsc_print_master()
!~     call MatNorm(p_tmp3,NORM_FROBENIUS,norm,ierr)
!~     write(pstr_out,fmt='(A,e24.12)') "p_tmp3: ",norm ; call petsc_print_master()

    call MatDestroy(p_tmp1, ierr)
    call MatDestroy(p_tmp2, ierr)
!~     call MatDestroy(p_tmp3,ierr)
    call MatDestroy(p_ev, ierr)
    call VecDestroy(p_ew, ierr)

  end subroutine petsc_get_sqrt_mat


  subroutine petsc_call_solver(p_A, p_B, p_X, matsolvertype, isolver)
#include <petsc/finclude/petscmat.h>
    use petscmat
    use petsc_mod
    use globals, only : nsim_rhs
    implicit none

    Mat :: p_A, p_B, p_X
    MatSolverType :: matsolvertype
    integer :: isolver

    if (isolver .eq. 1) then
!~       call petsc_solve_onsubcoms(p_A, p_B, p_X, matsolvertype)
      call petsc_solve_sparse_input_onsubcoms(p_A, p_B, p_X, nsim_rhs, matsolvertype)
    else if (isolver .eq. 2) then
      call petsc_solve_cols(p_A, p_B, p_X, nsim_rhs, matsolvertype)
    else if (isolver .eq. 3) then
      call petsc_solve_by_col_iter(p_A, p_B, p_X, matsolvertype)
    else if (isolver .eq. 4) then
      call petsc_solve_cols_sparse_input(p_A, p_B, p_X, nsim_rhs, matsolvertype)
    else
      write (pstr_out, fmt='(A,i8)') "unknown solver mode", isolver; call petsc_print_master()
    end if

  end subroutine petsc_call_solver
  
  subroutine get_mat_norm(p_A, norm)
#include <petsc/finclude/petscmat.h>
    use petscmat
    
    Mat :: p_A  
    PetscReal :: norm
    integer :: ierr
    
    call MatNorm(p_A, NORM_FROBENIUS, norm, ierr)
    
  end subroutine get_mat_norm

end module petsc_wrapper
