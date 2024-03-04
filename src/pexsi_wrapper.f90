!~ use psellinv from PEXSI library to obtain only the necessary (nz) elements of
!~ the density matrix  D_ij=int_dE G_ij(E) , from G=[H-SE-Sigma]^-1
!~ For this we need to copy the nz elements from PETSC to a temporary array, along
!~ with the CSC information.
!~ PEXSI returns (Ainv^T)_ij so to get the correct G we must set options%transpose=1.
!~ If inv_sel is used for something else in the future, one has to keep this in mind.

!~ PEXSI uses 1 indexing and PETSC 0 indexing !!!!
module pexsi_wrapper
  implicit none

  contains

    subroutine inv_sel(p_A, p_Ainv)
#include <petsc/finclude/petscmat.h>
      use petsc
      use kinds
      use petsc_mod
      use error_handler
      use f_ppexsi_interface
      use iso_c_binding
      use globals, only: inode, nprocs

      implicit none

      Mat :: p_A, p_Ainv

      integer(c_intptr_t) :: plan
      integer(c_int), allocatable ::  colptrloc(:), rowindloc(:)
      integer(c_int):: nprow, npcol, npSymbFact, outputFileIndex
      type(f_ppexsi_options) :: options
      PetscScalar, allocatable :: Aloc(:), Ainvloc(:), Sloc(:)
      integer :: nnz_local, nrow, ncol, ierr, irow, icol, nzcol, ioff, nl1, &
        nl2, il, nnz_global, nrow_local, ncol_local, ioff1, ioff2, jj(1), nz_found



!~       call MatTranspose(p_A,  MAT_INPLACE_MATRIX, p_A, ierr)

      call MatGetSize(p_A, nrow, ncol, ierr)
      nnz_local = petsc_get_nnz(p_A, MAT_LOCAL)
      call MPI_BARRIER(PETSC_COMM_WORLD, ierr)
      call MatGetLocalSize(p_A, nrow_local, ncol_local, ierr)

      call MatSetOption(p_Ainv, MAT_IGNORE_OFF_PROC_ENTRIES, PETSC_TRUE, ierr)
      call MatSetOption(p_Ainv, MAT_NO_OFF_PROC_ZERO_ROWS, PETSC_TRUE, ierr)
      call MatSetOption(p_Ainv, MAT_NO_OFF_PROC_ENTRIES, PETSC_TRUE, ierr)
      call MatSetOption(p_Ainv, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE, ierr)
!~       pstr_out = "XX"; call petsc_print_master()
!~       call petsc_mat_info(p_A, "p_A ", ierr)
!~       call petsc_mat_info(p_Ainv, "p_Ainv ", ierr)

      allocate(colptrloc(nrow_local+1), rowindloc(nnz_local), Aloc(nnz_local), &
        Ainvloc(nnz_local), Sloc(1), stat = ierr)
      if (ierr .ne. 0) then
        write (errormsg, fmt='(A,i8)') "allocation error ", ierr
        call error()
      end if

      call MatGetOwnershipRange(p_A, nl1, nl2, ierr)

      colptrloc(1) = 1 ! assuming a non-zero row
      irow = 0
      nz_found = 0
      do il = nl1, nl2 - 1
        irow = irow + 1
        ioff = colptrloc(irow)
        call MatGetRow(p_A, il, nzcol, rowindloc(ioff:nnz_local), Aloc(ioff:nnz_local), ierr)
        nz_found = nz_found + nzcol
        colptrloc(irow + 1) = colptrloc(irow)  + nzcol
        call MatRestoreRow(p_A, il, nzcol, rowindloc(ioff:nnz_local), Aloc(ioff:nnz_local), ierr)
        nz_found = nz_found + nzcol
        if (nz_found.eq.nnz_local) exit
      end do
      rowindloc = rowindloc + 1


      if (inode .eq. 0) then
        outputFileIndex = 0
      else
        outputFileIndex = -666
      end if

!~       nprow = nprocs
!~       npcol = 1
      call multiplicative_partition(nprocs, npcol, nprow)

      plan = f_ppexsi_plan_initialize(PETSC_COMM_WORLD, nprow, npcol, outputFileIndex, ierr)
      if (ierr .ne. 0) then
        write (errormsg, fmt='(A,i8)') "f_ppexsi_plan_initialize error ", ierr
        call error()
      end if

      call f_ppexsi_set_default_options(options)

      options%ordering=1
      options%rowOrdering=0
      options%transpose=1
      options%verbosity=0


      nnz_global = petsc_get_nnz(p_A, MAT_GLOBAL_SUM)

      call f_ppexsi_load_complex_hs_matrix(plan, options, nrow, nnz_global, &
        nnz_local, nrow_local, colptrloc, rowindloc, Aloc, 1, Sloc, ierr)

      if (ierr .ne. 0) then
        write (errormsg, fmt='(A,i8)') "f_ppexsi_load_complex_hs_matrix error ", ierr
        call error()
      end if

      call f_ppexsi_symbolic_factorize_complex_unsymmetric_matrix(plan, options, &
        Aloc, ierr)
      if (ierr .ne. 0) then
        write (errormsg, fmt='(A,i8)') &
          "f_ppexsi_symbolic_factorize_complex_unsymmetric_matrix error ", ierr
        call error()
      end if

      call f_ppexsi_selinv_complex_unsymmetric_matrix(plan, options, Aloc, Ainvloc, ierr)
      if (ierr .ne. 0) then
        write (errormsg, fmt='(A,i8)') &
          "f_ppexsi_selinv_complex_unsymmetric_matrix error ", ierr
        call error()
      end if

      call f_ppexsi_plan_finalize(plan, ierr)
      if (ierr .ne. 0) then
        write (errormsg, fmt='(A,i8)') "f_ppexsi_plan_finalize error ", ierr
        call error()
      end if

      irow = 0
      rowindloc = rowindloc - 1

      do il = nl1, nl2 - 1
        irow = irow + 1
        nzcol = colptrloc(irow + 1) - colptrloc(irow)
        ioff1 = colptrloc(irow)
        ioff2 = ioff1 + nzcol - 1
        jj = il
        call MatSetValues(p_Ainv, 1, jj, nzcol, rowindloc(ioff1:ioff2), &
          Ainvloc(ioff1:ioff2), INSERT_VALUES, ierr)
      end do
      call MatAssemblyBegin(p_Ainv, MAT_FINAL_ASSEMBLY, ierr)
      call MatAssemblyEnd(p_Ainv, MAT_FINAL_ASSEMBLY, ierr)


    end subroutine inv_sel

    subroutine multiplicative_partition(n, i1, i2)
      implicit none

      integer :: n, i1 ,i2, i

      i1 = n
      i2 = 1

      do i = 1, n / 2
        if (mod(n, i) .eq. 0) then
          if (((n / i) .le. i1).and. ( n / i .gt. i2)) then
            i1=n / i
            i2=n / i1
          else
            exit
          end if
        end if
      end do

    end subroutine multiplicative_partition

end module pexsi_wrapper
