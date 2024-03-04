module petsc_test
#include <slepc/finclude/slepceps.h>
  use petscmat
  use slepceps
  use slepc_mod
  use petsc_mod
  use petsc_wrapper
  implicit none

contains

  subroutine test_petsc()
    use globals
    use petsc_mod
    implicit none

    Mat :: a, c, f, y, e, d, g, h

    Vec :: b, x

    KSP :: ksp
    PC :: pc
    PetscScalar, pointer ::  xx_v(:), bb_v(:), mat_p(:, :)

    MatSolverType :: matsolvertype1
    MatType ::  mattype
    PetscReal :: norm, tol

    integer :: ierr, nl1, nl2, irow

!~     matsolvertype1=MATSOLVERMUMPS
!~     matsolvertype1=MATSOLVERPASTIX
!~     matsolvertype1=MATSOLVERSUPERLU_DIST
    matsolvertype1 = MATSOLVERMKL_CPARDISO

    call petsc_mat_info(p_h00_cc(0, 0), "p_h00_cc(0,0) ", ierr)

    call MatDuplicate(p_h00_cc(0, 0), MAT_COPY_VALUES, a, ierr)

    call MatGetOwnershipRange(a, nl1, nl2, ierr); CHKERRQ(ierr)

    call MatCreateVecs(a, b, PETSC_NULL_VEC, ierr)
    call MatCreateVecs(a, x, PETSC_NULL_VEC, ierr)

    irow = 0

    if ((irow .ge. nl1) .and. (irow .le. nl2 - 1)) then
      call VecGetArrayF90(b, xx_v, ierr)
      xx_v(irow - nl1 + 1) = 1d0
      call VecRestoreArrayF90(b, xx_v, ierr)
    end if

    call KSPCreate(PETSC_COMM_WORLD, ksp, ierr)
    call KSPSetOperators(ksp, a, a, ierr)
    call KSPSetType(ksp, KSPPREONLY, ierr)
    call KSPGetPC(ksp, pc, ierr)
    call PCSetType(pc, PCLU, ierr)
    call PCFactorSetMatSolverType(pc, matsolvertype1, ierr)
    call KSPSetUp(ksp, ierr)
    call PCFactorGetMatrix(pc, F, ierr)
    call KSPGetPC(ksp, pc, ierr)
    call KSPSolve(ksp, b, x, ierr)

    if ((irow .ge. nl1) .and. (irow .le. nl2 - 1)) then
      call VecGetArrayF90(b, xx_v, ierr)
      xx_v(irow - nl1 + 1) = 0d0
      call VecRestoreArrayF90(b, xx_v, ierr)
    end if

  end subroutine test_petsc

end module petsc_test
