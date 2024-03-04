module surf_gf_mod
  implicit none
contains

!   subroutine decimate2(w,t1,eps,lr)
! exploiting block structure to obtain all matricies in one shot. that is making one LU factorization,
! calculation X=A^-1*B by solving the linear equation.

  subroutine surf_gf3(z, i, c)
    implicit none
    double complex :: z
    integer :: i
    character(*) :: c

  end subroutine surf_gf3

  subroutine load_surf_gf(energy, ik, lr, p_gf, l_loaded)
#include <petsc/finclude/petscmat.h>
    use petscmat
    use kinds
    use petsc_mod
    use globals, only: l_ionode, kp_r

    implicit none

    complex(dp) :: energy
    character(*) :: lr
    integer :: ik
    Mat :: p_gf
    logical :: l_loaded

    integer :: ierr
    character(strln) :: gf_file, s_kx, s_ky, s_real, s_imag
    logical :: l_exist
    PetscScalar :: pp
    PetscReal :: pnorm

    l_loaded = .false.
  
    write (s_kx, fmt='(e16.9)') kp_r(1, ik)
    write (s_ky, fmt='(e16.9)') kp_r(2, ik)
    write (s_real, fmt='(e16.9)') real(energy, 8)
    write (s_imag, fmt='(e16.9)') aimag(energy)
    gf_file = trim(adjustl(s_real))//"_"//trim(adjustl(s_imag))
    gf_file = trim(adjustl(gf_file))//"_"//trim(adjustl(s_kx))//"_"//trim(adjustl(s_ky))
    gf_file = "./reaktor/gf_"//trim(adjustl(lr))//"_"//trim(adjustl(gf_file))//".dat"

    inquire (file=trim(gf_file), exist=l_exist)

    if (.not. l_exist) return

    call petsc_mat_direct_load(p_gf, gf_file, ierr)
    l_loaded = 0 .eq. ierr

  end subroutine load_surf_gf

  subroutine save_surf_gf(energy, ik, lr, p_gf)
#include <petsc/finclude/petscmat.h>
    use petscmat
    use kinds
    use petsc_mod
    use globals, only: l_ionode, kp_r

    implicit none

    complex(dp) :: energy
    character(*) :: lr
    integer :: ik
    Mat :: p_gf

    integer :: ierr
    character(strln) :: gf_file, s_kx, s_ky, s_real, s_imag
    logical :: l_exist
    PetscScalar :: pp
    PetscReal :: pnorm

    write (s_kx, fmt='(e16.9)') kp_r(1, ik)
    write (s_ky, fmt='(e16.9)') kp_r(2, ik)
    write (s_real, fmt='(e16.9)') real(energy, 8)
    write (s_imag, fmt='(e16.9)') aimag(energy)
    gf_file = trim(adjustl(s_real))//"_"//trim(adjustl(s_imag))
    gf_file = trim(adjustl(gf_file))//"_"//trim(adjustl(s_kx))//"_"//trim(adjustl(s_ky))
    gf_file = "./reaktor/gf_"//trim(adjustl(lr))//"_"//trim(adjustl(gf_file))//".dat"

    inquire (file=trim(gf_file), exist=l_exist)

    if (l_exist) return

    call petsc_mat_direct_save(p_gf, gf_file, ierr)

  end subroutine save_surf_gf

  subroutine decimate3_petsc(w1, w2, w3, w4, t1, t2, eps, lr)
#include <petsc/finclude/petscmat.h>
    use petscmat
    use petsc_mod
    use petsc_wrapper
    use globals, only: inode, matsolvertype_surf, mattype_surf
    use error_handler
    implicit none

    double precision, intent(in) :: eps

    Mat :: w1, w2, w3, w4, t1, t2

    MatStructure :: matstruct

    character(1), intent(in) :: lr

    Mat :: x1, x2, x3, x4, tmp, tmp2, tmp3, f
    MatReuse :: mat_reuse

    double precision :: conv, convold, conv_norm
    PetscReal :: norm
    integer :: n, step, ierr

    call petsc_get_densemat(w4, x1, mattype_surf)
    call petsc_get_densemat(w4, x2, mattype_surf)
    call petsc_get_densemat(w4, x3, mattype_surf)
    call petsc_get_densemat(w4, x4, mattype_surf)
!~     call petsc_get_densemat(w4,tmp,mattype_surf)
!~     call petsc_get_densemat(w4,tmp2,mattype_surf)

    step = 0
    conv = 666d0
    mat_reuse = MAT_INITIAL_MATRIX

    do

      step = step + 1
!~       if (inode.eq.0) write(0,*) "step ",step,conv

      call petsc_lu_fac(w4, f, matsolvertype_surf)  !f=lu(w4)
      call petsc_solve_direct(f, w3, x1)
      call petsc_solve_direct(f, t1, x2)
      call MatDestroy(f, ierr)

      call MatMatMult(w2, x1, mat_reuse, PETSC_DEFAULT_REAL, tmp, ierr)
      call MatAYPX(tmp, p_minus1, w1, SAME_NONZERO_PATTERN, ierr)

      call petsc_lu_fac(tmp, f, matsolvertype_surf)
      call petsc_solve_direct(f, t2, x3)
      call petsc_solve_direct(f, w2, x4)
      call MatDestroy(f, ierr)

      call MatMatMult(t1, x3, mat_reuse, PETSC_DEFAULT_REAL, tmp2, ierr)
      if (lr .eq. "r") then
        call MatNorm(tmp2, NORM_FROBENIUS, conv_norm, ierr)
      end if
      call MatAXPY(w4, p_minus1, tmp2, SAME_NONZERO_PATTERN, ierr) ! W4

      call MatMatMult(x4, x2, MAT_REUSE_MATRIX, PETSC_DEFAULT_REAL, tmp2, ierr)

      call MatMatMult(t1, tmp2, mat_reuse, PETSC_DEFAULT_REAL, tmp3, ierr)
      call MatCopy(tmp3, t1, SAME_NONZERO_PATTERN, ierr)

      call MatMatMult(x1, tmp2, MAT_REUSE_MATRIX, PETSC_DEFAULT_REAL, tmp, ierr) ! tmp=x1*x4*x2
      call MatAXPY(tmp, p_one, x2, SAME_NONZERO_PATTERN, ierr)

      call MatMatMult(t2, tmp, MAT_REUSE_MATRIX, PETSC_DEFAULT_REAL, tmp2, ierr)
      if (lr .eq. "l") then
        call MatNorm(tmp2, NORM_FROBENIUS, conv_norm, ierr)
      end if
      call MatAXPY(w1, p_minus1, tmp2, SAME_NONZERO_PATTERN, ierr)

      call MatMatMult(x1, x3, MAT_REUSE_MATRIX, PETSC_DEFAULT_REAL, tmp2, ierr)

      call MatMatMult(t2, tmp2, MAT_REUSE_MATRIX, PETSC_DEFAULT_REAL, tmp3, ierr)
      call MatCopy(tmp3, t2, SAME_NONZERO_PATTERN, ierr)

      if (lr .eq. "r") then
        convold = conv
        call MatNorm(w4, NORM_FROBENIUS, norm, ierr)
        conv = conv_norm/norm
      else if (lr .eq. "l") then
        convold = conv
        call MatNorm(w1, NORM_FROBENIUS, norm, ierr)
        conv = conv_norm/norm
      end if

      if (conv .le. eps) exit

      if ((step .ge. 1000) .or. (isnan(conv))) then
        write (errormsg, fmt='(A,i8,2e24.12,2A)') "step limit ", step, conv, convold, lr
        call error()
      end if

      mat_reuse = MAT_REUSE_MATRIX

    end do

    call MatDestroy(x1, ierr)
    call MatDestroy(x2, ierr)
    call MatDestroy(x3, ierr)
    call MatDestroy(x4, ierr)
    call MatDestroy(tmp, ierr)
    call MatDestroy(tmp2, ierr)
    call MatDestroy(tmp3, ierr)

  end subroutine decimate3_petsc

  subroutine decimate3_v2_petsc(w1, w2, w3, w4, t1, t2, eps, lr)
#include <petsc/finclude/petscmat.h>
    use petscmat
    use petsc_mod
    use petsc_wrapper
    use globals, only: inode, matsolvertype_surf, mattype_surf
    use error_handler
    implicit none

    double precision, intent(in) :: eps

    Mat :: w1, w2, w3, w4, t1, t2

    MatStructure :: matstruct

    character(1), intent(in) :: lr

    Mat :: x1, x2, x3, x4, tmp, f, oneMat, invW4, &
      WW1, WW2, WW3, WW4, invW4W3
    MatReuse :: mat_reuse

    double precision :: conv, convold, conv_norm
    PetscReal :: norm
    integer :: n, step, ierr

!~     call petsc_get_densemat(w4, x1, mattype_surf)
!~     call petsc_get_densemat(w4, x2, mattype_surf)
!~     call petsc_get_densemat(w4, x3, mattype_surf)
!~     call petsc_get_densemat(w4, x4, mattype_surf)

    call petsc_get_densemat(w4, oneMat, mattype_surf)
    call petsc_get_densemat(w4, invw4, mattype_surf)
    call petsc_get_densemat(w4, WW1, mattype_surf)
!~     call petsc_get_densemat(w4, WW2, mattype_surf)
!~     call petsc_get_densemat(w4, WW3, mattype_surf)
!~     call petsc_get_densemat(w4, WW4, mattype_surf)
!~     call petsc_get_densemat(w4, tmp, mattype_surf)

    call MatZeroEntries(oneMat,ierr)
    call MatShift(oneMat, p_one, ierr)
!~     call petsc_get_densemat(w4,tmp,mattype_surf)
!~     call petsc_get_densemat(w4,tmp2,mattype_surf)

    step = 0
    conv = 666d0
    mat_reuse = MAT_INITIAL_MATRIX
!~     mat_reuse = MAT_REUSE_MATRIX

    do

      step = step + 1
!~       if (inode.eq.0) write(0,*) "step ",step,conv
      call petsc_lu_fac(w4, f, matsolvertype_surf)
      call petsc_solve_direct(f, oneMat, invw4) ! invw4=W4^-1
      call MatDestroy(f, ierr)
!~       call MatMatMatMult(w2, invw4, w3, mat_reuse, PETSC_DEFAULT_REAL, WW2, ierr) ! WW2=W2*W4^-1*W3
      call MatMatMult(invw4, w3, mat_reuse, PETSC_DEFAULT_REAL, invW4W3, ierr) ! WW2=W2*W4^-1*W3
      call MatMatMult(w2, invW4W3, mat_reuse, PETSC_DEFAULT_REAL, WW2, ierr) ! WW2=W2*W4^-1*W3

      call MatAYPX(WW2, p_minus1, w1, SAME_NONZERO_PATTERN, ierr) ! WW2=W1-W2*W4^-1*W3
      call petsc_lu_fac(WW2, f, matsolvertype_surf)
      call petsc_solve_direct(f, oneMat, WW1)  !WW1=(W1-W2*W4^-1*W3)^-1
      call MatDestroy(f, ierr)


!~       call MatMatMatMult(WW1, W2, invW4, MAT_REUSE_MATRIX, PETSC_DEFAULT_REAL, WW2, ierr) ! WW2=WW1*W2*invW4
      call MatMatMult(WW1, W2, mat_reuse, PETSC_DEFAULT_REAL, tmp, ierr) ! WW2=WW1*W2*invW4
      call MatMatMult(tmp, invW4, MAT_REUSE_MATRIX, PETSC_DEFAULT_REAL, WW2, ierr) ! WW2=WW1*W2*invW4

!~       call MatMatMatMult(invW4, W3, WW1, mat_reuse, PETSC_DEFAULT_REAL, WW3, ierr) ! WW3=invW4*W3*WW1
!~       call MatMatMult(invW4, W3, MAT_REUSE_MATRIX, PETSC_DEFAULT_REAL, tmp, ierr) ! WW3=invW4*W3*WW1
      call MatMatMult(invW4W3, WW1, mat_reuse, PETSC_DEFAULT_REAL, WW3, ierr) ! WW3=invW4*W3*WW1

!~       call MatMatMatMult(invW4, W4, WW3, mat_reuse, PETSC_DEFAULT_REAL, WW4, ierr) ! WW4=invW4*W3*W2
!~       call MatMatMult(invW4, W3, MAT_REUSE_MATRIX, PETSC_DEFAULT_REAL, tmp, ierr) ! WW4=invW4*W3*W2
      call MatMatMult(invW4W3, WW2, mat_reuse, PETSC_DEFAULT_REAL, WW4, ierr) ! WW4=invW4*W3*W2
      call MatAYPX(WW4, p_one, invW4, SAME_NONZERO_PATTERN, ierr) ! WW4=invW4+invW4*W3*W2

      ! W1 update
!~       call MatMatMatMult(t2, WW4, t1, mat_reuse, PETSC_DEFAULT_REAL, invW4, ierr) ! invW4=t2*WW4*t1
      call MatMatMult(t2, WW4, MAT_REUSE_MATRIX, PETSC_DEFAULT_REAL, tmp, ierr) ! invW4=t2*WW4*t1
      call MatMatMult(tmp, t1, mat_reuse, PETSC_DEFAULT_REAL, WW4, ierr) ! invW4=t2*WW4*t1
      if (lr .eq. "r") then
        call MatNorm(WW4, NORM_FROBENIUS, conv_norm, ierr)
      end if
      call MatAXPY(W1, p_minus1, WW4, SAME_NONZERO_PATTERN, ierr) ! W1=W1-t2*WW4*t1

      ! W4 update
!~       call MatMatMatMult(t1, WW1, t2, MAT_REUSE_MATRIX, PETSC_DEFAULT_REAL, invW4, ierr) ! invW4=t1*WW1*t2
      call MatMatMult(t1, WW1, MAT_REUSE_MATRIX, PETSC_DEFAULT_REAL, tmp, ierr) ! invW4=t1*WW1*t2
      call MatMatMult(tmp, t2, MAT_REUSE_MATRIX, PETSC_DEFAULT_REAL, WW4, ierr) ! invW4=t1*WW1*t2
      if (lr .eq. "l") then
        call MatNorm(WW4, NORM_FROBENIUS, conv_norm, ierr)
      end if
      call MatAXPY(W4, p_minus1, WW4, SAME_NONZERO_PATTERN, ierr) ! W4=W4-t1*WW1*t2

      !t1 update
!~       call MatMatMatMult(t1, WW2, t1, MAT_REUSE_MATRIX, PETSC_DEFAULT_REAL, invW4, ierr) ! invW4=t1*WW2*t1
      call MatMatMult(t1, WW2, MAT_REUSE_MATRIX, PETSC_DEFAULT_REAL, tmp, ierr) ! invW4=t1*WW2*t1
      call MatMatMult(tmp, t1, MAT_REUSE_MATRIX, PETSC_DEFAULT_REAL, WW2, ierr) ! invW4=t1*WW2*t1
      call MatCopy(WW2, t1, SAME_NONZERO_PATTERN, ierr)

      !t2 update
!~       call MatMatMatMult(t2, WW3, t2, MAT_REUSE_MATRIX, PETSC_DEFAULT_REAL, invW4, ierr) ! invW4=t2*WW3*t2
      call MatMatMult(t2, WW3, MAT_REUSE_MATRIX, PETSC_DEFAULT_REAL, tmp, ierr) ! invW4=t2*WW3*t2
      call MatMatMult(tmp, t2, MAT_REUSE_MATRIX, PETSC_DEFAULT_REAL, WW3, ierr) ! invW4=t2*WW3*t2
      call MatCopy(WW3, t2, SAME_NONZERO_PATTERN, ierr)


      if (lr .eq. "l") then
        convold = conv
        call MatNorm(w4, NORM_FROBENIUS, norm, ierr)
        conv = conv_norm/norm
      else if (lr .eq. "r") then
        convold = conv
        call MatNorm(w1, NORM_FROBENIUS, norm, ierr)
        conv = conv_norm/norm
      end if

      if (conv .le. eps) exit

      if ((step .ge. 1000) .or. (isnan(conv))) then
        write (errormsg, fmt='(A,i8,2e24.12,2A)') "step limit ", step, conv, convold, lr
        call error()
      end if

      mat_reuse = MAT_REUSE_MATRIX

    end do


    call MatDestroy(oneMat, ierr)
    call MatDestroy(invw4, ierr)
    call MatDestroy(WW1, ierr)
    call MatDestroy(WW2, ierr)
    call MatDestroy(WW3, ierr)
    call MatDestroy(WW4, ierr)
    call MatDestroy(tmp, ierr)

  end subroutine decimate3_v2_petsc

  subroutine surf_gf3_petsc(z, ik, lr)
#include <petsc/finclude/petscmat.h>
    use petscmat
    use petsc_mod
    use petsc_wrapper
    use globals
    use mathstuff
    use kinds
    use timing
    use error_handler
    implicit none

    character(1) :: lr
    complex(dp) :: z
    integer :: ik

    double complex, parameter :: zone = 1d0

    Mat :: w1, w2, w3, w4, t1, t2, f, p_tmp

    PetscScalar :: zdos
    PetscReal :: norm

    integer, allocatable :: ipiv(:)
    integer :: ierr, imu1
    integer(8) :: counti, count_rate, countf

    logical :: l_loaded

    call system_clock(counti, count_rate)

    select case (lr)
    case ("r")

!~       if (l_loadsavegf.and.(.not.l_k_on_demand)) then
      if (l_loadsavegf) then
        call load_surf_gf(z, ik, "r", p_grrr, l_loaded)
        if (l_loaded) then
          call MatGetTrace(p_grrr, zdos, ierr)
          zdosr(1) = zdosr(1) + zdos*wkp_r(ik)
          return
        end if
      end if
!~       call petsc_get_densemat(p_h00k_r(ik),w1,mattype_surf)
!~       call petsc_get_densemat(p_h00k_r(ik),w2,mattype_surf)
!~       call petsc_get_densemat(p_h00k_r(ik),w3,mattype_surf)
!~       call petsc_get_densemat(p_h00k_r(ik),w4,mattype_surf)
!~       call petsc_get_densemat(p_h00k_r(ik),t1,mattype_surf)
!~       call petsc_get_densemat(p_h00k_r(ik),t2,mattype_surf)

      call MatDuplicate(p_h00k_r(ik), MAT_COPY_VALUES, w1, ierr)
      call MatScale(w1, p_minus1, ierr)
      call MatAXPY(w1, z, p_s00k_r(ik), SAME_NONZERO_PATTERN, ierr)
      call MatConvert(w1, mattype_surf, MAT_INPLACE_MATRIX, w1, ierr)

      call MatDuplicate(p_h01k_r(ik), MAT_COPY_VALUES, w2, ierr)
      call MatScale(w2, p_minus1, ierr)
      call MatAXPY(w2, z, p_s01k_r(ik), SAME_NONZERO_PATTERN, ierr)
      call MatConvert(w2, mattype_surf, MAT_INPLACE_MATRIX, w2, ierr)

      call MatDuplicate(p_h10k_r(ik), MAT_COPY_VALUES, w3, ierr)
      call MatScale(w3, p_minus1, ierr)
      call MatAXPY(w3, z, p_s10k_r(ik), SAME_NONZERO_PATTERN, ierr)
      call MatConvert(w3, mattype_surf, MAT_INPLACE_MATRIX, w3, ierr)

      call MatDuplicate(w1, MAT_COPY_VALUES, p_tmp, ierr)
      call MatDuplicate(w1, MAT_COPY_VALUES, w4, ierr)
      call MatDuplicate(w2, MAT_COPY_VALUES, t1, ierr)
      call MatDuplicate(w3, MAT_COPY_VALUES, t2, ierr)

      call decimate3_petsc(w1, w2, w3, w4, t1, t2, conv_dez, "r")

      call petsc_lu_fac(w4, f, matsolvertype_surf)
      call petsc_solve_direct(f, w3, t1)

      call MatDestroy(t2, ierr)
      call MatMatMult(w2, t1, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, t2, ierr)
      call MatAYPX(t2, p_minus1, p_tmp, SAME_NONZERO_PATTERN, ierr)

      call petsc_invert(t2, p_grrr, matsolvertype_surf, mattype_surf)

      if (l_loadsavegf) call save_surf_gf(z, ik, "r", p_grrr)

      call MatGetTrace(p_grrr, zdos, ierr)
      zdosr(1) = zdosr(1) + zdos*wkp_r(ik)

      call MatDestroy(w1, ierr)
      call MatDestroy(w2, ierr)
      call MatDestroy(w3, ierr)
      call MatDestroy(w4, ierr)
      call MatDestroy(t1, ierr)
      call MatDestroy(t2, ierr)
      call MatDestroy(f, ierr)
      call MatDestroy(p_tmp, ierr)

    case ("l")

!~       if (l_loadsavegf.and.(.not.l_k_on_demand)) then
      if (l_loadsavegf) then
        call load_surf_gf(z, ik, "l", p_gllr, l_loaded)
        if (l_loaded) then
          call MatGetTrace(p_gllr, zdos, ierr)
          zdosl(1) = zdosl(1) + zdos*wkp_l(ik)
          return
        end if
      end if

!~       call petsc_get_densemat(p_h00k_l(ik),w1,mattype_surf)
!~       call petsc_get_densemat(p_h00k_l(ik),w2,mattype_surf)
!~       call petsc_get_densemat(p_h00k_l(ik),w3,mattype_surf)
!~       call petsc_get_densemat(p_h00k_l(ik),w4,mattype_surf)
!~       call petsc_get_densemat(p_h00k_l(ik),t1,mattype_surf)
!~       call petsc_get_densemat(p_h00k_l(ik),t2,mattype_surf)

      call MatDuplicate(p_h00k_l(ik), MAT_COPY_VALUES, w1, ierr)
      call MatScale(w1, p_minus1, ierr)
      call MatAXPY(w1, z, p_s00k_l(ik), SAME_NONZERO_PATTERN, ierr)
      call MatConvert(w1, mattype_surf, MAT_INPLACE_MATRIX, w1, ierr)

      call MatDuplicate(p_h01k_l(ik), MAT_COPY_VALUES, w2, ierr)
      call MatScale(w2, p_minus1, ierr)
      call MatAXPY(w2, z, p_s01k_l(ik), SAME_NONZERO_PATTERN, ierr)
      call MatConvert(w2, mattype_surf, MAT_INPLACE_MATRIX, w2, ierr)

      call MatDuplicate(p_h10k_l(ik), MAT_COPY_VALUES, w3, ierr)
      call MatScale(w3, p_minus1, ierr)
      call MatAXPY(w3, z, p_s10k_l(ik), SAME_NONZERO_PATTERN, ierr)
      call MatConvert(w3, mattype_surf, MAT_INPLACE_MATRIX, w3, ierr)

      call MatDuplicate(w1, MAT_COPY_VALUES, p_tmp, ierr)
      call MatDuplicate(w1, MAT_COPY_VALUES, w4, ierr)
      call MatDuplicate(w2, MAT_COPY_VALUES, t1, ierr)
      call MatDuplicate(w3, MAT_COPY_VALUES, t2, ierr)

      call decimate3_petsc(w1, w2, w3, w4, t1, t2, conv_dez, "l")

      call petsc_lu_fac(w1, f, matsolvertype_surf)
      call petsc_solve_direct(f, w2, t1)

      call MatDestroy(t2, ierr)
      call MatMatMult(w3, t1, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, t2, ierr)
      call MatAYPX(t2, p_minus1, p_tmp, SAME_NONZERO_PATTERN, ierr)

      call petsc_invert(t2, p_gllr, matsolvertype_surf, mattype_surf)

      if (l_loadsavegf) call save_surf_gf(z, ik, "l", p_gllr)

      call MatGetTrace(p_gllr, zdos, ierr)
      zdosl(1) = zdosl(1) + zdos*wkp_r(ik)

      call MatDestroy(w1, ierr)
      call MatDestroy(w2, ierr)
      call MatDestroy(w3, ierr)
      call MatDestroy(w4, ierr)
      call MatDestroy(t1, ierr)
      call MatDestroy(t2, ierr)
      call MatDestroy(f, ierr)
      call MatDestroy(p_tmp, ierr)

    case default

      write (errormsg, fmt='(3A)') "caluclation surface GF for ", lr, " ... not !"
      call error()

    end select

    call system_clock(countf)
    timing_surf = timing_surf + real(countf - counti, 8)/real(count_rate, 8)

  end subroutine surf_gf3_petsc

end module surf_gf_mod
