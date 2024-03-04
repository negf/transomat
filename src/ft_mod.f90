module ft_mod
  implicit none

contains

! transfrom y(k)=sum_ij exp(k*R_ij)*x(i,j)
  subroutine fourier_trans(x, y, kp, nk, brav, maxi, maxj, same)
#include <petsc/finclude/petsc.h>
    use petscmat
    use petsc_mod
    use kinds
!~       use blas95
!~       use globals
    implicit none

    Mat, intent(inout) :: x(-maxi:maxi, -maxj:maxj)
    Mat, intent(inout) :: y(1:nk)

    integer :: maxi, maxj, nk
    integer :: same
    real(dp) :: brav(3, 3)
    real(dp) :: kp(:, :)

    real(dp) :: phi, r(2)
    integer :: i1, i2, ik, ierr
    PetscScalar :: phase

    do ik = 1, nk
      do i2 = -maxj, maxj
        do i1 = -maxi, maxi
          r = real(i1, 8)*brav(1:2, 1) + real(i2, 8)*brav(1:2, 2)
          phi = dot_product(r, kp(1:2, ik))
          phase = zexp(cmplx(0, phi, 8))
          if (same .eq. 0) then
            call MatAXPY(y(ik), phase, x(i1, i2), DIFFERENT_NONZERO_PATTERN, ierr)
          else
            call MatAXPY(y(ik), phase, x(i1, i2), SAME_NONZERO_PATTERN, ierr)
          end if
        end do
      end do
    end do

  end subroutine fourier_trans

! transfrom y(k)=sum_ij exp(k*R_ij)*x(i,j)
  subroutine fourier_trans_ik(x, y, kp, iik, brav, maxi, maxj, same)
#include <petsc/finclude/petsc.h>
    use petscmat
    use petsc_mod
    use kinds
!~       use blas95
!~       use globals
    implicit none

    Mat, intent(inout) :: x(-maxi:maxi, -maxj:maxj)
    Mat, intent(inout) :: y

    integer :: maxi, maxj, iik
    integer :: same
    real(dp) :: brav(3, 3)
    real(dp) :: kp(:, :)

    real(dp) :: phi, r(2)
    integer :: i1, i2, ik, ierr
    PetscScalar :: phase

    ik = iik
    do i2 = -maxj, maxj
      do i1 = -maxi, maxi
        r = real(i1, 8)*brav(1:2, 1) + real(i2, 8)*brav(1:2, 2)
        phi = dot_product(r, kp(1:2, ik))
        phase = zexp(cmplx(0, phi, 8))
        if (same .eq. 0) then
          call MatAXPY(y, phase, x(i1, i2), DIFFERENT_NONZERO_PATTERN, ierr)
        else
          call MatAXPY(y, phase, x(i1, i2), SAME_NONZERO_PATTERN, ierr)
        end if
      end do
    end do

  end subroutine fourier_trans_ik


  subroutine fourier3d_trans_ik(x_in, y, kp, iik, brav, maxi, maxj, maxk, same)
#include <petsc/finclude/petsc.h>
    use petscmat
    use petsc_mod
    use kinds
!~       use blas95
!~       use globals
    implicit none

!~     type(Mat_pointer_array), intent(in) :: x(-maxi:maxi, -maxj:maxj, -maxk:maxk)
!~ this does not preserve the correct bounds, -maxi:maxi -> 1:2*maxi+1.
!~ use a pointer for bound remapping.
    type(Mat_pointer_array), target, intent(in), contiguous :: x_in(:,:,:)
    Mat, intent(inout) :: y

    integer :: maxi, maxj, maxk, iik
    integer :: same
    real(dp) :: brav(3, 3)
    real(dp) :: kp(:, :)

    type(Mat_pointer_array), pointer :: x(:,:,:)
    real(dp) :: phi, r(3)
    integer :: i1, i2, i3, ik, ierr, n1, n2 ,n3
    PetscScalar :: phase

    n1 = size(x_in,1)/2
    n2 = size(x_in,2)/2
    n3 = size(x_in,3)/2
    x(-n1:n1,-n2:n2,-n3:n3)=>x_in(:,:,:)

    ik = iik
    do i3 = -maxk, maxk
      do i2 = -maxj, maxj
        do i1 = -maxi, maxi
          r = real(i1, 8)*brav(1:3, 1) + real(i2, 8)*brav(1:3, 2) + real(i3, 8)*brav(1:3, 3)
          phi = dot_product(r(1:3), kp(1:3, ik))
          phase = zexp(cmplx(0, phi, 8))
          if (same .eq. 0) then
            call MatAXPY(y, phase, x(i1, i2, i3)%p, DIFFERENT_NONZERO_PATTERN, ierr)
          else
            call MatAXPY(y, phase, x(i1, i2, i3)%p, SAME_NONZERO_PATTERN, ierr)
          end if
        end do
      end do
    end do

  end subroutine fourier3d_trans_ik

!~ ! transfrom y(i,j)=sum_k exp(k*R_ij)*x(k)
!~     subroutine fourier_back(x,y,kp,wkp,nk,brav,maxi,maxj,same)
!~ #include <slepc/finclude/slepceps.h>
!~       use slepceps
!~       use slepc_mod
!~       use petsc_mod
!~       use kinds
!~       use blas95
!~       implicit none

!~       Mat, allocatable, intent(inout) :: x(:)
!~       Mat, allocatable, intent(inout) :: y(:,:)
!~       real(dp) :: kp(:,:),brav(3,3)
!~       integer :: wkp(:),maxi,maxj,same

!~       real(dp) :: phi,r(2)
!~       integer :: i1,i2,ik,nk,ierr
!~       PetscScalar :: phase

!~       do i2=-maxj,maxj
!~         do i1=-maxi,maxi
!~           call MatZeroEntries(y(i1,i2),ierr)
!~         end do
!~       end do

!~       do ik=1,nk
!~         do i2=-maxj,maxj
!~           do i1=-maxi,maxi
!~             r=real(i1,8)*brav(1:2,1)+real(i2,8)*brav(1:2,2)
!~             phi=dot(r,kp(1:2,ik))
!~             phase=zexp(cmplx(0,-phi,8))
!~             if (same.eq.0) then
!~               call  MatAXPY(y(i1,i2),phase,x(ik),DIFFERENT_NONZERO_PATTERN,ierr)
!~               if (wkp(ik).eq.2) then
!~                 call MatTranspose(x(ik),MAT_INPLACE_MATRIX,x(ik),ierr)
!~                 phase=-phase
!~                 call  MatAXPY(y(i1,i2),phase,x(ik),DIFFERENT_NONZERO_PATTERN,ierr)
!~               end if
!~             else
!~               call  MatAXPY(y(i1,i2),phase,x(ik),SAME_NONZERO_PATTERN,ierr)
!~               if (wkp(ik).eq.2) then
!~                 call MatTranspose(x(ik),MAT_INPLACE_MATRIX,x(ik),ierr)
!~                 phase=-phase
!~                 call  MatAXPY(y(i1,i2),phase,x(ik),SAME_NONZERO_PATTERN,ierr)
!~               end if
!~             end if
!~           end do
!~         end do
!~       end do

!~     end subroutine fourier_back

! transfrom y(i,j)=sum_k exp(k*R_ij)*x(k)
! no y=0 on entry
! input matrices are changed to x->x^T if wkp.eq.2
  subroutine fourier_back_add(x, y, kp, wkp, nk, brav, maxi, maxj, ri, same, fac_in)
#include <petsc/finclude/petsc.h>
    use petscmat
    use petsc_mod
    use kinds
!~       use blas95
    use globals, only: pi
    implicit none

    PetscScalar, optional :: fac_in
    Mat, allocatable, target :: x(:)
    Mat :: y(-maxi:maxi, -maxj:maxj)
    Mat, pointer :: px
    real(dp) :: kp(:, :), brav(3, 3)
    integer :: wkp(:), maxi, maxj, same, ri

    real(dp) :: phi, r(2)
    integer :: i1, i2, ik, nk, ierr
    integer(8) :: i8
    PetscScalar :: phase, fac
    MatStructure :: matstr

! this will not work for restricted matrices, i.e., when we want to write out only what conquest actually reads.
! we need to write a custom MatAXPY. scheisse

    i8 = 0
    matstr = SAME_NONZERO_PATTERN
    if (same .ne. 0) matstr = DIFFERENT_NONZERO_PATTERN
    fac = 1d0
    if (present(fac_in)) fac = fac_in
! rotate Im->Re and -Re->Im, this messes up imaginary part of y!
! we could add Re/Im to petsc_aXpY
    if (ri .eq. 2) fac = -p_zione * fac

    do i2 = -maxj, maxj
      do i1 = -maxi, maxi
        r = real(i1, 8) * brav(1:2, 1) + real(i2, 8) * brav(1:2, 2)
        do ik = 1, nk
          px => x(ik)
          phi = -dot_product(r, kp(1:2, ik))
          phase = zexp(cmplx(0, phi, 8)) * fac
          call petsc_aXpY(y(i1, i2), x(ik), phase, PETSC_FALSE, .false.)

          if (wkp(ik) .eq. 2) then
!              call PetscBarrier(i8,ierr)
            call MPI_Barrier(MPI_COMM_WORLD, ierr)
            phi = dot_product(r, kp(1:2, ik))
            phase = zexp(cmplx(0, phi, 8)) * fac
            call petsc_aXpY(y(i1, i2), x(ik), phase, PETSC_FALSE, .true.)
          end if
        end do
      end do
    end do

  end subroutine fourier_back_add
  
!~ takes fortran arrays (1D) from each process and do the FT on each process independently  
  subroutine fourier_back_add_array(x, y, kp, wkp, nk, brav, maxi, maxj, ri)
#include <petsc/finclude/petsc.h>
    use petscmat
    use petsc_mod
    use kinds
!~       use blas95
    use globals, only: pi
    implicit none

    complex(dp) :: x(:, :)
    complex(dp) :: y(:, :, :, :) ! -maxi:maxi, -maxj:maxj)
    real(dp) :: kp(:, :), brav(3, 3)
    integer :: wkp(:), maxi, maxj, ri

    real(dp) :: phi, r(2)
    integer :: i1, i2, ik, nk, ierr, ii1, ii2
    PetscScalar :: phase, fac

    fac = 1d0
! rotate Im->Re and -Re->Im, this messes up imaginary part of y!
! we could add Re/Im to petsc_aXpY
    if (ri .eq. 2) fac = -p_zione

    do i2 = -maxj, maxj
      ii2 = i2 + maxj + 1
      do i1 = -maxi, maxi
        ii1 = i1 + maxi + 1
        r = real(i1, 8)*brav(1:2, 1) + real(i2, 8)*brav(1:2, 2)
        do ik = 1, nk          
          phi = -dot_product(r, kp(1:2, ik))
          phase = zexp(cmplx(0, phi, 8))*fac
          y(:, :, ii1, ii2) = y(:, :, ii1, ii2) + phase * x(:, :)
          if (wkp(ik) .eq. 2) then
            phi = dot_product(r, kp(1:2, ik))
            phase = zexp(cmplx(0, phi, 8))*fac
            y(:, :, ii1, ii2) = y(:, :, ii1, ii2) + phase * transpose(x(:, :))
          end if
        end do
      end do
    end do

  end subroutine fourier_back_add_array

  subroutine fourier3d_back_add(x, y_in, kp, wkp, nk, brav, maxi, maxj, maxk, ri, &
   same, lnewnz_in)
#include <petsc/finclude/petsc.h>
    use petscmat
    use petsc_mod
    use kinds
!~       use blas95
    use globals, only: pi
    implicit none

    Mat, allocatable :: x(:)
!~ this does not preserve the correct bounds, -maxi:maxi -> 1:2*maxi+1.
!~ use a pointer for bound remapping.
    Mat, target, contiguous :: y_in(:,:,:)
!~     type(Mat_pointer_array) :: y(-maxi:maxi, -maxj:maxj, -maxk:maxk)
    real(dp) :: kp(:, :), brav(3, 3)
    logical, optional :: lnewnz_in
    integer :: wkp(:), maxi, maxj, maxk, same, ri

    Mat, pointer :: y(:,:,:)
    real(dp) :: phi, r(3)
    integer :: i1, i2, i3, ik, nk, ierr, n1, n2, n3
    integer(8) :: i8
    PetscScalar :: phase, fac
    MatStructure :: matstr
    logical :: lnewnz

! this will not work for restricted matrices, i.e., when we want to write out only what conquest actually reads.
! we need to write a custom MatAXPY. scheisse

    lnewnz = PETSC_TRUE
    if (present(lnewnz_in)) lnewnz = lnewnz_in

    n1 = size(y_in,1)/2
    n2 = size(y_in,2)/2
    n3 = size(y_in,3)/2
    y(-n1:n1,-n2:n2,-n3:n3)=>y_in(:,:,:)


    i8 = 0
    matstr = SAME_NONZERO_PATTERN
    if (same .ne. 0) matstr = DIFFERENT_NONZERO_PATTERN
    fac = 1d0
! rotate Im->Re and -Re->Im, this messes up imaginary part of y!
! we could add Re/Im to petsc_aXpY
    if (ri .eq. 2) fac = -p_zione
    do i3 = -maxk, maxk
      do i2 = -maxj, maxj
        do i1 = -maxi, maxi
          r = real(i1, 8)*brav(1:3, 1) + real(i2, 8)*brav(1:3, 2) + real(i3, 8)*brav(1:3, 3)
          do ik = 1, nk
            phi = -dot_product(r(1:3), kp(1:3, ik))
            phase = zexp(cmplx(0, phi, 8))*fac
            call petsc_aXpY(y(i1, i2, i3), x(ik), phase, lnewnz, .false.)

            if (wkp(ik) .eq. 2) then
  !              call PetscBarrier(i8,ierr)
              call MPI_Barrier(MPI_COMM_WORLD, ierr)
              phi = dot_product(r, kp(1:3, ik))
              phase = zexp(cmplx(0, phi, 8))*fac
              call petsc_aXpY(y(i1, i2, i3), x(ik), phase, lnewnz, .true.)
            end if
          end do
        end do
      end do
    end do

  end subroutine fourier3d_back_add

end module ft_mod
