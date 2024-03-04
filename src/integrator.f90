module integrator_mod
#include <petsc/finclude/petscmat.h>
  use petscmat
  use kinds
  use misc
  use globals
!~   use blas95
!~   use lapack95
  use mathstuff
  use surf_gf_mod
  use petsc_mod, only: pstr_out, petsc_print_master
  implicit none

  integer :: intlr, n_fermipols
  real(dp) :: r_eq, e0, rpfix, ipfix, phi_low, eu, el
  complex(dp) :: y, x_low, y_low, dx_low, dy_low
  
  
  integer :: nint_nested, nsubdivs, error_reim_select
  real(dp) :: h_min
  real(dp), allocatable :: nodes(:), weight1(:), weight2(:), error_tracker(:)

  Mat :: p_tmp1_int1, p_tmp1_int2, p_tmp1_int3, p_tmp1_int4
  Mat :: p_tmp2_int1, p_tmp2_int2, p_tmp2_int3
    
  Mat, pointer :: p_zint_out
  Mat, allocatable :: p_dmat_error(:,:)

contains

  subroutine get_alpha(p_Dl, p_Dr, dc_weights)
#include <petsc/finclude/petscmat.h>
    use petsc_mod
    implicit none

    Mat :: p_Dl, p_Dr
    PetscScalar, allocatable :: dc_weights(:)

    PetscScalar :: a_l, a_r, p_val
    integer :: imu1, imu2, iat1, ierr, imu3
    Mat :: p_tmp_l, p_tmp_r    
    
    imu2 = 0
    imu3 = 0
    do iat1 = 1, nat_ecc
      a_l = 0d0
      a_r = 0d0
      do imu1 = 1, imu_ecc(iat1)
        imu2 = imu2 + 1
       call petsc_mat_getvalue(imu2 - 1, imu2 - 1, p_val, p_Dl, 1, PETSC_COMM_WORLD)
        a_l = a_l + p_val
       call petsc_mat_getvalue(imu2 - 1 , imu2 - 1 , p_val, p_Dr, 1, PETSC_COMM_WORLD)
        a_r = a_r + p_val
      end do
      p_val = real(a_l, 8) / real((a_l + a_r), 8)
      do imu1 = 1, imu_ecc(iat1)
        imu3 = imu3 + 1       
        dc_weights(imu3) = p_val
      end do
    end do
    
!~     do imu2 = 1, nmu_c
!~       write(0, fmt='(i8,2e24.12)') imu2, dc_tmp(imu2)
!~       do imu1 = 1, nmu_c
!~         dc_weights(imu1, imu2) = dc_tmp(imu1) * dc_tmp(imu2)
!~       end do
!~     end do
!~     write(0, fmt='(4e24.12)') &
!~    &  maxval(real(dc_weights)),minval(real(dc_weights)),&
!~    &  maxval(aimag(dc_weights)),minval(aimag(dc_weights))

  end subroutine get_alpha

  subroutine get_weight(wl, wr, alphal, alphar)

    implicit none

    complex(dp), allocatable :: wl(:, :), wr(:, :), alphal(:), alphar(:)
    complex(dp) :: znorm
    integer :: iat1, iat2

    do iat2 = 1, nat_ecc
      do iat1 = 1, nat_ecc
        znorm = zsqrt(alphal(iat1)*alphal(iat2)) + zsqrt(alphar(iat1)*alphar(iat2))
        wl(iat1, iat2) = zsqrt(alphal(iat1))*zsqrt(alphal(iat2))/znorm
        wr(iat1, iat2) = zsqrt(alphar(iat1))*zsqrt(alphar(iat2))/znorm
      end do
    end do

  end subroutine get_weight

  subroutine init_eq_int(mu1, mu2, lr)

    implicit none

    real(dp) :: mu1, mu2
    integer :: lr
    real(dp) :: eoff_l, eoff_u, mu_l, mu_u

    mu_l = min(mu1, mu2)
    mu_u = max(mu1, mu2)
    intlr = lr
    eoff_l = log(epsfermi/(1d0 - epsfermi))*kb*temperature_el + mu_l
    eoff_u = log(1d0/epsfermi - 1d0)*kb*temperature_el + mu_u
    r_eq = abs(elow - eoff_l)*0.5d0
    e0 = eoff_l - r_eq
    phi_low = pi - asin(delta_imag/r_eq)
    x_low = contour(phi_low)
    dx_low = real(eoff_u - x_low)
    n_fermipols = (delta_imag/(pi*kb*temperature_el) - 1)*0.5d0

    lnoneq_int = .false.

  end subroutine init_eq_int

  subroutine init_neq_int(lr)
#include <petsc/finclude/petscmat.h>
    use globals, only: l_ionode
    implicit none
    integer :: lr
    complex(dp) :: z1, z2
    
    integer :: ix, nx, iunit, ierr
    real(dp) :: dx
    complex(dp) :: xx

    intlr = lr
    el = log(epsfermi/(1d0 - epsfermi))*kb*temperature_el + min(mul, mur)
    eu = log(1d0/epsfermi - 1d0)*kb*temperature_el + max(mul, mur)
    z1 = el - mul
    z2 = el - mur
    
    write (pstr_out, fmt='(A, i3)') "intlr", intlr; call petsc_print_master()
    write (pstr_out, fmt='(A,4e24.12)') "el, eu", el, eu; call petsc_print_master()
    write (pstr_out, fmt='(A,2e24.12)') "low :", fermi(z1, temperature_el) - fermi(z2, temperature_el); call petsc_print_master()
    z1 = eu - mul
    z2 = eu - mur
    write (pstr_out, fmt='(A,2e24.12)') "up :", fermi(z1, temperature_el) - fermi(z2, temperature_el); call petsc_print_master()

    lnoneq_int = .true.
    
    if (l_ionode) then
      nx = 10000
      dx = (eu - el) / (real(nx, 8))
      xx = el
      open(newunit=iunit, file="fermitest.dat", action="write", status="replace")
      do ix = 0, nx
        xx = xx + dx
        write(iunit, fmt = '(6e24.12)') xx, fermi(xx - mul, temperature_el), fermi(xx - mur, temperature_el)
      end do
      close(iunit)
    end if

  end subroutine init_neq_int

  subroutine scale_Dne_dc(p_Dneq_dc, dc_weights)
#include <petsc/finclude/petscmat.h>
    use globals
      
    implicit none
    
    Mat :: p_Dneq_dc
    PetscScalar :: dc_weights(:)
    
    integer :: ierr, imu1, imu2,  nl1, nl2, i1, nlc1, nlc2, &
      nrow, ncol, nzrow, ii(1), nz
    integer, allocatable :: cols(:)
    PetscScalar, allocatable :: row_vals(:)
    PetscScalar :: p_val
    PetscReal :: weight
      
    call MatGetSize(p_Dneq_dc, nrow, ncol, ierr)
    call MatGetOwnershipRange(p_Dneq_dc, nl1, nl2, ierr)
    
    allocate(cols(ncol), row_vals(ncol), stat = ierr)
    
    do imu1 = 0, nmu_c - 1
      if ((imu1 .ge. nl1) .and. (imu1 .le. nl2 - 1)) then
        call MatGetRow(p_Dneq_dc, imu1, nzrow, cols, row_vals, ierr)
        do i1 = 1, nzrow
          imu2 = cols(i1)
!~           write(0, fmt='(2i8, 2e24.12)') imu1, imu2, dc_weights(imu2 + 1, imu1 + 1)
!~           row_vals(i1) = row_vals(i1) * dc_weights(imu2 + 1, (imu1 + 1))
          if (((real(dc_weights(imu1), 8).lt.0d0).or.(real(dc_weights(imu2), 8).lt.0d0)).and.&
         &  ((abs((aimag(dc_weights(imu1)))).ge.1d-10).or.&
         &   (abs((aimag(dc_weights(imu2)))).ge.1d-10))) then
            write(0, fmt='(A,2i8,4e24.12)') "warning: ", imu1, imu2, dc_weights(imu1), dc_weights(imu2)
          end if
          weight = dsqrt(real(dc_weights(imu1), 8) * real(dc_weights(imu2), 8))
          if (weight.lt.0d0) weight = 0d0
          row_vals(i1) = row_vals(i1) * weight
        end do
        ii(1) = imu1
        nz = nzrow
        call MatRestoreRow(p_Dneq_dc, imu1, nz, PETSC_NULL_INTEGER, PETSC_NULL_SCALAR, ierr)        
        call MatSetValues(p_Dneq_dc, 1, ii, nzrow, cols, row_vals, INSERT_VALUES, ierr)
      end if
      call MatAssemblyBegin(p_Dneq_dc, MAT_FINAL_ASSEMBLY, ierr)
      call MatAssemblyEnd(p_Dneq_dc, MAT_FINAL_ASSEMBLY, ierr)
    end do
    
  end subroutine scale_Dne_dc

  subroutine blockscale_mat(a, w, nat, imuecc)
    implicit none

    complex(dp), intent(inout) :: a(:, :)
    complex(dp), intent(in) :: w(:, :)
    integer, intent(in) :: nat, imuecc(:)
    integer :: imu1, imu2, iat1, iat2, i1, i2, j1, j2

    i2 = 1
    j2 = 1
    do iat2 = 1, nat
      j2 = i2 + imuecc(iat2) - 1
      i1 = 1
      j1 = 1
      do iat1 = 1, nat
        j1 = i1 + imuecc(iat1) - 1
        a(i1:j1, i2:j2) = w(iat1, iat2)*a(i1:j1, i2:j2)
        i1 = j1 + 1
      end do
      i2 = j2 + 1
    end do

  end subroutine blockscale_mat

  subroutine hpsort_aa(ra, rb, reverse)

    implicit none

! input,output:
    real(dp) :: ra(:)
    real(dp) :: rb(:, :)
! local :
    real(dp), allocatable :: rrb(:)
    real(dp) :: rra
    integer :: l, ir, n, j, i, m, minus
    logical, optional :: reverse

    minus = 1
    if (present(reverse)) then
      if (reverse) minus = -1
    end if

    n = size(ra, 1)
    m = size(rb, 1)
    allocate (rrb(m))
    l = n/2 + 1
    ir = n
    ra = ra*minus
    do
      if (l > 1) then
        l = l - 1
        rra = ra(l)
        rrb(:) = rb(:, l)
      else
        rra = ra(ir)
        rrb(:) = rb(:, ir)
        ra(ir) = ra(1)
        rb(:, ir) = rb(:, 1)
        ir = ir - 1
        if (ir .eq. 1) then
          ra(1) = rra
          rb(:, 1) = rrb(:)
          ra = ra*minus
          return
        end if
      end if
      i = l
      j = l + l
      do while (j .le. ir)
        if (j .lt. ir) then
          if (ra(j) .lt. ra(j + 1)) j = j + 1
        end if
        if (rra < ra(j)) then
          ra(i) = ra(j)
          rb(:, i) = rb(:, j)
          i = j
          j = j + j
        else
          j = ir + 1
        end if

      end do
      ra(i) = rra
      rb(:, i) = rrb(:)
    end do
    ra = ra*minus
  end subroutine hpsort_aa

  subroutine get_quadrule(integrator, xi, v1, v2, nint_order, nlow)
    use kinds
    use integrations_weights
    use error_handler

    implicit none

    integer :: integrator
    integer :: nint_order, nlow
    real(dp), allocatable :: xi(:), v1(:), v2(:)

    integer :: off, i, ierr
    real(dp), allocatable :: x2(:), w1(:), w2(:), vdummy(:, :)
    real(dp) :: dx

    if (integrator .eq. 1) then
      if (mod(nint_order, 2) .eq. 0) nint_order = nint_order + 1
      nlow = (nint_order - 1)/2
      if (l_output_progress) then
        write (pstr_out, fmt='(A,2i8)') " Adaptive Gauss Kronrod Integrator ", nlow, nint_order; call petsc_print_master()
      end if
      allocate (xi(nint_order), v1(nint_order), v2(nint_order), &
                x2(nlow + 1), w1(nlow + 1), w2(nlow + 1), vdummy(1, nint_order), stat=ierr)
      if (ierr .ne. 0) then
        write (errormsg, fmt='(A,i8)') "allocation error ", ierr
        call error()
      end if
      call kronrod(nlow, 1d-12, x2, w1, w2)
      off = 0

      do i = 1, nlow
        xi(i) = x2(i)
        v1(i) = w1(i)
        v2(i) = w2(i)
        xi(nint_order - i + 1) = -x2(i)
        v1(nint_order - i + 1) = w1(i)
        v2(nint_order - i + 1) = w2(i)
      end do
      xi(nlow + 1) = x2(nlow + 1)
      v1(nlow + 1) = w1(nlow + 1)
      v2(nlow + 1) = w2(nlow + 1)
            
    else if (integrator .eq. 2) then
      if ((mod(nint_order, 2)) .eq. 0) nint_order = nint_order + 1
      nlow = (nint_order - 1)/2 + 1
      if (l_output_progress) then
        write (pstr_out, fmt='(2i8,A)') nlow, nint_order, " Adaptive Clenshaw-Curtis integrator"; call petsc_print_master()
      end if
      allocate (xi(nint_order), v1(nint_order), v2(nint_order), &
                x2(nint_order), w1(nint_order), w2(nint_order), vdummy(1, nint_order), stat=ierr)
      x2 = 0d0
      w1 = 0d0
      xi = 0d0
      v1 = 0d0
      v2 = 0d0
      call clenshaw_curtis_compute(nlow, x2(1:nlow), w2(1:nlow))
      call clenshaw_curtis_compute(nint_order, x2, w1)
      do i = 1, nlow
        v1(i) = w1(2*(i - 1) + 1)
        xi(i) = x2(2*(i - 1) + 1)
        v2(i) = w2(i)
        if (i .eq. nlow) exit
        v1(nlow + i) = w1(2*i)
        xi(nlow + i) = x2(2*i)
      end do

    else if (integrator .eq. 3) then

      if (mod(nint_order, 2) .eq. 0) nint_order = nint_order + 1
      nlow = (nint_order)/2 + 1
      if (l_output_progress) then
        write (pstr_out, fmt='(2i8,A)') nlow, nint_order, " Adaptive Newton-Cotes integrator"; call petsc_print_master()
      end if
      allocate (xi(nint_order), v1(nint_order), v2(nint_order), &
                x2(nint_order), w1(nint_order), w2(nint_order), vdummy(1, nint_order), stat=ierr)
      if (ierr .ne. 0) then
        write (errormsg, fmt='(A,i8)') "allocation error ", ierr
        call error()
      end if
      x2 = 0d0
      w2 = 0d0
      call line_ncc_rule(nlow, -1d0, 1d0, x2(1:nlow), w2(1:nlow))
      call line_ncc_rule(nint_order, -1d0, 1d0, x2, w1)
      do i = 1, nlow
        v1(i) = w1(2*(i - 1) + 1)
        xi(i) = x2(2*(i - 1) + 1)
        v2(i) = w2(i)
        if (i .eq. nlow) exit
        v1(nlow + i) = w1(2*i)
        xi(nlow + i) = x2(2*i)
      end do

    else if (integrator .eq. 4) then

      if (mod(nint_order, 2) .eq. 0) nint_order = nint_order + 1
      nlow = (nint_order)/2 + 1
      if (l_output_progress) then
        write (pstr_out, fmt='(2i8,A)') nlow, nint_order, " Adaptive trapezoidal integrator"; call petsc_print_master()
      end if
      allocate (xi(nint_order), v1(nint_order), v2(nint_order), &
                x2(nint_order), w1(nint_order), w2(nint_order), vdummy(1, nint_order), stat=ierr)
      if (ierr .ne. 0) then
        write (errormsg, fmt='(A,i8)') "allocation error ", ierr
        call error()
      end if
      x2 = 0d0
      w2 = 0d0
      v1 = 0d0
      v2 = 0d0
      v1(1) = 1d0
      v1(nint_order) = 1d0
      v2(1) =1d0
      v2(nint_order) = 1d0
      xi(1) = -1d0
      xi(nint_order) =  1d0
      dx = 2d0 / real((nint_order - 1), 8)
      do i = 2, nint_order - 1
        v1(i) = 0.5d0
        xi(i) = xi(1) + dx * (i - 1)
        if (mod(i + 1, 2).eq.0) v2(i) = 0.5d0
      end do
      v1 = v1 * (xi(2) - xi(1))
      v2 = v2 * (xi(3) - xi(1))
    end if

    deallocate (x2, w1, w2)
    allocate (x2(nint_order))
    x2 = xi
    vdummy(1, :) = v1(:)
    call hpsort_aa(xi, vdummy)
    v1(:) = vdummy(1, :)
    vdummy(1, :) = v2(:)
    call hpsort_aa(x2, vdummy)
    v2(:) = vdummy(1, :)

!~    if (inode.eq.0) then
!~      do i=1,nint_order
!~        write(0,fmt='(i8,3e24.12)') i,xi(i),v1(i),v2(i)
!~      end do
!~    end if

    if (l_output_progress) then
      write (pstr_out, fmt='(A)') repeat("-", 4); call petsc_print_master()
    end if

  end subroutine get_quadrule

  function fermi(e, t)

    use kinds

    implicit none

    complex(dp) :: ii = dcmplx(0d0, 1d0)
    complex(dp) :: fermi, e, x, ze
    real(dp) :: t, kt, infinity

    kt = t*kb
    x = e/kt

    ze = zexp(x)
    infinity = HUGE(0_dp)
    if (real(ze, 8) .gt. infinity) then
      fermi = 1d-32
    else
      fermi = 1d0/(1d0 + zexp(x))
    end if

    if (fermi.ne.fermi) fermi=0d0

  end function fermi

  function bose(e, t)

    use kinds

    implicit none

    complex(dp) :: ii = dcmplx(0d0, 1d0)
    complex(dp) :: bose, e, x, ze, b
    real(dp) :: t, kt, infinity

    kt = t*kb
    x = e/kt

    ze = zexp(x)
    b = 1d0/(zexp(x) - 1d0)
    if (isnan(real(b))) b = 0d0
    if (isnan(aimag(b))) b = 0d0
    bose = b

  end function bose

  function contour(x)

    use kinds

    implicit none

    real(dp) :: x
    complex(dp) :: contour

    contour = cmplx(-cos(x), sin(x), 8)
    contour = r_eq*contour + e0

  end function contour

  function dcontour(x)

    use kinds

    implicit none

    real(dp) :: x
    complex(dp) :: dcontour

    dcontour = r_eq*cmplx(sin(x), cos(x), 8)

  end function dcontour

  function contour_x(x)

    use kinds

    implicit none

    real(dp) :: x
    complex(dp) :: contour_x

    contour_x = x_low + dx_low * x

  end function contour_x

  function dcontour_x(x)

    use kinds

    implicit none

    real(dp) :: x
    complex(dp) :: dcontour_x

    dcontour_x = dx_low

  end function dcontour_x

  function contour_y(x)

    use kinds

    implicit none

    real(dp) :: x
    complex(dp) :: contour_y

    contour_y = y_low + dy_low*x

  end function contour_y

  function dcontour_y(x)

    use kinds

    implicit none

    real(dp) :: x
    complex(dp) :: dcontour_y

    dcontour_y = dy_low

  end function dcontour_y

  subroutine get_fermipole_contribution(p_dmat, mu)
#include <petsc/finclude/petsc.h>
    use petscmat
    use petsc_mod
    use petsc_wrapper
    implicit none

    Mat :: p_dmat
    real(dp) :: mu

    complex(dp) :: z, fac
    integer :: ifp, ierr

    call MatZeroEntries(p_dmat, ierr)

    do ifp = 0, n_fermipols
      z = cmplx(mu, (2d0*ifp + 1)*pi*kb*temperature_el)
      write (pstr_out, fmt='(i8,2e24.12)') ifp, z; call petsc_print_master()
      ierr = eval_gr(z)
      fac = zione*2d0*pi*kb*temperature_el
      call MatAXPY(p_dmat, fac, p_tmpcc1, SAME_NONZERO_PATTERN, ierr)
!~         dmat=dmat+zione*2d0*pi*kb*temp*tmpcc
    end do

  end subroutine get_fermipole_contribution

  subroutine evalqr(f, a, b, x, w1, w2, p_zint, p_zint2, nhigh, nlow, errout, errout2, d1, d2)
#include <petsc/finclude/petsc.h>
    use petscmat
    use petsc_mod
    use petsc_wrapper
    use kinds
    use mathstuff
    use globals, only: ldouble_contour, lnoneq_int

    implicit none

    real(dp) ::a, b, errout, errout2
    real(dp), optional :: d1, d2
    real(dp), allocatable :: x(:), w1(:), w2(:)
    Mat :: p_zint, p_zint2
    integer :: nhigh, nlow, ierr
    integer, external :: f

    Mat :: p_zintg, p_zintg2
    integer :: i, n, j, i_low, i_high, i_max, k1, k2
    real(dp) :: xx, norm, m1, m2, timing_local
    integer(8) :: counti, count_rate, countf

    PetscScalar, allocatable :: ws1(:), ws2(:)

    allocate (ws1(nhigh), ws2(nhigh))

    ws1 = w1*(b - a)*0.5d0
    ws2 = w2*(b - a)*0.5d0

    call MatZeroEntries(p_zint, ierr)

    call MatDuplicate(p_zint, MAT_SHARE_NONZERO_PATTERN, p_zintg, ierr)

    if ((ldouble_contour) .and. (lnoneq_int)) then
      call MatDuplicate(p_zint2, MAT_SHARE_NONZERO_PATTERN, p_zintg2, ierr)
      call MatZeroEntries(p_zint2, ierr)
    end if
    timing_local = 0d0
    if (l_output_progress) then
      write (pstr_out, fmt='(i3,e16.6)') 0, timing_local; call PetscPrintf(PETSC_COMM_WORLD, trim(pstr_out), ierr)
    end if
    do i = 1, nhigh
      call system_clock(counti, count_rate)
      xx = x(i)*(b - a)*0.5d0 + (a + b)*0.5d0

      ierr = f(xx)

      call MatAXPY(p_zint, ws1(i), p_tmpcc1, SAME_NONZERO_PATTERN, ierr)

      if (ws2(i) .ne. 0d0) then
        call MatAXPY(p_zintg, ws2(i), p_tmpcc1, SAME_NONZERO_PATTERN, ierr)
      end if

      if ((ldouble_contour) .and. (lnoneq_int)) then
        call MatAXPY(p_zint2, ws1(i), p_tmpcc2, SAME_NONZERO_PATTERN, ierr)
        if (ws2(i) .ne. 0d0) call MatAXPY(p_zintg2, ws2(i), p_tmpcc2, SAME_NONZERO_PATTERN, ierr)
      end if

      int_counter = int_counter + 1
      call system_clock(countf)
      timing_local = timing_local + real(countf - counti, 8)/real(count_rate, 8)
      if (l_output_progress) then
        write (pstr_out, fmt='(A,i3,e16.6)') repeat(ACHAR(8), 19), i, timing_local; call PetscPrintf(PETSC_COMM_WORLD, trim(pstr_out), ierr)
      end if
    end do
    call PetscPrintf(PETSC_COMM_WORLD, " int eps", ierr)

    call MatAXPY(p_zintg, p_minus1, p_zint, SAME_NONZERO_PATTERN, ierr)
    call MatNorm(p_zintg, NORM_FROBENIUS, errout, ierr)
    if ((ldouble_contour) .and. (lnoneq_int)) then
      call MatAXPY(p_zintg2, p_minus1, p_zint2, SAME_NONZERO_PATTERN, ierr)
      call MatNorm(p_zintg2, NORM_FROBENIUS, errout2, ierr)
!~       write(0,*) errout, errout2
!~       errout = (errout + errout2)*0.5d0
      errout = max(errout, errout2)
      errout2 = errout
    else
      errout2 = errout
    end if

    call MatDestroy(p_zintg, ierr)

    if ((ldouble_contour) .and. (lnoneq_int)) then
      call MatDestroy(p_zintg2, ierr)
    end if

    if (.not. present(d2)) return

    d1 = a + (b - a)/3d0
    d2 = a + (b - a)*2d0/3d0

  end subroutine evalqr

  function get_gr_cc_neq(x)
#include <petsc/finclude/petscmat.h>
    use kinds
    use petsc_mod
    use petsc_wrapper
    use globals, only: eta_elec, iik, vb, nmu_l, nmu_r, nmu_c, wkp_r, eta_cc

    implicit none

    real(dp) :: x

    integer:: get_gr_cc_neq

    Mat :: p_gammal_block, p_gammar_block, p_x, b, p_g_block, p_rhs
    Mat, pointer :: p_gr_inv

    complex(dp) :: z

    integer :: ierr
    real(dp) :: vbb
    logical :: l_solvemode_5
    
    PetscScalar :: p_scal
    PetscReal :: norm

    get_gr_cc_neq = 0

    z = x

    p_gr_inv => p_invGr ! see eval_gr
    
    call init_gr(z, p_gr_inv)

    l_solvemode_5 = solver_mode .eq. 5
    
    if (intlr .eq. 2) then  ! GrGammaLGa*(f_L-f_R)
      if (l_solvemode_5) solver_mode = 4
      call petsc_get_a_with_b_c(p_g_block, p_gr_inv, p_gammal, mattype_dense)    
      if (solver_mode .eq. 2) then
        call petsc_get_a_with_b_c(p_rhs, p_gr_inv, p_gammal, mattype_dense)      
        call petsc_one_mat(p_rhs, 0, nmu_l - 1) ! prepare RHS nmu_c,nmu_l (1 0 0 .. 0,0 1 0 .. 0, ..) C style 1=0
        call petsc_call_solver(p_gr_inv, p_rhs, p_g_block, matsolvertype_cc, solver_mode)
      else if ((solver_mode .eq. 1).or.(solver_mode .eq. 3).or.(solver_mode .eq. 4)) then
        call petsc_get_a_with_b_c(p_rhs, p_gr_inv, p_gammal, mattype_sparse)
        call petsc_one_mat(p_rhs, 0, nmu_l - 1,2) ! prepare RHS nmu_c,nmu_l (1 0 0 .. 0,0 1 0 .. 0, ..) C style 1=0
        call petsc_call_solver(p_gr_inv, p_rhs, p_g_block, matsolvertype_cc, solver_mode)
      end if
      if (l_solvemode_5) solver_mode = 5


      call MatDestroy(p_rhs, ierr)

      call MatMatMult(p_g_block, p_gammal, MAT_INITIAL_MATRIX, 1d0, p_x, ierr) ! X_block=(G11*Gamma_L,G12*Gamma_L,G13*Gamma_L)
      call petsc_matmatconj_restricted(p_x, p_g_block, p_tmpcc1)      

!~ ! debug start --
!~         call MatNorm(p_gammal,NORM_FROBENIUS,norm,ierr)
!~         if (inode.eq.0)  write(0,fmt='(2e16.8,A,e16.8)') z, " Gammal ",norm
!~         call MatNorm(p_tmpcc2,NORM_FROBENIUS,norm,ierr)
!~         if (inode.eq.0)  write(0,fmt='(2e16.8,A,e16.8)') z, " G*Gammal*G' ",norm
!~ ! debug end --

      p_scal = 0.5d0*(fermi(z - mul, temperature_el) - fermi(z - mur, temperature_el)) ! add weight functions, fermi window and absorb factor 0.5d0
      call MatScale(p_tmpcc1, p_scal, ierr)
      
      call MatDestroy(p_x, ierr)
      call MatDestroy(p_rhs, ierr)
      call MatDestroy(p_g_block, ierr)             

    end if

    if (intlr .eq. 1) then  ! GrGammaRGa*(f_R-f_L)
    
      if (l_solvemode_5) solver_mode = 4
      call petsc_get_a_with_b_c(p_g_block, p_gr_inv, p_gammar, mattype_dense)    
      if (solver_mode .eq. 2) then
        call petsc_get_a_with_b_c(p_rhs, p_gr_inv, p_gammar, mattype_dense)      
        call petsc_one_mat(p_rhs, nmu_c - nmu_r + 1 - 1, nmu_c - 1) ! prepare RHS nmu_c,nmu_r ( .., 0 0 .. 1 0, 0 0 .. 0 1 )   C style 1=0
        call petsc_call_solver(p_gr_inv, p_rhs, p_g_block, matsolvertype_cc, solver_mode)
      else if ((solver_mode .eq. 1).or.(solver_mode .eq. 3).or.(solver_mode .eq. 4)) then
        call petsc_get_a_with_b_c(p_rhs, p_gr_inv, p_gammar, mattype_sparse)
        call petsc_one_mat(p_rhs, nmu_c - nmu_r + 1 - 1, nmu_c - 1, 2) ! prepare RHS nmu_c,nmu_r ( .., 0 0 .. 1 0, 0 0 .. 0 1 )   C style 1=0
        call petsc_call_solver(p_gr_inv, p_rhs, p_g_block, matsolvertype_cc, solver_mode)
      end if
      if (l_solvemode_5) solver_mode = 5
  
      call MatDestroy(p_rhs, ierr)
  
      call MatMatMult(p_g_block, p_gammar, MAT_INITIAL_MATRIX, 1d0, p_x, ierr) ! X_block=(G13*Gamma_R,G23*Gamma_R,G33*Gamma_R)    
      call petsc_matmatconj_restricted(p_x, p_g_block, p_tmpcc1)
!~      ! debug start --
!~        call MatNorm(p_gammar,NORM_FROBENIUS,norm,ierr)
!~        if (inode.eq.0)  write(0,fmt='(2e16.8,A,e16.8)') z," Gammar ",norm
!~        call MatNorm(p_tmpcc1,NORM_FROBENIUS,norm,ierr)
!~        if (inode.eq.0)  write(0,fmt='(2e16.8,A,e16.8)') z," G*Gammar*G' ",norm
!~  ! debug end --
  
      p_scal = 0.5d0*(fermi(z - mur, temperature_el) - fermi(z - mul, temperature_el)) ! add weight functions, fermi window and absorb factor 0.5d0
      call MatScale(p_tmpcc1, p_scal, ierr)
  
      call MatDestroy(p_x, ierr)
      call MatDestroy(p_rhs, ierr)
      call MatDestroy(p_g_block, ierr)    

    end if

!~     call MatDestroy(p_gr_inv, ierr)

  end function get_gr_cc_neq

  function get_gr_cc_eq(x)

    use kinds
    use petsc_mod
    use petsc_wrapper
    use globals, only: eta_elec, iik, vb, nmu_l, nmu_r, nmu_c, wkp_r, eta_cc, ngroups
    use error_handler
    use pexsi_wrapper

    implicit none

    real(dp) :: x

    integer:: get_gr_cc_eq

    Mat :: p_mat_one
    Mat, pointer :: p_gr_inv

    complex(dp) :: z, dc
    integer :: ierr
    PetscScalar :: p_scal
    PetscReal :: norm

    get_gr_cc_eq = 0

    if (contour_select .eq. 1) then
      z = contour(x)
      dc = dcontour(x)
    else if (contour_select .eq. 2) then
      z = contour_x(x)
      dc = dcontour_x(x)
    else
      write (errormsg, fmt='(A,i8)') "contour_select invalid", contour_select
      call error()
    end if

    p_gr_inv => p_invGr ! see eval_gr

    call init_gr(z, p_gr_inv)

    if (solver_mode .eq. 2) then
      call petsc_get_densemat(p_gr_inv, p_mat_one, mattype_dense)
      call petsc_one_mat(p_mat_one, 0, nmu_c - 1)
      call petsc_call_solver(p_gr_inv, p_mat_one, p_tmpcc1, matsolvertype_cc, solver_mode)
      call MatDestroy(p_mat_one, ierr)
    else if ((solver_mode .eq. 1).or.(solver_mode .eq. 3).or.(solver_mode .eq. 4)) then
      call petsc_get_a_with_b_c(p_mat_one, p_gr_inv, p_gr_inv, mattype_sparse)
      call petsc_one_mat(p_mat_one, 0, nmu_c - 1, 2)
      call petsc_call_solver(p_gr_inv, p_mat_one, p_tmpcc1, matsolvertype_cc, solver_mode)
      call MatDestroy(p_mat_one, ierr)
    else if (solver_mode .eq. 5) then
      call inv_sel(p_gr_inv, p_tmpcc1)
    end if

    !~ -1 comes from the -1/pi factor, 1/pi will be added in the add.
    if (intlr .eq. 1) then ! mul <= mur
      p_scal = -dc*fermi(z - mul, temperature_el)
    else if (intlr .eq. 2) then ! mur < mul
      p_scal = -dc*fermi(z - mur, temperature_el)
    else if (intlr .eq. 3) then
      p_scal = -dc*(fermi(z - mul, temperature_el) + fermi(z - mur, temperature_el))*0.5d0
    else if (intlr .eq. 4) then
      p_scal = -dc*(fermi(z - mul, temperature_el) - fermi(z - mur, temperature_el))
    else
      write (errormsg, fmt='(A,i8)') "only left=1 and right=2 electrodes", intlr
      call error()
    end if

    call MatScale(p_tmpcc1, p_scal, ierr)

!~     call MatDestroy(p_gr_inv, ierr)

  end function get_gr_cc_eq

!~     subroutine init_gr(z,p_tmp1)
!~ initialises G^r ,i.e. p_tmp1=E*S-H-Sigma_L-Sigma_R
  subroutine init_gr(z, p_tmp1)
#include <petsc/finclude/petscmat.h>
    use petsc_mod
    use petsc_wrapper
    use globals, only: eta_elec, iik
    use kinds
    use error_handler

    implicit none

    complex(dp) :: z
    Mat :: p_tmp1

    integer :: ierr
    real(dp) :: vbb
    complex(dp) :: zz
    PetscReal :: norm

    ierr = eval_sigma(z)
    call MatZeroEntries(p_tmp1, ierr)
    call MatAXPY(p_tmp1, p_minus1, p_h00k_cc(iik), DIFFERENT_NONZERO_PATTERN, ierr)
    zz = (z + eta_cc*zione)
    call MatAXPY(p_tmp1, zz, p_s00k_cc(iik), DIFFERENT_NONZERO_PATTERN, ierr)

    if (l_use_sigma_l) call petsc_add_sub_B_to_A(p_sigmalr, p_tmp1, 0, 0, &
   &  p_minus1, ADD_VALUES, PETSC_FALSE)
    if (l_use_sigma_r) call petsc_add_sub_B_to_A(p_sigmarr, p_tmp1, &
   &  nmu_c - nmu_r, nmu_c - nmu_r, p_minus1, ADD_VALUES, PETSC_FALSE)

!bias shift
     if (((vl.ne.0d0).or.(vr.ne.0d0)).and.l_diag_fixH) then
       zz = -vr
!~        call petsc_add_sub_B_to_A(p_s10k_r(iik), p_tmp1, nmu_c - nmu_r, nmu_c - nmu_r*2, zz, ADD_VALUES, PETSC_FALSE)       
!~        call petsc_add_sub_B_to_A(p_s01k_r(iik), p_tmp1, nmu_c - nmu_r*2, nmu_c - nmu_r, zz, ADD_VALUES, PETSC_FALSE)      
       call petsc_add_sub_B_to_A(p_s00k_r(iik), p_tmp1, nmu_c - nmu_r, nmu_c - nmu_r, zz, ADD_VALUES, PETSC_FALSE)            
!~        call petsc_add_sub_B_to_A(p_s00k_r(iik), p_tmp1, nmu_c - nmu_r*2, nmu_c - nmu_r*2, zz, ADD_VALUES, PETSC_FALSE)
       
     
       zz = -vl ! this shall be change to genral bias voltage in left and right if vl and vr       
       call petsc_add_sub_B_to_A(p_s00k_l(iik), p_tmp1, 0, 0, zz, ADD_VALUES, PETSC_FALSE)
!~        call petsc_add_sub_B_to_A(p_s00k_l(iik), p_tmp1, nmu_l, nmu_l, zz, ADD_VALUES, PETSC_FALSE)
!~        call petsc_add_sub_B_to_A(p_s10k_l(iik), p_tmp1, nmu_l, 0, zz, ADD_VALUES, PETSC_FALSE)
!~        call petsc_add_sub_B_to_A(p_s01k_l(iik), p_tmp1, 0, nmu_l, zz, ADD_VALUES, PETSC_FALSE)
     end if

  end subroutine init_gr

  function eval_gr(z)
#include <petsc/finclude/petscmat.h>
    use petsc_mod
    use petsc_wrapper
    use globals, only: nsim_inv
    use kinds
    use pexsi_wrapper

    implicit none

    complex(dp) :: z

    integer :: eval_gr

    Mat :: p_tmp1, p_tmp2

    Mat, pointer :: p_gr_inv
    Mat :: p_mat_one
    integer :: ierr

    eval_gr = 0

!~ Bug in PETSC 3.16 causes preallocator matrix to be messed up after first call.
!~ So we initialize p_invGr in init_sys.f90 and reuse it in each call to init_gr,
!~ where it's entries are set to 0 at the beginning. To keep the structure 
!~ of the code and in case we want to go back use pointer to p_invG and leave most
!~ of the code unchanged, save the call to preallocated in init_gr.
    p_gr_inv => p_invGr
    
    call init_gr(z, p_gr_inv)

    if (solver_mode .eq. 2) then
      call petsc_get_densemat(p_gr_inv, p_mat_one, mattype_dense)
      call petsc_one_mat(p_mat_one, 0, nmu_c - 1)
      call petsc_call_solver(p_gr_inv, p_mat_one, p_tmpcc1, matsolvertype_cc, solver_mode)
      call MatDestroy(p_mat_one, ierr)
    else if ((solver_mode .eq. 1).or.(solver_mode .eq. 3).or.(solver_mode .eq. 4)) then
      call petsc_get_a_with_b_c(p_mat_one, p_gr_inv, p_gr_inv, mattype_sparse)
      call petsc_one_mat(p_mat_one, 0, nmu_c - 1, 2)
      call petsc_call_solver(p_gr_inv, p_mat_one, p_tmpcc1, matsolvertype_cc, solver_mode)
      call MatDestroy(p_mat_one, ierr)
    else if (solver_mode .eq. 5) then
      call inv_sel(p_gr_inv, p_tmpcc1)
    end if

!~     call MatDestroy(p_gr_inv, ierr)

  end function eval_gr

  function eval_sigma(z)
#include <petsc/finclude/petscmat.h>
    use petsc_mod
    use petsc_wrapper
    use kinds
!~       use blas95
    use globals, only: eta_elec, iik

    implicit none

    complex(dp) :: z

    integer:: eval_sigma

    complex(dp) :: zz
    integer :: ierr
    PetscReal :: norml, normr, normgrr, normgll
    PetscScalar :: pl11,plnn,pr11,prnn

    Mat :: p_tmp1, p_tmp2, p_tmp3

    eval_sigma = 0

    
! sigma_l

!~       call petsc_get_densemat(p_h00k_l(iik),p_tmp1,mattype_surf)
!~       call petsc_get_densemat(p_h00k_l(iik),p_tmp2,mattype_surf)
!~       call petsc_get_densemat(p_h00k_l(iik),p_tmp3,mattype_surf)

    if (l_use_sigma_l) then
      call surf_gf3_petsc(z - vl + cmplx(0d0, eta_elec), iik, "l")
!~       call MatNorm(p_gllr, NORM_FROBENIUS, normgll, ierr)
      
      ! pull out -1 factor to avoid calling MatScale twice !!! (z*S_CL-H_CL)*g_LL*(z*S_LC-H_LC)
      zz = -(z - vl + cmplx(0d0, eta_elec))
  
      call MatDuplicate(p_s10k_l(iik), MAT_COPY_VALUES, p_tmp1, ierr)
      call MatScale(p_tmp1, zz, ierr)
      call MatAXPY(p_tmp1, p_one, p_h10k_l(iik), DIFFERENT_NONZERO_PATTERN, ierr)
  
      call MatDuplicate(p_s01k_l(iik), MAT_COPY_VALUES, p_tmp2, ierr)
      call MatScale(p_tmp2, zz, ierr)
      call MatAXPY(p_tmp2, p_one, p_h01k_l(iik), DIFFERENT_NONZERO_PATTERN, ierr)
  
      call MatConvert(p_tmp1, mattype_surf, MAT_INPLACE_MATRIX, p_tmp1, ierr)
      call MatConvert(p_tmp2, mattype_surf, MAT_INPLACE_MATRIX, p_tmp2, ierr)
      call MatMatMult(p_tmp1, p_gllr, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, p_tmp3, ierr)
      call MatDestroy(p_tmp1, ierr)
      call MatMatMult(p_tmp3, p_tmp2, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, p_tmp1, ierr) ! there is no ELEMENTAL*ELEMENTAL=DENSE so we need MatConvert
      call MatConvert(p_tmp1, mattype_dense, MAT_REUSE_MATRIX, p_sigmalr, ierr)

! gammaL=i*(sigmaLr-sigmaLa)
!~     call MatNorm(p_sigmalr, NORM_FROBENIUS, norml, ierr)
      call MatTransposeSetPrecursor(p_sigmalr, p_gammal, ierr)
      call MatTranspose(p_sigmalr, MAT_REUSE_MATRIX, p_gammal, ierr)
      call MatConjugate(p_gammal, ierr)
      call MatAYPX(p_gammal, p_minus1, p_sigmalr, SAME_NONZERO_PATTERN, ierr) !p_gammal=(sigmaLr-sigmaLa)
      call MatScale(p_gammal, p_zione, ierr)
  
      call MatDestroy(p_tmp1, ierr)
      call MatDestroy(p_tmp2, ierr)
      call MatDestroy(p_tmp3, ierr)
    end if

! sigma_r

!~       call petsc_get_densemat(p_h00k_r(iik),p_tmp1,mattype_surf)
!~       call petsc_get_densemat(p_h00k_r(iik),p_tmp2,mattype_surf)
!~       call petsc_get_densemat(p_h00k_r(iik),p_tmp3,mattype_surf)
    if (l_use_sigma_r) then
      call surf_gf3_petsc(z - vr + cmplx(0d0, eta_elec), iik, "r")
!~       call MatNorm(p_grrr, NORM_FROBENIUS, normgrr, ierr)
      ! pull out -1 factor to avoid calling MatScale twice !!! (z*S_CR-H_CR)*g_RR*(z*S_RC-H_RC)
      zz = -(z - vr + cmplx(0d0, eta_elec))
  
      call MatDuplicate(p_s01k_r(iik), MAT_COPY_VALUES, p_tmp1, ierr)
      call MatScale(p_tmp1, zz, ierr)
      call MatAXPY(p_tmp1, p_one, p_h01k_r(iik), DIFFERENT_NONZERO_PATTERN, ierr)
  
      call MatDuplicate(p_s10k_r(iik), MAT_COPY_VALUES, p_tmp2, ierr)
      call MatScale(p_tmp2, zz, ierr)
      call MatAXPY(p_tmp2, p_one, p_h10k_r(iik), DIFFERENT_NONZERO_PATTERN, ierr)
  
      call MatConvert(p_tmp1, mattype_surf, MAT_INPLACE_MATRIX, p_tmp1, ierr)
      call MatConvert(p_tmp2, mattype_surf, MAT_INPLACE_MATRIX, p_tmp2, ierr)
      call MatMatMult(p_tmp1, p_grrr, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, p_tmp3, ierr)
      call MatDestroy(p_tmp1, ierr)
      call MatMatMult(p_tmp3, p_tmp2, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, p_tmp1, ierr) ! there is no ELEMENTAL*ELEMENTAL=DENSE so we need MatConvert
      call MatConvert(p_tmp1, mattype_dense, MAT_REUSE_MATRIX, p_sigmarr, ierr)

! gammaR=i*(sigmaRr-sigmaRa)
!~     call MatNorm(p_sigmarr, NORM_FROBENIUS, normr, ierr)
      call MatTransposeSetPrecursor(p_sigmarr, p_gammar, ierr)
      call MatTranspose(p_sigmarr, MAT_REUSE_MATRIX, p_gammar, ierr)
      call MatConjugate(p_gammar, ierr)
      call MatAYPX(p_gammar, p_minus1, p_sigmarr, SAME_NONZERO_PATTERN, ierr) !p_gammar=(sigmaRr-sigmaRa)
      call MatScale(p_gammar, p_zione, ierr)
  
      call MatDestroy(p_tmp1, ierr)
      call MatDestroy(p_tmp2, ierr)
      call MatDestroy(p_tmp3, ierr)
!~     call petsc_mat_getvalue(0, 0, pl11, p_sigmalr, 1, PETSC_COMM_WORLD, lcast_in=.true.)
!~     call petsc_mat_getvalue(0, 0, pr11, p_sigmarr, 1, PETSC_COMM_WORLD, lcast_in=.true.)
!~     call petsc_mat_getvalue(nmu_l - 1 - 9 + 1, nmu_l - 1 - 9 + 1, plnn, p_sigmalr, 1, PETSC_COMM_WORLD, lcast_in=.true.)
!~     call petsc_mat_getvalue(nmu_r - 1 - 9 + 1, nmu_r - 1 - 9 + 1, prnn, p_sigmarr, 1, PETSC_COMM_WORLD, lcast_in=.true.)
!~     if (inode.eq.0) then
!~       write(0,fmt='(A,14e15.6)') "check sigma ",zz, normgll, normgrr, norml, normr, pl11, pr11, plnn, prnn
!~     end if
    end if
  end function eval_sigma
  

  
  subroutine integrate(f, a, b, p_zint1, eps, max_error)
#include <petsc/finclude/petsc.h>
    use petscmat
    use petsc_mod
    use petsc_wrapper
    use globals, only : maxsub
    use kinds

    implicit none  
    
    integer, external :: f
    real(dp) :: a, b, eps, max_error
    Mat, target :: p_zint1
    
    integer :: ierr
    
    p_zint_out => p_zint1
    
    allocate(error_tracker(maxsub))
    error_tracker = 0d0
    nsubdivs = 0
    
    call get_quadrule(integrator, nodes, weight1, weight2, nint_order, nint_nested)
    
    call MatDuplicate(p_zint_out, MAT_SHARE_NONZERO_PATTERN, p_tmp1_int1, ierr)
!~     call MatDuplicate(p_zint_out, MAT_SHARE_NONZERO_PATTERN, p_tmp1_int2, ierr)
!~     call MatDuplicate(p_zint_out, MAT_SHARE_NONZERO_PATTERN, p_tmp1_int3, ierr)  
!~     call MatDuplicate(p_zint_out, MAT_SHARE_NONZERO_PATTERN, p_tmp1_int4, ierr)  
  
    call MatZeroEntries(p_zint_out, ierr)
    
    int_counter = 0
    if (l_output_progress) then
      write (pstr_out, fmt='(A,2e24.12)') "integrate: ",a, b
      call petsc_print_master()
    end if    
    call adaptive_int2(f, a, b, 0, eps, max_error)
  
    deallocate(nodes, weight1, weight2, error_tracker)

    call MatDestroy(p_tmp1_int1, ierr)
!~     call MatDestroy(p_tmp1_int2, ierr)
!~     call MatDestroy(p_tmp1_int3, ierr)  

  
  end subroutine integrate
  
  
  subroutine evalqr_2(f, a, b, p_zint1, error_nested)
#include <petsc/finclude/petsc.h>
    use petscmat
    use petsc_mod
    use petsc_wrapper
    use globals, only : nint_order, ncell_c, wkp_c, dlat_c, iik
    use ft_mod
    use kinds  
    implicit none
    
    integer, external :: f
    real(dp) :: a, b, norm
    real(dp), optional :: error_nested
    Mat :: p_zint1
    
    integer :: ierr, i, i1 ,i2
    PetscScalar, allocatable :: ws1(:), ws2(:)
    real(dp), allocatable :: xs(:)
    logical :: l_nested_error
    Mat :: p_zint1_nested
    Mat, allocatable :: p_tmp1(:)
    real(dp) :: timing_local
    integer(8) :: counti, count_rate, countf
    PetscScalar :: p_fac
    
    
    l_nested_error = .false.
    if (present(error_nested)) l_nested_error = .true.
    
    allocate (ws1(nint_order), ws2(nint_order), xs(nint_order), p_tmp1(1))
    
    ws1 = weight1*(b - a)*0.5d0
    ws2 = weight2*(b - a)*0.5d0
    xs = nodes*(b - a)*0.5d0 + (a + b)*0.5d0
    
    call MatZeroEntries(p_zint1, ierr)
    if (l_nested_error) then
      call MatDuplicate(p_zint1, MAT_SHARE_NONZERO_PATTERN, p_zint1_nested, ierr)
      call MatZeroEntries(p_zint1_nested, ierr)
    end if 

    timing_local = 0d0
    if (l_output_progress) then
      write (pstr_out, fmt='(i3,e11.4)') 0, timing_local; call PetscPrintf(PETSC_COMM_WORLD, trim(pstr_out), ierr)
    end if
    
    do i = 1, nint_order   
     
      call system_clock(counti, count_rate)
      
      ierr = f(xs(i))
      
      call MatAXPY(p_zint1, ws1(i), p_tmpcc1, SAME_NONZERO_PATTERN, ierr)
      
!~       call MatNorm(p_tmpcc1, NORM_FROBENIUS, error_nested, ierr)
!~       if (inode.eq.0) write(0, fmt = '(i8,A,2e24.12, i8)') nsubdivs," XXX ",xs(i), error_nested, i
      
      if (present(error_nested)) then
        if (ws2(i) .ne. 0d0) then
          call MatAXPY(p_zint1_nested, ws2(i), p_tmpcc1, SAME_NONZERO_PATTERN, ierr)
        end if
      end if


      int_counter = int_counter + 1
      
      call system_clock(countf)
      timing_local = timing_local + real(countf - counti, 8)/real(count_rate, 8)
      if (l_output_progress) then
        write (pstr_out, fmt='(A,i3,e11.4)') repeat(ACHAR(8), 14), i, timing_local; call PetscPrintf(PETSC_COMM_WORLD, trim(pstr_out), ierr)
      end if
      
    end do
              
    if (l_nested_error) then
      if (error_reim_select .eq. -1) then
      
        call MatAXPY(p_zint1_nested, p_minus1, p_zint1, SAME_NONZERO_PATTERN, ierr)
        call MatNorm(p_zint1_nested, NORM_FROBENIUS, error_nested, ierr)
!~         call MatNorm(p_zint1, NORM_FROBENIUS, norm, ierr)
!~         error_nested = error_nested / norm
        
      else
      
        do i1 = -ncell_c(1), ncell_c(1)
          do i2 = -ncell_c(2), ncell_c(2)
            call MatZeroEntries(p_dmat_error(i1, i2), ierr)
          end do
        end do     
        
        p_tmp1(1) = p_zint1_nested
        call fourier_back_add(p_tmp1, p_dmat_error(:, :), kp_c(1:2, iik:iik), &
       &  wkp_c(iik:iik), 1, dlat_c, ncell_c(1), ncell_c(2), error_reim_select, 0)        
        call MatCopy(p_zint1, p_zint1_nested, SAME_NONZERO_PATTERN, ierr)
        p_fac = -1d0
        call fourier_back_add(p_tmp1, p_dmat_error(:, :), kp_c(1:2, iik:iik), &
       &  wkp_c(iik:iik), 1, dlat_c, ncell_c(1), ncell_c(2), error_reim_select, 0, p_fac)        
       
        error_nested = 0d0
        do i1 = -ncell_c(1), ncell_c(1)
          do i2 = -ncell_c(2), ncell_c(2)
            call MatRealPart(p_dmat_error(i1, i2), ierr)
            call MatNorm(p_dmat_error(i1, i2), NORM_FROBENIUS, norm, ierr)
            error_nested = max(norm, error_nested)
          end do
        end do        
         
      end if
      
      call MatDestroy(p_zint1_nested, ierr)
    end if
    
  end subroutine evalqr_2
   

  recursive subroutine adaptive_int2(f, a, b, level, eps, est_error)
#include <petsc/finclude/petsc.h>
    use petscmat
    use petsc_mod
    use petsc_wrapper
    use globals, only : maxsub
    use error_handler
    use kinds

    implicit none  
  
    integer, external :: f
    
    integer :: level
    real(dp) :: a, b, eps, est_error
    logical :: l_subfile, l_finished

    real(dp) :: h, error_ab, x1, x2, error_nested
    integer :: ierr, nsub, i

    nsubdivs = nsubdivs + 1

    if (level > maxsub) then
      write (errormsg, fmt = '(A, 2i8)') "level > maxsub ", level, maxsub
      call error()
    end if

    nsub = 2
    h =  (b - a) / real(nsub, 8)
    h_min = min(h, h_min)
        
    call evalqr_2(f, a, b, p_tmp1_int1, error_nested)          

!~ use subdiv as error estimate
!~     x1 = a    
!~     call MatZeroEntries(p_tmp1_int3, ierr)
!~     do i = 1, nsub
!~       x2 = x1 + h
!~       call evalqr_2(f, x1, x2, p_tmp1_int2)     
!~       call MatAYPX(p_tmp1_int3, p_one, p_tmp1_int2, SAME_NONZERO_PATTERN, ierr)
!~       x1 = x1 + h        
!~     end do
    !~ estimate error |integral_ab-(integral_ah+integral_hb)|
!~     call MatAXPY(p_tmp1_int1, p_minus1, p_tmp1_int3, SAME_NONZERO_PATTERN, ierr)     
!~     call MatNorm(p_tmp1_int1, NORM_FROBENIUS, est_error, ierr)
!~     error_ab = est_error
!~     call MatCopy(p_tmp1_int3, p_tmp1_int1, SAME_NONZERO_PATTERN, ierr)
!~ ----    

    error_ab = error_nested
    est_error = error_nested
    
    if (l_output_progress) then
!~       pstr_out=""
!~       call petsc_print_master()
      call MatNorm(p_tmp1_int1, NORM_FROBENIUS, est_error, ierr)
      write (pstr_out, fmt='(A,2i5, 5e16.8, i8)') " info: ", nsubdivs, level, a, b, h, error_ab, error_nested, int_counter
      call petsc_print_master()
    end if        
    
    
    if ((error_ab > eps) .and. (h .ge. 1d-7)) then
      x1 = a    
      do i = 1, nsub
        x2 = x1 + h
        call adaptive_int2(f, x1, x2, level + 1, eps, est_error)
        x1 = x1 + h        
      end do
    else 
      call MatAXPY(p_zint_out, p_one, p_tmp1_int1, SAME_NONZERO_PATTERN, ierr) 
    end if
  
    
  
  end subroutine adaptive_int2

  subroutine adaptive_int3(f, a, b, p_zint1, p_zint2, nsub, eps, abserr, subfile)
#include <petsc/finclude/petsc.h>
    use petscmat
    use petsc_mod
    use petsc_wrapper
    use kinds

    implicit none

    integer, external :: f
    integer :: nsub, i1, i2, i3, norder
    real(dp) :: abserr, a, b, eps, abserr_tot
    character(*), optional :: subfile
    logical :: l_subfile, l_finished
    Mat :: p_zint1, p_zint2

    integer :: ierr, isub, n, nlow, iunit, i, ii(1), istart, nsub_used, kk, nsub_used2, icheck

    real(dp), allocatable :: nodes(:), w1(:), w2(:), subint(:, :), suberr(:), suberr2(:), &
                             subdiv(:, :)
    real(dp) :: m, subcurrent(2, 3), globalnorm, errorg, error2, error_1, error_2, &
   &  norm, sub(2), fac

    int_counter = 0

    call get_quadrule(integrator, nodes, w1, w2, nint_order, nlow)

!~       n=size(zint1,1)

! setup necessary temporary matrices
    call MatDuplicate(p_zint1, MAT_SHARE_NONZERO_PATTERN, p_tmp1_int1, ierr)
    call MatDuplicate(p_zint1, MAT_SHARE_NONZERO_PATTERN, p_tmp1_int2, ierr)
    call MatDuplicate(p_zint1, MAT_SHARE_NONZERO_PATTERN, p_tmp1_int3, ierr)

! initialize matrices
    call MatZeroEntries(p_zint1, ierr)
    if (ldouble_contour .and. lnoneq_int) call MatZeroEntries(p_zint2, ierr)

    if (ldouble_contour .and. lnoneq_int) then
      call MatDuplicate(p_zint1, MAT_SHARE_NONZERO_PATTERN, p_tmp2_int1, ierr)
      call MatDuplicate(p_zint1, MAT_SHARE_NONZERO_PATTERN, p_tmp2_int2, ierr)
      call MatDuplicate(p_zint1, MAT_SHARE_NONZERO_PATTERN, p_tmp2_int3, ierr)
    end if

    allocate (subint(2, nsub), suberr(nsub), suberr2(nsub), subdiv(2, nsub))
    subint = 0d0
    suberr = 0d0
    suberr2 = 0d0
    subdiv = 0d0
    l_subfile = .false.
    if (present(subfile)) inquire (file=trim(subfile), exist=l_subfile)
    if (l_output_progress) then
      write (pstr_out, fmt='(A,l)') "read subdivision from file "//trim(subfile)//" ", l_subfile; call petsc_print_master()
    end if

    if (l_subfile) then
      suberr = 0d0
      if (l_ionode) then
        open (newunit=iunit, file=trim(subfile), action="read", status="old")
        read (iunit, *) isub
        do i1 = 1, isub
          read (iunit, *) i3, subint(1:2, i1), suberr(i1), subdiv(1:2, i1)
        end do
        close (iunit)
      end if
      call MPI_bcast(isub, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
      call MPI_bcast(subint(1:2, 1:isub), 2*isub, MPI_DOUBLE_PRECISION, 0, PETSC_COMM_WORLD, ierr)
      call MPI_bcast(suberr(1:isub), isub, MPI_DOUBLE_PRECISION, 0, PETSC_COMM_WORLD, ierr)
      call MPI_bcast(subdiv(1:2, 1:isub), 2*isub, MPI_DOUBLE_PRECISION, 0, PETSC_COMM_WORLD, ierr)
      do i2 = isub, 1, -1
!~           i1=isub-i2+1
        if (l_output_progress) then
          write (pstr_out, fmt='(A,i8,3e24.12)') "doing sub ", i2, subint(1:2, i2), suberr(i2); call petsc_print_master(.false.)
        end if
        call evalqr(f, subint(1, i2), subint(2, i2), nodes, w1, w2, p_tmp1_int1, p_tmp2_int1, nint_order, nlow, error_1, error_2, subdiv(1, i2), subdiv(2, i2))
        suberr(i2) = max(error_1, error_2)
        if (l_output_progress) then
          write (pstr_out, fmt='(2X,e24.12)') suberr(i2); call petsc_print_master()
        end if
        subcurrent(1:2, 1) = subint(1:2, i2)
        call MatAXPY(p_zint1, p_one, p_tmp1_int1, SAME_NONZERO_PATTERN, ierr)
        if (ldouble_contour .and. lnoneq_int) call MatAXPY(p_zint2, p_one, p_tmp2_int1, SAME_NONZERO_PATTERN, ierr)       !     p_zint2=p_tmp2_int1
      end do
      call MatNorm(p_zint1, NORM_FROBENIUS, norm, ierr)
      globalnorm = norm
      if (ldouble_contour .and. lnoneq_int) then
        call MatNorm(p_zint2, NORM_FROBENIUS, norm, ierr)
        globalnorm = (globalnorm + norm)*0.5d0
      end if
      subcurrent(1:2, 1) = subint(1:2, isub)
      subcurrent(1:2, 2) = 0d0
      subcurrent(1:2, 3) = 0d0
      suberr2 = suberr
      call hpsort_aa(suberr, subint, .true.)
      call hpsort_aa(suberr2, subdiv, .true.)
      isub = isub + 3
    else
      isub = 1
      subint(1, isub) = a
      subint(2, isub) = b
      suberr = 0d0

      i = 0
      call evalqr(f, subint(1, 1), subint(2, 1), nodes, w1, w2, p_tmp1_int1, p_tmp2_int1, nint_order, nlow, error_1, error_2, subdiv(1, 1), subdiv(2, 1))
      errorg = 0d0
      suberr(1) = max(error_1, error_2)
      subcurrent(1:2, 1) = subint(1:2, 1)
      call MatAXPY(p_zint1, p_one, p_tmp1_int1, SAME_NONZERO_PATTERN, ierr)
      if (ldouble_contour .and. lnoneq_int) call MatAXPY(p_zint2, p_one, p_tmp2_int1, SAME_NONZERO_PATTERN, ierr)       !     p_zint2=p_tmp2_int1

      call MatNorm(p_zint1, NORM_FROBENIUS, norm, ierr)
      globalnorm = norm
      if (ldouble_contour .and. lnoneq_int) then
        call MatNorm(p_zint2, NORM_FROBENIUS, norm, ierr)
        globalnorm = (globalnorm + norm)*0.5d0
      end if
      abserr = suberr(1)
    end if
!~       suberr(1)=0d0

    if (l_output_progress) then
      call PetscPrintf(PETSC_COMM_WORLD, NEW_LINE('A'), ierr)
    end if
    do

      abserr = maxval(suberr) !/globalnorm
      if (l_output_progress) then
        write (pstr_out, fmt='(i4,7e16.8,e10.3,i8)') &
!~           isub, suberr(1:3)/globalnorm, subint(1:2, 1), globalnorm, abserr, eps, int_counter; call petsc_print_master()
          isub, suberr(1:3), subint(1:2, 1), globalnorm, abserr, eps, int_counter; call petsc_print_master()
      end if
      if ((abserr .gt. eps) .and. (isub*3 + 3 .le. nsub)) then
!~           error2=error2-suberr(1)
        suberr(1) = 0d0
        m = subint(2, 1) - subint(1, 1)
        m = m/3d0
        i1 = isub + 1
        i2 = isub + 2
        i3 = isub + 3
        subint(1, i1) = subint(1, 1)
        subint(2, i1) = subdiv(1, 1) !subint(1,1)+1d0*m

        subint(1, i2) = subdiv(1, 1) !subint(1,1)+1d0*m
        subint(2, i2) = subdiv(2, 1) !subint(1,1)+2d0*m

        subint(1, i3) = subdiv(2, 1) !subint(1,1)+2d0*m
        subint(2, i3) = subint(2, 1)
        suberr(1) = 0d0
        isub = isub + 3
!~           write(0,fmt='(A,3e24.12)') "new devision ",subint(1:2,i1),subint(2,i1)-subint(1,i1)
!~           write(0,fmt='(A,3e24.12)') "new devision ",subint(1:2,i2),subint(2,i2)-subint(1,i2)
!~           write(0,fmt='(A,3e24.12)') "new devision ",subint(1:2,i3),subint(2,i3)-subint(1,i3)
      else
        exit
      end if

      if (all(subcurrent(1:2, 1) - subint(1:2, 1) .eq. 0)) then
        call MatAXPY(p_zint1, p_minus1, p_tmp1_int1, SAME_NONZERO_PATTERN, ierr)
        if (ldouble_contour .and. lnoneq_int) call MatAXPY(p_zint2, p_minus1, p_tmp2_int1, SAME_NONZERO_PATTERN, ierr)
      else if (all(subcurrent(1:2, 2) - subint(1:2, 1) .eq. 0)) then
        call MatAXPY(p_zint1, p_minus1, p_tmp1_int2, SAME_NONZERO_PATTERN, ierr)
        if (ldouble_contour .and. lnoneq_int) call MatAXPY(p_zint2, p_minus1, p_tmp2_int2, SAME_NONZERO_PATTERN, ierr)
      else if (all(subcurrent(1:2, 3) - subint(1:2, 1) .eq. 0)) then
        call MatAXPY(p_zint1, p_minus1, p_tmp1_int3, SAME_NONZERO_PATTERN, ierr)
        if (ldouble_contour .and. lnoneq_int) call MatAXPY(p_zint2, p_minus1, p_tmp2_int3, SAME_NONZERO_PATTERN, ierr)
      else
        call evalqr(f, subint(1, 1), subint(2, 1), nodes, w1, w2, p_tmp1_int1, p_tmp2_int1, nint_order, nlow, error_1, error_2, subdiv(1, 1), subdiv(2, 1))
        call MatAXPY(p_zint1, p_minus1, p_tmp1_int1, SAME_NONZERO_PATTERN, ierr)
        if (ldouble_contour .and. lnoneq_int) call MatAXPY(p_zint2, p_minus1, p_tmp2_int1, SAME_NONZERO_PATTERN, ierr)
      end if

! a..a+1/3
      call evalqr(f, subint(1, i1), subint(2, i1), nodes, w1, w2, p_tmp1_int1, p_tmp2_int1, nint_order, nlow, error_1, error_2, subdiv(1, i1), subdiv(2, i1))
      suberr(i1) = max(error_1, error_2)
      call MatAXPY(p_zint1, p_one, p_tmp1_int1, SAME_NONZERO_PATTERN, ierr)
      if (ldouble_contour .and. lnoneq_int) call MatAXPY(p_zint2, p_one, p_tmp2_int1, SAME_NONZERO_PATTERN, ierr)

! a+1/3*d..a+2/3*d
      call evalqr(f, subint(1, i2), subint(2, i2), nodes, w1, w2, p_tmp1_int2, p_tmp2_int2, nint_order, nlow, error_1, error_2, subdiv(1, i2), subdiv(2, i2))
      suberr(i2) = max(error_1, error_2)
      call MatAXPY(p_zint1, p_one, p_tmp1_int2, SAME_NONZERO_PATTERN, ierr)
      if (ldouble_contour .and. lnoneq_int) call MatAXPY(p_zint2, p_one, p_tmp2_int2, SAME_NONZERO_PATTERN, ierr)

! a+2/3d..a+3/3*d=b
      call evalqr(f, subint(1, i3), subint(2, i3), nodes, w1, w2, p_tmp1_int3, p_tmp2_int3, nint_order, nlow, error_1, error_2, subdiv(1, i3), subdiv(2, i3))
      suberr(i3) = max(error_1, error_2)
      call MatAXPY(p_zint1, p_one, p_tmp1_int3, SAME_NONZERO_PATTERN, ierr)
      if (ldouble_contour .and. lnoneq_int) call MatAXPY(p_zint2, p_one, p_tmp2_int3, SAME_NONZERO_PATTERN, ierr)

      call MatNorm(p_zint1, NORM_FROBENIUS, norm, ierr)
      globalnorm = norm
      if (ldouble_contour .and. lnoneq_int) then
        call MatNorm(p_zint2, NORM_FROBENIUS, norm, ierr)
        globalnorm = (globalnorm + norm)*0.5d0
      end if

      if (l_output_progress) call PetscPrintf(PETSC_COMM_WORLD, NEW_LINE('A'), ierr)
      if (l_output_progress) then
        write (pstr_out, fmt='(i8,6e24.12)') i1, subint(1:2, i1),&
        &subint(2, i1) - subint(1, i1), suberr(i1), subdiv(1:2, i1)
        call petsc_print_master()
      end if
      if (l_output_progress) then
        write (pstr_out, fmt='(i8,6e24.12)') i2, subint(1:2, i2),&
        &subint(2, i2) - subint(1, i2), suberr(i2), subdiv(1:2, i2)
        call petsc_print_master()
      end if
      if (l_output_progress) then
        write (pstr_out, fmt='(i8,6e24.12)') i3, subint(1:2, i3),&
        &subint(2, i3) - subint(1, i3), suberr(i3), subdiv(1:2, i3)
        call petsc_print_master()
      end if

      if ((abs((subint(2, i1) - subint(1, i1))) .le. 1d-9).and. &
     & suberr(i1).gt.eps)then
          if (l_ionode) write(0, fmt='(A,4e24.12)') "subdiv very small", abs(subint(2, i1) - subint(1, i1)), subint(2, i1), subint(1, i1),suberr(i1)
!~         suberr(i1) = eps
      end if
      if ((abs((subint(2, i2) - subint(1, i2))) .le. 1d-9).and. &
     & suberr(i2).gt.eps) then
        if (l_ionode) write(0, fmt='(A,4e24.12)') "subdiv very small", abs(subint(2, i2) - subint(1, i2)), subint(2, i2), subint(1, i2),suberr(i2)
!~         suberr(i2) = eps
      end if
      if ((abs((subint(2, i3) - subint(1, i3))) .le. 1d-9).and. &
     & suberr(i3).gt.eps) then
        if (l_ionode) write(0, fmt='(A,4e24.12)') "subdiv very small", abs(subint(2, i3) - subint(1, i3)), subint(2, i3), subint(1, i3),suberr(i3)
!~         suberr(i3) = eps
      end if

      subint(1:2, 1) = -1d0
      suberr(1) = -1d0
      subcurrent(1:2, 1) = subint(1:2, i1)
      subcurrent(1:2, 2) = subint(1:2, i2)
      subcurrent(1:2, 3) = subint(1:2, i3)
      suberr2 = suberr
      call hpsort_aa(suberr, subint, .true.)
      call hpsort_aa(suberr2, subdiv, .true.)

    end do
    

    call MatDestroy(p_tmp1_int1, ierr)
    call MatDestroy(p_tmp1_int2, ierr)
    call MatDestroy(p_tmp1_int3, ierr)
    if (ldouble_contour .and. lnoneq_int) then
      call MatDestroy(p_tmp2_int1, ierr)
      call MatDestroy(p_tmp2_int2, ierr)
      call MatDestroy(p_tmp2_int3, ierr)
    end if

    if (present(subfile)) then
      if (l_ionode) then
!~           suberr=-suberr

!~           call hpsort_aa (suberr(1:i3),subint_s(1:2,1:i3),.true.)
!~           call hpsort_aa (suberr2(1:i3),subdiv(1:2,1:i3),.true.)
!~           suberr=-suberr

        ii = minloc(abs(suberr))
        i3 = ii(1) - 1
        if (i3 .eq. 1) return
        do i1 = 1, i3 - 1
          do i2 = i1 + 1, i3
            if (subint(1, i1) .ge. subint(1, i2)) then
              sub = subint(1:2, i1)
              subint(1:2, i1) = subint(1:2, i2)
              subint(1:2, i2) = sub
              sub = subdiv(1:2, i1)
              subdiv(1:2, i1) = subdiv(1:2, i2)
              subdiv(1:2, i2) = sub
              error2 = suberr(i1)
              suberr(i1) = suberr(i2)
              suberr(i2) = error2
            end if
          end do
        end do
        
!~         do i1 = 1, i3
!~           write (pstr_out, fmt='(i4,5e24.12,A)') i1, subint(1:2, i1),suberr(i1), subdiv(1:2, i1), " CHECK"
!~           call petsc_print_master()
!~         end do


        kk = 0
        istart = 2
        nsub_used = i3
        nsub_used2 = i3
        do 
          kk = kk + 1
          i3 = nsub_used
          l_finished = .true.
          fac = 1d-4
          do i=istart,i3 -1    
            if (nsub_used.le.i) exit
            if (abs(suberr(i)).le.eps*fac) then
              if (abs(suberr(i-1)).le.abs(suberr(i+1))) then
                i2 = i - 1
              else 
                i2 = i + 1
              end if
              if ((abs(suberr(i))+abs(suberr(i2))).le.eps*fac) then        
                l_finished = .false.
                if (i2 .lt. i) then
                  suberr(i2) = suberr(i)+suberr(i2)            
                  subint(2,i2) = subint(2,i)          
                  subint(1:2,i:nsub_used-1)=subint(1:2,i+1:nsub_used)
                  suberr(i:nsub_used-1)=suberr(i+1:nsub_used)
                  nsub_used = nsub_used - 1            
                else
                  suberr(i) = suberr(i)+suberr(i2)
                  subint(2,i) = subint(2,i2)
                  subint(1:2,i2:nsub_used-1)=subint(1:2,i2+1:nsub_used)
                  suberr(i2:nsub_used-1)=suberr(i2+1:nsub_used)
                  nsub_used = nsub_used - 1
                end if   
!~                 write(6,fmt='(A,4i4)') "-x-x-x-x-",kk,i,i2,nsub_used
                do icheck=1,nsub_used
                  if (suberr(i).le.0) cycle
!~                   write(6,fmt='(i4,3e24.12)') icheck, subint(1:2,icheck), suberr(icheck)
                  if (icheck.ge.2) then
                    if (subint(2,icheck-1).ne.subint(1,icheck)) then
                      if (l_ionode) write(0,*) "disconnected segment ",icheck,subint(1:2,icheck-1),subint(1:2,icheck)
                      stop
                    end if
                  end if            
                end do             
!~                 write(6,*) "-x-x-x-x-",i,nsub_used
              end if
            end if    
          end do
          if (l_finished) exit
        end do


        suberr2 = suberr
        call hpsort_aa(suberr(1:nsub_used), subint(1:2, 1:nsub_used), .true.)
        call hpsort_aa(suberr2(1:nsub_used), subdiv(1:2, 1:nsub_used), .true.)
!~         ii = minloc(abs(suberr)) ! that's probably not necessary anymore. remove ?
!~         i3 = ii(1) - 1
        open (newunit=iunit, file=trim(subfile), action="write", status="replace")
        write (iunit, fmt='(2i8, 5e24.12)') nsub_used, nsub_used2, globalnorm, a, b
        do i1 = nsub_used, 1, -1
          write (iunit, fmt='(i8,5es45.24e5)') i1, subint(1:2, i1), suberr(i1), subdiv(1:2, i1)
        end do
        close (iunit)
      end if
    end if

!~       nullify(p_tmp1_int1_ptr)
!~       nullify(p_tmp1_int2_ptr)
!~       nullify(p_tmp1_int3_ptr)
!~       nullify(p_tmp2_int1_ptr)
!~       nullify(p_tmp2_int2_ptr)
!~       nullify(p_tmp2_int3_ptr)

  end subroutine adaptive_int3

end module integrator_mod
