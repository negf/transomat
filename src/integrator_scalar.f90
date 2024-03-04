module integrator_scalar
  use kinds
  implicit none

  complex(dp) :: zint_scalar_output

contains

  subroutine evalqr_scalar(f, a, b, x, w1, w2, p_zint, p_zint2, nhigh, nlow, errout, errout2, d1, d2)
    use petsc_mod
    use kinds
    use mathstuff
    use globals, only: l_output_progress, int_counter

    implicit none

    real(dp) ::a, b, errout, errout2
    real(dp), optional :: d1, d2
    real(dp), allocatable :: x(:), w1(:), w2(:)
    complex(dp) :: p_zint, p_zint2
    integer :: nhigh, nlow, ierr
    integer, external :: f

    complex(dp) :: p_zintg
    integer :: i, n, j, i_low, i_high, i_max, k1, k2
    real(dp) :: xx, norm, m1, m2, timing_local
    integer(8) :: counti, count_rate, countf

    complex(dp), allocatable :: ws1(:), ws2(:)

    allocate (ws1(nhigh), ws2(nhigh))

    ws1 = w1*(b - a)*0.5d0
    ws2 = w2*(b - a)*0.5d0

    p_zint = 0d0
    p_zintg = 0d0

    timing_local = 0d0
    if (l_output_progress) then
      write (pstr_out, fmt='(i3,e16.6)') 0, timing_local; call PetscPrintf(PETSC_COMM_WORLD, trim(pstr_out), ierr)
    end if
    do i = 1, nhigh
      call system_clock(counti, count_rate)
      xx = x(i)*(b - a)*0.5d0 + (a + b)*0.5d0

      ierr = f(xx)

      p_zint = p_zint + ws1(i)*zint_scalar_output

      if (ws2(i) .ne. 0d0) then
        p_zintg = p_zintg + ws2(i)*zint_scalar_output
      end if

      int_counter = int_counter + 1
      call system_clock(countf)
      timing_local = timing_local + real(countf - counti, 8)/real(count_rate, 8)
      if (l_output_progress) then
        write (pstr_out, fmt='(A,i3,e16.6)') repeat(ACHAR(8), 19), i, timing_local; call PetscPrintf(PETSC_COMM_WORLD, trim(pstr_out), ierr)
      end if
    end do
!~       call PetscPrintf(PETSC_COMM_WORLD,"          ",ierr)

    p_zintg = p_zintg - p_zint

    errout = abs(p_zintg)
    errout2 = errout

    if (.not. present(d2)) return

    d1 = a + (b - a)/3d0
    d2 = a + (b - a)*2d0/3d0

  end subroutine evalqr_scalar

  subroutine adaptive_int3_scalar(f, a, b, p_zint1, p_zint2, nsub, eps, abserr)
    use petsc_mod
    use integrator_mod
    use globals
    use kinds

    implicit none

    integer, external :: f
    integer :: nsub, i1, i2, i3, norder
    real(dp) :: abserr, a, b, eps, abserr_tot
    complex(dp) :: p_zint1, p_zint2, s_tmp1_int1, s_tmp1_int2, s_tmp1_int3,&
    &s_tmp2_int1, s_tmp2_int2, s_tmp2_int3

    integer :: ierr, isub, n, nlow, iunit, i, ii(1)

    real(dp), allocatable :: nodes_int(:), w1(:), w2(:), subint(:, :), suberr(:), suberr2(:), &
                             subdiv(:, :)
    real(dp) :: m, subcurrent(2, 3), globalnorm, errorg, error2, error_1, error_2, norm, sub(2)

    int_counter = 0
    call get_quadrule(integrator, nodes_int, w1, w2, nint_order, nlow)

    p_zint1 = 0d0
    p_zint2 = 0d0

    allocate (subint(2, nsub), suberr(nsub), suberr2(nsub), subdiv(2, nsub))
    subint = 0d0
    suberr = 0d0
    suberr2 = 0d0
    subdiv = 0d0

    isub = 1
    subint(1, isub) = a
    subint(2, isub) = b
    suberr = 0d0

    i = 0
    call evalqr_scalar(f, subint(1, 1), subint(2, 1), nodes_int, w1, w2, s_tmp1_int1, s_tmp2_int1, nint_order, nlow, error_1, error_2, subdiv(1, 1), subdiv(2, 1))
    errorg = 0d0
    suberr(1) = max(error_1, error_2)
    subcurrent(1:2, 1) = subint(1:2, 1)
    p_zint1 = p_zint1 + s_tmp1_int1
    globalnorm = abs(p_zint1)
    abserr = suberr(1)

    if (l_output_progress) then
      call PetscPrintf(PETSC_COMM_WORLD, NEW_LINE('A'), ierr)
    end if
    do

      abserr = maxval(suberr)/globalnorm
      if (l_output_progress) then
        write (pstr_out, fmt='(i4,7e16.8,e10.3,i8)') &
          isub, suberr(1:3)/globalnorm, subint(1:2, 1), globalnorm, abserr, eps, int_counter; call petsc_print_master()
      end if
      if ((abserr .gt. eps) .and. (isub*3 + 3 .le. nsub)) then

        suberr(1) = 0d0
        m = subint(2, 1) - subint(1, 1)
        m = m/3d0
        i1 = isub + 1
        i2 = isub + 2
        i3 = isub + 3
        subint(1, i1) = subint(1, 1)
        subint(2, i1) = subdiv(1, 1)

        subint(1, i2) = subdiv(1, 1)
        subint(2, i2) = subdiv(2, 1)

        subint(1, i3) = subdiv(2, 1)
        subint(2, i3) = subint(2, 1)
        suberr(1) = 0d0
        isub = isub + 3

      else
        exit
      end if

      if (all(subcurrent(1:2, 1) - subint(1:2, 1) .eq. 0)) then
        p_zint1 = p_zint1 - s_tmp1_int1
      else if (all(subcurrent(1:2, 2) - subint(1:2, 1) .eq. 0)) then
        p_zint1 = p_zint1 - s_tmp1_int2
      else if (all(subcurrent(1:2, 3) - subint(1:2, 1) .eq. 0)) then
        p_zint1 = p_zint1 - s_tmp1_int3
      else
        call evalqr_scalar(f, subint(1, 1), subint(2, 1), nodes_int, w1, w2, s_tmp1_int1, s_tmp2_int1, nint_order, nlow, error_1, error_2, subdiv(1, 1), subdiv(2, 1))
        p_zint1 = p_zint1 - s_tmp1_int1
      end if

! a..a+1/3
      call evalqr_scalar(f, subint(1, i1), subint(2, i1), nodes_int, w1, w2, s_tmp1_int1, s_tmp2_int1, nint_order, nlow, error_1, error_2, subdiv(1, i1), subdiv(2, i1))
      suberr(i1) = max(error_1, error_2)
      p_zint1 = p_zint1 + s_tmp1_int1

! a+1/3*d..a+2/3*d
      call evalqr_scalar(f, subint(1, i2), subint(2, i2), nodes_int, w1, w2, s_tmp1_int2, s_tmp2_int2, nint_order, nlow, error_1, error_2, subdiv(1, i2), subdiv(2, i2))
      suberr(i2) = max(error_1, error_2)
      p_zint1 = p_zint1 + s_tmp1_int2

! a+2/3d..a+3/3*d=b
      call evalqr_scalar(f, subint(1, i3), subint(2, i3), nodes_int, w1, w2, s_tmp1_int3, s_tmp2_int3, nint_order, nlow, error_1, error_2, subdiv(1, i3), subdiv(2, i3))
      suberr(i3) = max(error_1, error_2)
      p_zint1 = p_zint1 + s_tmp1_int3

      norm = abs(p_zint1)
      globalnorm = norm

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

      if (abs((subint(2, i1) - subint(1, i1))) .le. 1d-3) then
        suberr(i1) = 0d0
      end if
      if (abs((subint(2, i2) - subint(1, i2))) .le. 1d-3) then
        suberr(i2) = 0d0
      end if
      if (abs((subint(2, i3) - subint(1, i3))) .le. 1d-3) then
        suberr(i3) = 0d0
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

!~     write(pstr_out,fmt='(A,i8)') "ints ",int_counter ; call petsc_print_master()

  end subroutine adaptive_int3_scalar

end module
