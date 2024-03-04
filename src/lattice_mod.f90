module lattice_mod
  implicit none
contains

  subroutine read_cell(brav, iunit)
#include <petsc/finclude/petscmat.h>
    use petscmat
    use kinds
    use globals, only: eps_real
!~       use blas95
    use error_handler
    implicit none

    real(dp), intent(out) :: brav(3, 3)
    integer, intent(in) :: iunit

    integer :: i, j, ierr
    character(256) :: str

    do i = 1, 3
      read (iunit, *) brav(1:3, i)
      write (str, fmt='(3f24.12)') brav(1:3, i)
      call PetscPrintf(PETSC_COMM_WORLD, trim(str)//New_Line('A'), ierr)
    end do

    do i = 1, 2
      do j = i + 1, 3
      if (.not. abs(dot_product(brav(1:3, i), brav(1:3, j))) .lt. eps_real) then
        write (errormsg, *) "super cell vector ", i, " and ", j, " are not orthogonal"
        call error()
      end if
      end do
    end do

  end subroutine read_cell

  subroutine read_kpoint_file(kps, wkps, nkp)
    use petscmat
    use kinds
    use globals, only: rlat_l
    use error_handler
    implicit none

    real(dp), allocatable :: kps(:, :)
    integer, allocatable :: wkps(:)
    integer :: nkp

    integer :: iunit, ik, ierr
    character(256) :: instr, frac, str
    logical :: lfrac

    open (newunit=iunit, file="k_points.dat", status="old", action="read")

    read (iunit, fmt='(A)') instr
    call PetscPrintf(PETSC_COMM_WORLD, trim(instr)//New_Line('A'), ierr)
    read (instr, *) nkp, frac
    frac = adjustl(frac)
    lfrac = .false.
    if (index(frac, "frac") .ge. 1) lfrac = .true.

    allocate (kps(2, nkp), wkps(nkp), stat=ierr)
    if (ierr .ne. 0) then
      write (errormsg, *) "allocation error ", ierr
      call error()
    end if

    write (str, fmt='(A)') "k-points and weights from file:"
    call PetscPrintf(PETSC_COMM_WORLD, trim(str)//New_Line('A'), ierr)
    do ik = 1, nkp
      read (iunit, *) kps(1:2, ik), wkps(ik)
      if (lfrac) kps(1:2, ik) = kps(1, ik)*rlat_l(1:2, 1) + kps(2, ik)*rlat_l(1:2, 2)
      write (str, fmt='(2e24.12,i8)') kps(1:2, ik), wkps(ik)
      call PetscPrintf(PETSC_COMM_WORLD, trim(str)//New_Line('A'), ierr)
    end do

    close (iunit)

  end subroutine read_kpoint_file

  subroutine get_rlat(dlat, rlat, kps, wkps, nkp, nkx, nky, nkz, makeks)
    use petscmat
    use kinds
    use globals, only: pi, ksym
    use mathstuff, only: vector_product
    use error_handler
!~       use blas95
    implicit none

    real(dp), intent(in) :: dlat(3, 3)
    integer, intent(in) :: nkx, nky, nkz
    logical, intent(in) :: makeks
    real(dp), allocatable :: kps(:, :)
    integer, allocatable :: wkps(:)
    real(dp), intent(out) :: rlat(3, 3)
    integer :: nkp
    integer :: i, j, nxk, nyk, nk1, nk2, nk3, ierr, n1, &
      ik1, ik2, ik3, i1, in1, in2, in3
    integer, allocatable :: ikps(:, :)
    real(dp) :: dk(3), deltak(3, 3)
    character(256) :: str
    logical :: l_withgamma

    rlat(1:3, 1) = vector_product(dlat(1:3, 2), dlat(1:3, 3))
    rlat(1:3, 1) = rlat(1:3, 1)/(dot_product(dlat(1:3, 1), rlat(1:3, 1)))

    rlat(1:3, 2) = vector_product(dlat(1:3, 3), dlat(1:3, 1))
    rlat(1:3, 2) = rlat(1:3, 2)/(dot_product(dlat(1:3, 2), rlat(1:3, 2)))

    rlat(1:3, 3) = vector_product(dlat(1:3, 1), dlat(1:3, 2))
    rlat(1:3, 3) = rlat(1:3, 3)/(dot_product(dlat(1:3, 3), rlat(1:3, 3)))

    rlat = rlat*2d0*pi
    call PetscPrintf(PETSC_COMM_WORLD, "reciprocal lattice:"//New_Line('A'), ierr)
    do i = 1, 3
      write (str, fmt='(3f24.12)') rlat(1:3, i)
      call PetscPrintf(PETSC_COMM_WORLD, trim(str)//New_Line('A'), ierr)
    end do
    deltak(1:3, 1) = rlat(1:3, 1)/real(max(1, nkx), dp)
    deltak(1:3, 2) = rlat(1:3, 2)/real(max(1, nky), dp)
    deltak(1:3, 3) = rlat(1:3, 3)/real(max(1, nkz), dp)

    if (.not. makeks) return

    call PetscPrintf(PETSC_COMM_WORLD, "steps reciprocal lattice:"//New_Line('A'), ierr)
    do i = 1, 3
      write (str, fmt='(3f24.12)') deltak(1:3, i)
      call PetscPrintf(PETSC_COMM_WORLD, trim(str)//New_Line('A'), ierr)
    end do

    l_withgamma = .true.
    if (ksym) then
      nkp = floor(nkx*nky*nkz*0.5d0) + 1
      if (mod(nkx*nky*nkz, 2) .eq. 0) then
        nkp = nkp - 1
        l_withgamma = .false.
      end if
    else
      nkp = nkx*nky*nkz
    end if

    dk = 0d0
    if (mod(nkx, 2) .eq. 0) dk = dk + deltak(1:3, 1)*0.5d0
    if (mod(nky, 2) .eq. 0) dk = dk + deltak(1:3, 2)*0.5d0
    if (mod(nkz, 2) .eq. 0) dk = dk + deltak(1:3, 3)*0.5d0

    allocate (kps(3, nkp), wkps(nkp), ikps(3, nkp), stat=ierr)
    if (ierr .ne. 0) then
      write (errormsg, *) "allocation error ", ierr
      call error()
    end if

    kps = 0d0
    wkps = 0d0
    n1 = 0d0

    in1 = -1
    ik1_loop: do ik1 = floor(nkx*0.5d0), -floor(nkx*0.5d0), -1
      if ((mod(nkx, 2) .eq. 0) .and. (ik1 .eq. 0)) cycle ik1_loop
      in1 = in1 + 1
      in2 = -1
      ik2_loop: do ik2 = floor(nky*0.5d0), -floor(nky*0.5d0), -1
        if ((mod(nky, 2) .eq. 0) .and. (ik2 .eq. 0)) cycle ik2_loop
        in2 = in2 + 1
        in3 = -1        
        ik3_loop: do ik3 = floor(nkz*0.5d0), -floor(nkz*0.5d0), -1
          if ((mod(nkz, 2) .eq. 0) .and. (ik3 .eq. 0)) cycle ik3_loop
  
          if (ksym) then
            do i1 = 1, n1
              if (all((/ik1, ik2, ik3/) + ikps(1:3, i1) .eq. 0)) then
                wkps(i1) = wkps(i1) + 1
                cycle ik3_loop
              end if
            end do
          end if
          
          in3 = in3 + 1
          n1 = n1 + 1
          kps(1:3, n1) = -dk(1:3) + (floor(nkx*0.5d0) - in1)*deltak(1:3, 1) +&
            (floor(nky*0.5d0) - in2)*deltak(1:3, 2) + (floor(nkz*0.5d0) - in3)*deltak(1:3, 3)
          ikps(1:3, n1) = (/ik1, ik2, ik3/)
          wkps(n1) = 1
        end do ik3_loop
      end do ik2_loop
    end do ik1_loop

! ----- this gives for even number of k-points non symmetric BZ sampling
! i think it would be better to shift dk/2 to get a symmetric sampling (see above)
!~       nk1=nkx-mod(nkx-1,2)
!~       nk2=nky-mod(nky-1,2)

!~       if (ksym) then
!~         nkp=nkx*nky-nk1*nk2+nk1*nk2/2+1
!~       else
!~         nkp=nkx*nky
!~       end if

!~       allocate(kps(2,nkp),wkps(nkp),ikps(2,nkp),stat=ierr)
!~       if (ierr.ne.0) then
!~         write(0,*) "allocation error ",ierr
!~         stop
!~       end if

!~       n1=0
!~       wkps=0d0

!~       ik2_loop: do ik2=-ceiling(real(nky-1)/2d0),(nky-1)/2
!~         ik1_loop: do ik1=-ceiling(real(nkx-1)/2d0),(nkx-1)/2

!~           if (ksym) then
!~             do i1=1,n1-1
!~               if (all((/ ik1, ik2 /)+ikps(1:2,i1).eq.0)) then
!~                 wkps(i1)=wkps(i1)+1
!~                 cycle ik1_loop
!~               end if
!~             end do
!~           end if

!~           n1=n1+1
!~           kps(1:2,n1)=ik1*rlat(1:2,1)+ik2*rlat(1:2,2)+rlat(1:2,1)*0.5d0*mod(nkx+1,2)+rlat(1:2,2)*0.5d0*mod(nky+1,2)
!~           ikps(1:2,n1)=(/ ik1, ik2 /)
!~           wkps(n1)=1

!~         end do ik1_loop
!~       end do ik2_loop
! ---------------------------

    call PetscPrintf(PETSC_COMM_WORLD, "k-point shift"//New_Line('A'), ierr)
    write (str, fmt='(3e24.12)') dk(1:3)
    call PetscPrintf(PETSC_COMM_WORLD, trim(str)//New_Line('A'), ierr)
    call PetscPrintf(PETSC_COMM_WORLD, "k-points and weights:"//New_Line('A'), ierr)
    do ik1 = 1, nkp
      write (str, fmt='(4i8,3e24.12,i4)') ik1, ikps(1:3, ik1), kps(1:3, ik1), wkps(ik1)
      call PetscPrintf(PETSC_COMM_WORLD, trim(str)//New_Line('A'), ierr)
    end do

  end subroutine get_rlat

end module lattice_mod
