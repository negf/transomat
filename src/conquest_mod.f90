module conquest_mod
  implicit none

contains

  subroutine read_conquest_info(data_dir, xyz, cell, nat, tomat, tomatr, inao, neigh, nlist, ndim, atoms)
    use MPI
    use kinds
    use error_handler

    use globals, only: imu_ecc, l_ionode, ecc_dir
    implicit none

    character(*) :: data_dir
    real(dp) :: xyz(:, :), tomatr(:, :)
    real(dp) :: cell(3, 3)
    integer, allocatable :: atoms(:)
    integer, allocatable :: tomat(:, :), inao(:), neigh(:), ndim(:), nlist(:), tmp_nlist(:)
    integer :: nat, nfiles
    integer :: i1, nneigh, iat1, iat2, idum, natfile, &
      ifile, iunit, ioff, ierr, imat, imu1, imu2, &
      inu1, inu2, iat3, intvec(3)
    integer, allocatable :: natms(:), nnao(:)
    integer, allocatable :: niks(:, :), miks(:, :)
    real(dp), allocatable :: riks(:,:)
    real(dp) :: vec(3), dist(3), vec2(3)
    
    if (l_ionode) then
      allocate (nlist(0))
      ioff = 0; inao = 0; neigh = 0; ndim = 0
      write(0,*) trim(data_dir)//"/mat_info.dat"
      open (newunit=iunit, file=trim(data_dir)//"/mat_info.dat", status="old")
      read (iunit, *) nfiles

      iat1 = 0; ioff = 0; imat = 0; nneigh = 0
      do ifile = 1, nfiles

        read (iunit, *) idum, natfile
        ioff = ioff + natfile

        read (iunit, *) inao(ioff - natfile + 1:ioff)
        read (iunit, *) neigh(ioff - natfile + 1:ioff)
        read (iunit, *) ndim(ioff - natfile + 1:ioff)

        nneigh = nneigh + sum(neigh)
        allocate (tmp_nlist(nneigh), stat=ierr)
        if (ierr .ne. 0) then
          write (errormsg, *) "allocation error ", ierr
          call error()
        end if
        tmp_nlist(1:size(nlist)) = nlist(1:size(nlist))
        call move_alloc(tmp_nlist, nlist)

        do i1 = 1, natfile
          iat1 = iat1 + 1
          read (iunit, *) iat2
          atoms(iat1) = iat2
          allocate (nnao(neigh(iat1)), natms(neigh(iat1)), niks(3, neigh(iat1)), &
         &  riks(3, neigh(iat1)),  miks(3, neigh(iat1)), stat=ierr)
          if (ierr .ne. 0) then
            write (errormsg, *) "allocation error ", ierr
            call error()
          end if
          read (iunit, *) nnao(1:neigh(iat1))
          read (iunit, *) natms(1:neigh(iat1))
          nlist(1 + sum(neigh(1:iat1 - 1)):sum(neigh(1:iat1))) = natms(1:neigh(iat1))
          do iat2 = 1, neigh(iat1)
            read (iunit, *) vec2(1:3), intvec(1:3)
! now  Conquest dumps in reduced coordinates but still relative distances
            vec = vec2(1)*cell(1:3, 1) + vec2(2)*cell(1:3, 2) + vec2(3)*cell(1:3, 3)
            dist = vec - xyz(1:3, atoms(iat1)) + xyz(1:3, natms(iat2))
            niks(1, iat2) = nint(dist(1)/cell(1, 1))
            niks(2, iat2) = nint(dist(2)/cell(2, 2))
            niks(3, iat2) = nint(dist(3)/cell(3, 3))
            miks(1:3, iat2) = intvec(1:3)            
            riks(1:3, iat2) = vec2(1:3)
          end do
!~            write(6,*) "neigh(iat1)",neigh(iat1)
          do iat2 = 1, neigh(iat1)
            do imu2 = 1, nnao(iat2)
              do imu1 = 1, inao(iat1)
                inu1 = 0
                do iat3 = 1, atoms(iat1) - 1
                  inu1 = inu1 + imu_ecc(iat3)
                end do
                inu1 = inu1 + imu1
                inu2 = 0
                do iat3 = 1, natms(iat2) - 1
                  inu2 = inu2 + imu_ecc(iat3)
                end do
                inu2 = inu2 + imu2
                imat = imat + 1
                tomat(1, imat) = inu1
                tomat(2, imat) = inu2
                tomat(3:5, imat) = niks(1:3, iat2)
                tomat(6:8, imat) = miks(1:3, iat2)
                tomatr(1:3, imat) = riks(1:3, iat2)
              end do
            end do
          end do

          deallocate (nnao, natms, niks, miks, riks)

        end do
      end do
      close (iunit)
    end if

    call MPI_Bcast(nneigh, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (.not. l_ionode) then
      allocate (nlist(nneigh), stat=ierr)
      if (ierr .ne. 0) then
        write (errormsg, *) "allocation error ", ierr
        call error()
      end if
    end if
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call MPI_Bcast(nlist, size(nlist), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(inao, size(inao), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(ndim, size(ndim), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(tomat, size(tomat), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  end subroutine read_conquest_info

end module conquest_mod
