module readsys
contains

  subroutine read_sysinfo(sys_xyz, sys_species, sys_nmu, natoms, dlat, ncell, nmat, &
    ef, ne, nspecies, species, sysfile)
    use petsc
    use kinds
    use lattice_mod, only: read_cell
    use globals, only: t_species, inode, nprocs
    use petsc_mod, only: pstr_out, petsc_print_master
    use error_handler
    use MPI

    implicit none

    real(dp), allocatable, intent(out) :: sys_xyz(:, :)
    real(dp) :: dlat(3, 3), ef, ne(2)
    integer :: ncell(3), nmat
    integer, allocatable, intent(out) :: sys_nmu(:)
    integer, intent(out) :: natoms, nspecies
    integer, allocatable, intent(out)  :: sys_species(:)
    type(t_species), allocatable :: species(:)

    character(*), intent(in) :: sysfile

    integer :: iunit, io, ierr, i, iat, ityp, j, nbf, l, l1, ip, i1, irow, ispecies
    character(256) ::  strin

    if (inode .eq. 0) then
      open (newunit=iunit, file=trim(sysfile), action="read", status="old")
      call read_cell(dlat, iunit)
      read (iunit, *) ncell(1:3)
      read (iunit, *) nmat
      read (iunit, *) ef
      read (iunit, *) ne(1)
      read (iunit, *) ne(2)
      read (iunit, *) natoms
      allocate (sys_xyz(3, natoms), sys_nmu(natoms), sys_species(natoms), stat=ierr)
      if (ierr .ne. 0) then
        write (errormsg, fmt='(A,i8)') "allocation error ", ierr
        call error()
      end if

      do i = 1, natoms
        read (iunit, *) iat ,sys_species(i), sys_nmu(i), sys_xyz(1:3, i)
      end do
      read(iunit, *) strin
      read(iunit, *) nspecies
      allocate(species(nspecies), stat = ierr)
      if (ierr .ne. 0) then
        write (errormsg, fmt='(A,i8)') "allocation error ", ierr
        call error()
      end if
      
      do i = 1, nspecies
        read(iunit, *) j, j, species(i)%nbf
        allocate(species(i)%l(species(i)%nbf), species(i)%ibf(species(i)%nbf), &
       &  stat = ierr)
        if (ierr .ne. 0) then
          write (errormsg, fmt='(A,i8)') "allocation error ", ierr
          call error()
        end if        
        nbf = 0
        do
          read(iunit, *) j, l
          do l1 = 1, 2*l + 1
            nbf = nbf + 1
            species(i)%l(nbf) = l
            species(i)%ibf(nbf) = j
          end do
          if (nbf .eq. species(i)%nbf) exit
        end do
        
      end do
      
      close (iunit)
    end if
    
    call MPI_Bcast(dlat, 9, MPI_DOUBLE_PRECISION, 0, PETSC_COMM_WORLD, ierr)
    call MPI_Bcast(ncell, 3, MPI_DOUBLE_PRECISION, 0, PETSC_COMM_WORLD, ierr)
    call MPI_Bcast(nmat, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
    call MPI_Bcast(ef,13, MPI_DOUBLE_PRECISION, 0, PETSC_COMM_WORLD, ierr)
    call MPI_Bcast(ne, 2, MPI_DOUBLE_PRECISION, 0, PETSC_COMM_WORLD, ierr)
    call MPI_Bcast(natoms, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
    call MPI_Bcast(nspecies, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
    if (inode .ne. 0) then
      allocate (sys_xyz(3, natoms), sys_nmu(natoms), sys_species(natoms), &
     &  species(nspecies), stat=ierr)
      if (ierr .ne. 0) then
        write (errormsg, fmt='(A,i8)') "allocation error ", ierr
        call error()
      end if       
    end if        
    call MPI_Bcast(sys_species, natoms, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
    call MPI_Bcast(sys_nmu, natoms, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
    call MPI_Bcast(sys_xyz, 3*natoms, MPI_DOUBLE_PRECISION, 0, PETSC_COMM_WORLD, ierr)
    call MPI_Bcast(species(:)%nbf, nspecies, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)    
    do i = 1, nspecies
      if (inode .ne. 0) then
        allocate(species(i)%l(species(i)%nbf),species(i)%ibf(species(i)%nbf), &
       &  stat = ierr)
        if (ierr .ne. 0) then
          write (errormsg, fmt='(A,i8)') "allocation error ", ierr
          call error()
        end if 
      end if
      call MPI_Bcast(species(i)%l(:), species(i)%nbf, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
      call MPI_Bcast(species(i)%ibf(:), species(i)%nbf, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
    end do
  
!~     do ip = 0, nprocs - 1 
!~       if (ip .eq. inode) then
!~         irow = 0    
!~         do iat = 1, natoms
!~           ispecies = sys_species(iat)
!~           do i1 = 1, species(ispecies)%nbf
!~             irow = irow + 1
!~             write(0, fmt='(6i8)') ip, iat, ispecies, irow, species(ispecies)%ibf(i1), &
!~            &  species(ispecies)%l(i1)
!~           end do
!~         end do
!~       end if
!~       call MPI_BARRIER(PETSC_COMM_WORLD, ierr)
!~     end do
      
!~    ncell(3) = 1
  end subroutine read_sysinfo

end module readsys
