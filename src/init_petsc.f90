module petsc_control

  implicit none

contains

  subroutine petsc_solver_init()
#include <petsc/finclude/petscmat.h>
    use petscmat
    use globals, only: mattype_surf, matsolvertype_surf, mattype_cc, matsolvertype_cc,&
   &mattype_dense, mattype_sparse, mattype_electrode
    implicit none
    integer :: ierr
    PetscErrorCode  :: p_ierr
    character(256) :: instr
    logical :: lflag
    PetscBool :: p_lflag

! default mattyps and solvers
    mattype_surf = MATELEMENTAL
!~     mattype_surf = MATSCALAPACK
    mattype_electrode = MATMPIAIJ
!~       mattype_electrode=MATDENSE
    mattype_cc = MATMPIAIJ
    mattype_dense = MATMPIDENSE
    mattype_sparse = MATMPIAIJ

    matsolvertype_surf = MATSOLVERELEMENTAL ! this must be a solver for dense matrices
!~     matsolvertype_surf = MATSOLVERSCALAPACK                     ! this must be a solver for dense matrices


    matsolvertype_cc="petsc_default"

    call PetscOptionsGetString(PETSC_NULL_OPTIONS,&
      PETSC_NULL_CHARACTER, '-pc_factor_mat_solver_type', instr, p_lflag, p_ierr)
    if (p_lflag) then
      matsolvertype_cc = trim(adjustl(instr))
    end if
    call PetscPrintf(PETSC_COMM_WORLD, "Matsolver_CC  : "//trim(matsolvertype_cc)//New_Line('A'), ierr)
    call PetscPrintf(PETSC_COMM_WORLD, "Matsolver_Surf: "//trim(matsolvertype_surf)//New_Line('A'), ierr)

  end subroutine petsc_solver_init

  subroutine petsc_init()
#include <slepc/finclude/slepceps.h>
    use slepceps
    use petscmat
    use globals, only: inode, nprocs, l_ionode, l_loadsavegf

    implicit none

    integer :: rank, np, iilen, ierr, i, ilow, iup, nup, nlow
    PetscErrorCode :: p_ierr
    character(128) :: my_name, outstr, option_file
    PetscInt:: major, minor, subminor, release

!~       call PetscInitialize(PETSC_NULL_CHARACTER,p_ierr)

!~     call SlepcInitialize(PETSC_NULL_CHARACTER, ierr)
    option_file="transp.petsc"
!~     call PetscInitialize(option_file, ierr)
    call SlepcInitialize(option_file, ierr)
!~     call PetscOptionsInsertFile(PETSC_COMM_WORLD, PETSC_NULL_OPTIONS, "./trans.petsc", PETSC_TRUE, ierr)


    call PetscPrintf(PETSC_COMM_WORLD, "Init PETSC \n", p_ierr)
    call PetscGetVersionNumber(major, minor, subminor, release, ierr)
    write (outstr, fmt='(4i8)'), major, minor, subminor, release
    call PetscPrintf(PETSC_COMM_WORLD, trim(outstr)//new_line('A'), p_ierr)

    call MPI_COMM_SIZE(PETSC_COMM_WORLD, np, ierr)
    nprocs = np
    call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)
    inode = rank

    l_ionode = 0 .eq. inode

  end subroutine petsc_init

  subroutine petsc_subcomm_init()
#include <slepc/finclude/slepceps.h>
    use slepceps
    use petscmat
    use globals, only: inode, nprocs, l_ionode, inode_sub, psubcomm, &
      inode_group, ngroups, nodes_group, group_range

    implicit none

    integer :: rank, np, iilen, ierr, i, ilow, iup, nup, nlow
    PetscErrorCode :: p_ierr
    character(256) :: my_name

    nup = ceiling(real(nprocs)/real(ngroups))
    nlow = nup - 1
    iup = nprocs - ngroups*(nup - 1)
    ilow = ngroups - iup

    allocate (nodes_group(ngroups))

    if ((inode + 1) .le. (iup*nup)) then
      inode_group = ceiling(real(inode + 1)/real(nup))
    else
      inode_group = iup + ceiling(real(inode + 1 - nup*iup)/real(nlow))
    end if

    if (inode_group .le. iup) then
      group_range(1) = (inode_group - 1) * nup
      group_range(2) = inode_group * nup - 1
    else
      group_range(1) = iup * nup
      group_range(2) = iup * nup  + (inode_group - iup) * nlow - 1
    end if


    do i = 1, ngroups
      if (i .le. iup) then
        nodes_group(i) = nup
      else
        nodes_group(i) = nlow
      end if
    end do

    call MPI_Comm_split(PETSC_COMM_WORLD, inode_group, inode, psubcomm, ierr);
    call MPI_Comm_rank(psubcomm, inode_sub, ierr)

    call MPI_Get_processor_name(my_name, iilen, ierr)
    if (l_ionode) write (0, fmt='(11X,A,7X,A,5X,A,5X,A,X,A,X,A)') "inode", "inode_sub", "inode_group", "nodes_group", "group range", "host name"
    call flush (0)
    call MPI_BARRIER(PETSC_COMM_WORLD, ierr)
    do i = 0, nprocs - 1
      if (i .eq. inode) write (0, fmt='(4i16,X,i12,i12,XA)') inode, inode_sub, inode_group, nodes_group(inode_group), group_range, trim(adjustl(my_name))
      call flush (0)
      call MPI_BARRIER(PETSC_COMM_WORLD, ierr)
    end do

  end subroutine petsc_subcomm_init

!~     subroutine petsc_set_locals()

!~       use globals, only : inode,nprocs,nmu_ecc,nmu_elec_l,nmu_elec_r,nmu_c,nmu_l,nmu_r

!~       implicit none

!~ #include "finclude/petsc.h90"

!~       integer :: np,ip
!~       character(256) :: numstr,hello
!~       PetscErrorCode :: ierr

!~       if (inode.eq.0) write(6,*) "number of nodes",nprocs

!~       do ip=0,nprocs-1

!~         if (ip.eq.inode) then
!~           write(numstr,*) inode
!~           numstr=adjustl(numstr)
!~           hello = "Init node "//trim(numstr)
!~           write(numstr,*) nprocs-1
!~           numstr=adjustl(numstr)
!~           hello = trim(hello)// " / "//numstr//" ! "
!~           call PetscPrintf( PETSC_COMM_SELF, trim(hello)//achar(10), ierr )

!~           nmu_ecc_loc=petsc_getlocal_nmu(nmu_ecc,nprocs,inode)
!~           nmu_elec_l_loc=petsc_getlocal_nmu(nmu_elec_l,nprocs,inode)
!~           nmu_elec_r_loc=petsc_getlocal_nmu(nmu_elec_r,nprocs,inode)
!~           nmu_c_loc=petsc_getlocal_nmu(nmu_c,nprocs,inode)
!~           nmu_l_loc=petsc_getlocal_nmu(nmu_l,nprocs,inode)
!~           nmu_r_loc=petsc_getlocal_nmu(nmu_r,nprocs,inode)

!~           write(6,*) "nmu_ecc_loc",nmu_ecc,nmu_ecc_loc
!~           write(6,*) "nmu_elec_l_loc",nmu_elec_l,nmu_elec_l_loc
!~           write(6,*) "nmu_elec_r_loc",nmu_elec_r,nmu_elec_r_loc
!~           write(6,*) "nmu_c_loc",nmu_c,nmu_c_loc
!~           write(6,*) "nmu_l_loc",nmu_l,nmu_l_loc
!~           write(6,*) "nmu_r_loc",nmu_r,nmu_r_loc

!~         end if
!~         call mpi_barrier(MPI_COMM_WORLD,ierr)
!~       end do

!~     end subroutine petsc_set_locals

  function petsc_getlocal_nmu(n, np, inode)
    implicit none

    integer :: petsc_getlocal_nmu, n, np, inode

    integer :: nl

    nl = ceiling(real(n)/np)
    if (nl*(inode + 1) .gt. n) then
      nl = n - nl*inode
    end if

    petsc_getlocal_nmu = nl

  end function

  subroutine set_up_petsc_mat(n1l, n1u, n2l, n2u, infofile)
    use globals
    implicit none

    character(*) :: infofile
    integer :: n1l, n1u, n2l, n2u
    PetscInt, allocatable :: d_nnz(:), o_nnz(:)

  end subroutine set_up_petsc_mat

end module petsc_control
