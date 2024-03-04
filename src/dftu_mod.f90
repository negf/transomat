module dftu_mod
#include <petsc/finclude/petsc.h>
  use petsc  
  
  implicit none
  
  character(256) :: dftu_projector_info(6)=[ "P.S+S.P", "S1/2*P*S1/2", "S*P*S" , &
 &                                           "S*P+P*S", "Tr_sym(\hat{D})", "Tr{PDP'}" ]
  
  type t_Uoperator
    logical :: l_apply_U = .false.
    Mat :: Pk
    PetscScalar :: U = 0d0   
  end type t_Uoperator
  
  integer :: n_Ushells
  type(t_Uoperator), allocatable :: Uoperator(:)
  
  contains
  
    subroutine init_Pkmat(p_Pk_pre)

#include <petsc/finclude/petsc.h>
      use petsc  
      use petsc_mod
      use globals
      use k_on_demand_mod
      
      implicit none    
      
      Mat :: p_Pk_pre
      
      integer :: nrow, ncol, nl1, nl2, nlc1, nlc2, ierr
      
      call MatGetSize(p_k00_ik_cc(1), nrow, ncol, ierr)
      call MatGetOwnershipRange(p_k00_ik_cc(1), nl1, nl2, ierr)
      call MatGetOwnershipRangeColumn(p_k00_ik_cc(1), nlc1, nlc2, ierr)
      
      call MatCreate(PETSC_COMM_WORLD, p_Pk_pre, ierr)
      call MatSetSizes(p_Pk_pre, nl2 - nl1, nlc2 - nlc1, nrow, ncol, ierr)
      call MatSetType(p_Pk_pre, MATPREALLOCATOR, ierr)
      call MatSetUp(p_Pk_pre, ierr)
    
    end subroutine init_Pkmat
  
    subroutine init_Uoperator()

#include <petsc/finclude/petsc.h>
      use petsc  
      use petsc_mod
      use globals
      use k_on_demand_mod    
      
      implicit none  
      
      Mat :: p_Pk_pre
      integer :: iat, ishell1, ispecies, ierr, i1, maxnbf, nz, i2, ishell_at, &
     &  irow, nl1, nl2, nzdiag, jrow, maxl
      integer, allocatable :: cols(:), nzloc(:)
      character(256) :: sshell, ss
      logical :: l_Ufound
      real(dp) :: U_ishell

      call MatGetOwnershipRange(p_k00_ik_cc(1), nl1, nl2, ierr)
      nl2 = nl2 - 1

      n_Ushells = 0
      maxnbf = 0
      maxl = 0
      do iat = 1, nat_ecc
        ispecies = species_ecc(iat)        
        n_Ushells = n_Ushells + species_c(ispecies)%nUshells_at
        maxl = max(maxl, maxval(species_c(ispecies)%l(:)))
        maxnbf = max(maxl, species_c(ispecies)%nbf)
      end do
      maxl = 2 * maxl + 1
      allocate(Uoperator(n_Ushells), cols(maxl * maxnbf), nzloc(nl2-nl1+1))
      
!~       if (inode.eq.0) write(0,fmt='(12A11)') "inode", "iat", "ispecies", "ishell_at", "i1", "ishell1", "Ushell(i1)", &
!~      &  "ishell_at", "nl1", "nl2", "jrow", "irow", "nzdiag"
      call MPI_BARRIER(PETSC_COMM_WORLD, ierr)
        irow = -1        
        nzloc = 0
        jrow = -1
        l_Ufound = .false.
        ishell1 = 0
        
        do iat = 1, nat_ecc
        
          nzdiag = 0
          cols = -1
          ispecies = species_ecc(iat)
          
          do ishell_at = 1, maxval(species_c(ispecies)%U_shell(:))
          
            irow = jrow
            do i1 = 1, species_c(ispecies)%nbf
              irow = irow + 1
              if ((species_c(ispecies)%l_U(i1)) .and. (species_c(ispecies)%U_shell(i1) .eq. ishell_at)) then            
                l_Ufound = .true.
                U_ishell = species_c(ispecies)%U(i1)
                if ((irow .ge. nl1) .and. (irow .le. nl2)) then
                  nzdiag = nzdiag + 1
                  cols(nzdiag) = irow
                  nzloc(irow - nl1 + 1) = 1
                end if
              end if         
!~               write(0,fmt='(13i11,X,3l)') inode, iat, ispecies,ishell_at, i1, ishell1, species_c(ispecies)%U_shell(i1), ishell_at, nl1, nl2, jrow, irow, &
!~              &  nzdiag, species_c(ispecies)%l_U(i1), species_c(ispecies)%U_shell(i1).eq.ishell_at, l_Ufound    
            end do 

          if (l_Ufound) then
          
            ishell1 = ishell1 + 1
            Uoperator(ishell1)%l_apply_U = .true.
            Uoperator(ishell1)%U = U_ishell
            
            call init_Pkmat(p_Pk_pre)
            do i2 = 1, nzdiag
              call MatSetValue(p_Pk_pre, cols(i2), cols(i2), p_zero, INSERT_VALUES, ierr)   
            end do
            call MatAssemblyBegin(p_Pk_pre, MAT_FINAL_ASSEMBLY, ierr)
            call MatAssemblyEnd(p_Pk_pre, MAT_FINAL_ASSEMBLY, ierr)            
            
            call petsc_get_a_with_b_c(Uoperator(ishell1)%Pk, p_Pk_pre, p_Pk_pre, mattype_sparse)
            call MatPreallocatorPreallocate(p_Pk_pre, PETSC_TRUE, Uoperator(ishell1)%Pk, ierr)
            call MatDestroy(p_Pk_pre, ierr)
              
            do i2 = 1, nzdiag
              call MatSetValue(Uoperator(ishell1)%Pk, cols(i2), cols(i2), p_one, INSERT_VALUES, ierr)   
            end do            
            call MatAssemblyBegin(Uoperator(ishell1)%Pk, MAT_FINAL_ASSEMBLY, ierr)
            call MatAssemblyEnd(Uoperator(ishell1)%Pk, MAT_FINAL_ASSEMBLY, ierr)
            
!~             write(sshell, fmt='(i16)')  ishell1
!~             write(ss, fmt='(e24.12)') real(Uoperator(ishell1)%U, 8)
!~             sshell = "Pk_"//trim(adjustl(sshell))//"_"//trim(adjustl(ss))
!~             call petsc_mat_info(Uoperator(ishell1)%Pk,trim(sshell) ,ierr)  
!~             call dump_nonzero_structure(Uoperator(ishell1)%Pk, trim(sshell), PETSC_COMM_WORLD, 0)

            l_Ufound = .false.
            nzloc = 0
            cols = -1
            
          end if
          
        end do
        jrow = jrow + species_c(ispecies)%nbf                                        
      end do
      
    end subroutine init_Uoperator
    
    subroutine test_ne2(p_mat, nek)

#include <petsc/finclude/petsc.h>
    use petsc  
    use petsc_mod
    use globals
    use k_on_demand_mod    
    implicit none      
    
    Mat :: p_mat(1)
    Mat :: p_SD,p_DS,p_S,p_D
    integer :: ik, ierr
    PetscScalar :: ne, nek
    
    
    
    p_S=p_s00_ik_cc(1)
    p_D=p_mat(1)
!~     call dump_nonzero_structure(p_D, strin, PETSC_COMM_WORLD, 0, 1)
    call MatMatMult(p_D, p_S, MAT_INITIAL_MATRIX, &
     &  PETSC_DEFAULT_REAL, p_SD, ierr)
    call MatGetTrace(p_SD, nek, ierr)
!~     write(pstr_out, fmt = '(A,i8,2e24.12,i8)') "Nek ",iik, nek, wkp_c(iik)
!~     call petsc_print_master()
    nek = nek * wkp_c(iik)
    call MatDestroy(p_SD, ierr)
  
    end subroutine test_ne2
    
    subroutine test_ne()
#include <petsc/finclude/petsc.h>
    use petsc  
    use petsc_mod
    use globals
    use k_on_demand_mod    
    implicit none      
    
    Mat :: p_mat(1)
    Mat :: p_SD,p_DS,p_S,p_D
    integer :: ik, ierr
    PetscScalar :: ne, nek
    
    !~     call MatTranspose(p_k00_ik_cc(1),  MAT_INITIAL_MATRIX, p_D, ierr)
    
    p_S=p_s00_ik_cc(1)
    p_D=p_k00_ik_cc(1)
    
    
    call k_on_demand(1, 2)
    call MatMatMult(p_S, p_D, MAT_INITIAL_MATRIX, &
     &  PETSC_DEFAULT_REAL, p_SD, ierr)
    call MatMatMult(p_D, p_S, MAT_INITIAL_MATRIX, &
     &  PETSC_DEFAULT_REAL, p_DS, ierr)
     
     
      ne = 0d0
      do ik = nk_c, 1, -1
      
        iik  = ik
        call k_on_demand(ik, 2)

        call MatMatMult(p_S, p_D, MAT_REUSE_MATRIX, &
       &  PETSC_DEFAULT_REAL, p_SD, ierr)
        call MatMatMult(p_D, p_S, MAT_REUSE_MATRIX, &
       &  PETSC_DEFAULT_REAL, p_DS, ierr)
        call MatAXPY(p_SD, p_one, p_DS,  SAME_NONZERO_PATTERN, ierr)
        call MatGetTrace(p_SD, nek, ierr)
        nek = 0.5d0*nek
        write(pstr_out, fmt = '(A,i8,2e24.12,i8)') "nek ",iik, nek, wkp_c(iik)
        call petsc_print_master()
        nek = nek * wkp_c(iik)
        ne = ne + nek
        
  !~       call MatDestroy(p_DS, ierr)
  !~       call MatDestroy(p_S, ierr)
      end do    
      call MatDestroy(p_SD, ierr)
      call MatDestroy(p_DS, ierr)
      ne = ne / (nktot)
      write(pstr_out, fmt = '(2e24.12)') ne
      call petsc_print_master()
    end subroutine test_ne
  
    subroutine read_dftu(dftu_file, species)
      use petsc
      use globals, only: inode, t_species, l_ionode, eh, nprocs, nat_ecc, species_ecc, &
     &  dftu_projector
      use kinds
      use petsc_mod, only: pstr_out, petsc_print_master
      use MPI
      
      implicit none
      
      character(*) :: dftu_file
      type(t_species) :: species(:)
      integer :: ispecies, ierr, iunit, ibf, nspecies, i, j, ip, i1, irow, iat, n, &
     &  ishell, jshell
      real(dp) :: U
      logical :: lexist
      
      pstr_out="DFT+U projector "//trim(dftu_projector_info(dftu_projector))
      call petsc_print_master()
      
      nspecies=size(species)
      
      do i = 1, nspecies
        allocate(species(i)%l_U(species(i)%nbf), species(i)%U(species(i)%nbf), &
       &   species(i)%U_shell(species(i)%nbf), stat = ierr)
        species(i)%l_U = .false.
        species(i)%U = 0d0
        species(i)%U_shell = 0
        species(i)%nUshells_at = 0
      end do
      
      
      if (l_ionode) then
        inquire(file = trim(dftu_file), exist = lexist)
        if (lexist) then
          open(newunit = iunit, file = trim(dftu_file), action = "read", status = "old")
          ishell = 0
          do 
            read(unit = iunit, fmt = *, iostat = ierr) n, U
            if (ierr .ne. 0) exit              
            ishell = ishell + 1
            do j = 1, n            
              read(unit = iunit, fmt = *, iostat = ierr) ispecies, ibf
              if (ierr .ne. 0) exit            
              if (j.eq.1) species(ispecies)%nUshells_at = species(ispecies)%nUshells_at + 1
              do i = 1, species(ispecies)%nbf
                if (species(ispecies)%ibf(i) .eq. ibf) then
                  species(ispecies)%l_U(i) = .true.
                  species(ispecies)%U_shell(i) = ishell
                  species(ispecies)%U(i) = U / eh         
                end if
              end do  
            end do
            if (ierr .ne. 0) exit
          end do
        end if
        
      end if
      
      
      pstr_out = "dft+U: "//trim(dftu_file)
      call petsc_print_master()
      do i = 1, nspecies
        call MPI_Bcast(species(i)%l_U(:), species(i)%nbf, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr)
        call MPI_Bcast(species(i)%nUshells_at, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
        call MPI_Bcast(species(i)%U_shell(:), species(i)%nbf, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
        call MPI_Bcast(species(i)%U(:), species(i)%nbf, MPI_DOUBLE_PRECISION, 0, & 
       &  PETSC_COMM_WORLD, ierr)
        do j = 1, species(i)%nbf
          write(pstr_out, fmt = '(4i4,e24.12,X,l)') i, j, species(i)%l(j), species(i)%U_shell(j), &
         & species(i)%U(j), species(i)%l_U(j)
          call petsc_print_master()
        end do
      end do 
!~     call MPI_Bcast(lexist, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr)  
!~     if (.not.lexist) return  
!~     do ip = 0, nprocs - 1 
!~       if (ip .eq. inode) then
!~         irow = 0    
!~         do iat = 1, nat_ecc
!~           ispecies = species_ecc(iat)
!~           do i1 = 1, species(ispecies)%nbf
!~             irow = irow + 1
!~             write(0, fmt='(i8,A,5i8, l)') ip, " test ", iat, ispecies, irow, &
!~            &  species(ispecies)%ibf(i1), species(ispecies)%l(i1), species(ispecies)%l_U(i1)
!~           end do
!~         end do
!~       end if
!~       call MPI_BARRIER(PETSC_COMM_WORLD, ierr)
!~     end do      
      
    end subroutine read_dftu
    
  
    subroutine add_dftu_to_H(l_save_to_disk)
#include <petsc/finclude/petsc.h>
      use petsc  
      use petsc_mod
      use k_on_demand_mod
      use ft_mod
      use globals
      use error_handler
      use misc
      implicit none  
      
      logical :: l_save_to_disk
      
      Mat, allocatable, target :: p_H_tmp(:, :, :), p_mat_tmp(:)
      Mat :: A
      PetscScalar :: fac
      PetscReal :: norm01, norm02, norm03
      integer :: ik, ierr, i1 ,i2 ,i3
      character(256) :: hfile
      
      
      l_dftu_sparse_projector = .true.
      
      allocate(p_H_tmp(-ncell_c(1):ncell_c(1), -ncell_c(2):ncell_c(2), -ncell_c(3):ncell_c(3)),&
     &  p_mat_tmp(1), stat = ierr)
      
      
      do i1 = -ncell_c(1), ncell_c(1)
        do i2 = -ncell_c(2), ncell_c(2)
          do i3 = -ncell_c(3), ncell_c(3)
            write(pstr_out, fmt = '(3i8)') i1, i2, i3
            call petsc_print_master()            
            call MatDuplicate(p_h00_diag(i1, i2 ,i3)%p, MAT_SHARE_NONZERO_PATTERN, &
           &  p_H_tmp(i1, i2, i3), ierr)
            call MatZeroEntries(p_H_tmp(i1, i2 ,i3), ierr)
          end do
        end do
      end do
      
      do ik = nk_c, 1, -1
        iik = ik
        call k_on_demand(ik, 3)
        call MatNorm(p_h00k_cc(ik), NORM_FROBENIUS, norm01, ierr)
        
        if (dftu_projector .eq. 1) then
          call get_VU_SpD_DpS_v2(p_h00k_cc(ik), p_s00k_cc(ik), p_k00k_cc(ik))
        else if (dftu_projector .eq. 2) then
          call get_VU_TDT(p_h00k_cc(ik), p_s00k_cc(ik), p_k00k_cc(ik), 2)
        else if (dftu_projector .eq. 3) then
          call get_VU_TDT(p_h00k_cc(ik), p_s00k_cc(ik), p_k00k_cc(ik), 1)
        else if (dftu_projector .eq. 4) then
          call get_VU_SD_DS(p_h00k_cc(ik), p_s00k_cc(ik), p_k00k_cc(ik))
        else if (dftu_projector .eq. 5) then
          call get_VU_Tr_symD_v2(p_h00k_cc(ik), p_s00k_cc(ik), p_k00k_cc(ik))
        else if (dftu_projector .eq. 6) then
          call get_VU_Tr_PkDPk_bar(p_h00k_cc(ik), p_s00k_cc(ik), p_k00k_cc(ik))
        else 
          write (errormsg, fmt='(A,i8)') "unknown DFTU projector ", dftu_projector
          call error()
        end if
        call MatNorm(p_h00k_cc(ik), NORM_FROBENIUS, norm02, ierr)
        call MatHermitianTranspose(p_h00k_cc(ik), MAT_INITIAL_MATRIX, A, ierr)
        call MatAXPY(A,p_minus1, p_h00k_cc(ik),SAME_NONZERO_PATTERN, ierr)
        call MatNorm(A, NORM_FROBENIUS, norm03, ierr)        
        call MatDestroy(A,ierr)
        write (pstr_out, fmt='(A,i8, 3e18.8)') "Norm(H(k)-H(k)')", ik, norm01, norm02, norm03
        call petsc_print_master()
        p_mat_tmp(1) = p_h00k_cc(ik)
        call fourier3d_back_add(p_mat_tmp, p_H_tmp, kp_c(1:3, ik:ik), &
          wkp_c(ik:ik), 1, dlat_c, ncell_c(1), ncell_c(2), ncell_c(3), 0, 0, PETSC_FALSE)
      
      end do
      
      fac = 1d0 / real(nktot, 8)
      do i1 = -ncell_c(1), ncell_c(1)
        do i2 = -ncell_c(2), ncell_c(2)
          do i3 = -ncell_c(3), ncell_c(3)
            call MatScale(p_H_tmp(i1, i2, i3), fac, ierr)
            
            if (l_save_to_disk) then
              call MatNorm(p_H_tmp(i1, i2, i3), NORM_FROBENIUS, norm01, ierr)
              call MatNorm(p_h00_diag(i1, i2 ,i3)%p, NORM_FROBENIUS, norm02, ierr)
              call MatDuplicate(p_h00_diag(i1, i2 ,i3)%p, MAT_COPY_VALUES, &
            &  A, ierr)            
              call MatAXPY(A, p_minus1, p_H_tmp(i1, i2, i3),SAME_NONZERO_PATTERN, ierr)
              call MatNorm(A, NORM_FROBENIUS, norm03, ierr)                   
              write (pstr_out, fmt='(A,3i6, 3e18.8)') "test ",i1,i2,i3,norm01,norm02,norm03
              call petsc_print_master()
              hfile = "H_"//trim(spin_str)//"_"//trim(int2str(i1))//"_"//trim(int2str(i2))//&
            &  "_"//trim(int2str(i3))//"_petsc.dat"
              call petsc_mat_direct_save(p_H_tmp(i1 ,i2, i3), hfile, ierr)
              call MatDestroy(A,ierr)
            else
              call MatCopy(p_H_tmp(i1, i2, i3), p_h00_diag(i1, i2 ,i3)%p, &
             &  SAME_NONZERO_PATTERN, ierr)
              
            end if
            
            call MatDestroy(p_H_tmp(i1, i2, i3), ierr)
            
          end do
        end do
      end do      
    
    end subroutine add_dftu_to_H
  

  
  subroutine get_VU(p_H, p_S, p_D)
#include <petsc/finclude/petsc.h>
    use petsc  
    use petsc_mod
    use globals
    use k_on_demand_mod
    use error_handler
    implicit none  
    
    Mat :: p_H, p_S, p_D
    
    PetscScalar, allocatable :: S_v(:), D_v(:), HU_v(:)
    PetscScalar :: qij, Vij_1, Vij_2
    integer, allocatable :: cols_S(:)
    integer :: ierr, ncols, nrow_global, ncol_global, iat, i1, &
   &  nrow1, nrow2, ispecies, irow, j1, jspecies, jcol, jat, jpos, ij
    
    
    call MatGetOwnershipRange(p_S, nrow1, nrow2, ierr)
    nrow1 = nrow1 + 1 ! 1 indexed
    allocate(S_v(nzmax_row), D_v(nzmax_row), HU_v(nzmax_row), cols_S(nzmax_row), &
   &  stat = ierr)
    if (ierr .ne. 0) then
      write (errormsg, fmt='(A,i8)') "allocation error ", ierr
      call error()
    end if 
    
!~     call MatTranspose(p_S, MAT_INPLACE_MATRIX, p_S, ierr)
!~    qij = 0.5 ( Dij x Sji + Sij x Dji) = 0.5 ( Dij x Sij* + Sij x Dij*)

    irow = 0    
    do iat = 1, nat_ecc
      ispecies = species_ecc(iat)
      do i1 = 1, species_c(ispecies)%nbf
        irow = irow + 1
        if ((species_c(ispecies)%l_U(i1)) .and. (irow .ge. nrow1) .and. &
       &  (irow.le.nrow2)) then
       
          jcol = 0
          ncols = 0
          do jat = 1, nat_ecc
            jspecies = species_ecc(jat)
            do j1 = 1, species_c(jspecies)%nbf
              jcol = jcol + 1
              if ((species_c(ispecies)%U_shell(i1) .ne. species_c(jspecies)%U_shell(j1)) &
             &  .or. (iat .ne. jat))  cycle
              ncols = ncols + 1
              cols_S(ncols) = jcol - 1
            end do !ji
          end do !jat
          
          if (ncols .eq. 0) cycle
          
          call MatGetValues(p_S, 1, [irow -1], ncols, cols_S(1:ncols), S_v(1:ncols), ierr)
          call MatGetValues(p_D, 1, [irow -1], ncols, cols_S(1:ncols), D_v(1:ncols), ierr)
          
          do ij = 1, ncols
            qij  = 0.5d0 * (D_v(ij) * dconjg(S_v(ij)) + S_v(ij)* dconjg(D_v(ij)))
            Vij_1 = -(dconjg(S_v(ij) * qij) + dconjg(S_v(ij)) * qij) 
            Vij_2 = 0d0
            if ((irow - 1) .eq. cols_S(ij)) then
              Vij_2 = dconjg(S_v(ij))
            end if
             HU_v(ij) = 0.5d0 * species_c(ispecies)%U(i1) * (Vij_1 + Vij_2)
             write(0, fmt = '(2i4,A,3i6,9e14.6)') inode, iik," XX ",iat, irow, cols_S(ij) + 1, species_c(ispecies)%U(i1),qij , Vij_1, Vij_2, HU_v(ij) !,  species_c(ispecies)%l_U(i1)
          end do  !ij          
          call MatSetValues(p_H, 1, [irow -1], ncols, cols_S(1:ncols), HU_v(1:ncols), &
         &  ADD_VALUES, ierr)
        end if         
      end do !i1
      call MatAssemblyBegin(p_H, MAT_FINAL_ASSEMBLY, ierr)
      call MatAssemblyEnd(p_H, MAT_FINAL_ASSEMBLY, ierr)
    end do !iat

  end subroutine get_VU
  
  
  subroutine get_VU_SpD_DpS(p_H, p_S, p_D)
#include <petsc/finclude/petsc.h>
    use petsc  
    use petsc_mod
    use globals
    use k_on_demand_mod
    use error_handler
    implicit none  
    
    Mat :: p_H, p_S, p_D
    
    Mat :: p_PkQPk, p_PkSPk, p_Q
    Vec :: v_Diag
    PetscScalar :: p_scal   
    integer :: ierr, ishell1
    PetscInt, pointer :: colmapS(:), colmapD(:), colmapQ(:)
    character(256) :: ss
    
    
    call MatDuplicate(p_H,  MAT_SHARE_NONZERO_PATTERN, p_Q, ierr)
    call MatZeroEntries(p_Q, ierr)
  
    p_scal = 0.5d0
    call petsc_A_px_B(p_D, p_S, p_Q, p_scal, .false., .true., INSERT_VALUES)    
    call petsc_A_px_B(p_S, p_D, p_Q, p_scal, .false., .true., ADD_VALUES)
    
!~     p_scal = 2d0
!~     call MatScale(p_Q, p_scal, ierr)
    
    call MatCreateVecs(p_H, v_Diag, PETSC_NULL_VEC, ierr)


     
    do ishell1 = 1, n_Ushells
    
      
    
      call MatMatMatMult(Uoperator(ishell1)%Pk, p_S, Uoperator(ishell1)%Pk, &
     &  MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, p_PkSPk, ierr)
     
     
      call MatMatMatMult(Uoperator(ishell1)%Pk, p_Q, Uoperator(ishell1)%Pk, &
     &  MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, p_PkQPk, ierr)
      
      call MatGetDiagonal(p_PkSPk, v_Diag, ierr)      
      
      p_scal = 0.5d0 * Uoperator(ishell1)%U
      call VecScale(v_Diag, p_scal , ierr)
      call MatDiagonalSet(p_H, v_Diag, ADD_VALUES, ierr)
      
      p_scal = -0.5d0 * Uoperator(ishell1)%U
      call petsc_A_px_B(p_PkQPk, p_PkSPk, p_H, p_scal, .true., .true., ADD_VALUES)    
      call petsc_A_px_B(p_PkQPk, p_PkSPk, p_H, p_scal, .false., .true., ADD_VALUES)    
         
      
      call MatDestroy(p_PkSPk, ierr)
      call MatDestroy(p_PkQPk, ierr)
    end do
    
    call VecDestroy(v_Diag, ierr)
    call MatDestroy(p_Q, ierr)

  end subroutine get_VU_SpD_DpS
  
  subroutine get_VU_SpD_DpS_v2(p_H, p_S, p_D)
#include <petsc/finclude/petsc.h>
    use petsc  
    use petsc_mod
    use globals
    use k_on_demand_mod
    use error_handler
    implicit none  
    
    Mat :: p_H, p_S, p_D
    
    Mat :: p_PkQPk, p_PkSPk, p_Q
    Vec :: v_Diag
    PetscScalar :: p_scal   
    integer :: ierr, ishell1
    PetscInt, pointer :: colmapS(:), colmapD(:), colmapQ(:)
    character(256) :: ss
    
    
    call MatDuplicate(p_H,  MAT_SHARE_NONZERO_PATTERN, p_Q, ierr)
    call MatZeroEntries(p_Q, ierr)
  
    p_scal = 0.5d0
    call petsc_A_px_B(p_D, p_S, p_Q, p_scal, .false., .true., INSERT_VALUES)    
    call petsc_A_px_B(p_S, p_D, p_Q, p_scal, .false., .true., ADD_VALUES)
    
!~     p_scal = 2d0
!~     call MatScale(p_Q, p_scal, ierr)
    
    call MatCreateVecs(p_H, v_Diag, PETSC_NULL_VEC, ierr)


     
    do ishell1 = 1, n_Ushells
    
      
    
      call MatMatMatMult(Uoperator(ishell1)%Pk, p_S, Uoperator(ishell1)%Pk, &
     &  MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, p_PkSPk, ierr)
     
      
      call MatGetDiagonal(p_PkSPk, v_Diag, ierr)      
      
      p_scal = 0.5d0 * Uoperator(ishell1)%U
      call VecScale(v_Diag, p_scal , ierr)
      call MatDiagonalSet(p_H, v_Diag, ADD_VALUES, ierr)
      
      call MatConjugate(p_PkSPk, ierr)
      call MatMatMatMult(Uoperator(ishell1)%Pk, p_Q, Uoperator(ishell1)%Pk, &
     &  MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, p_PkQPk, ierr)
     
      
      p_scal = -Uoperator(ishell1)%U      
      call petsc_A_px_B(p_PkQPk, p_PkSPk, p_H, p_scal, .false., .true., ADD_VALUES)    
         
      
      call MatDestroy(p_PkSPk, ierr)
      call MatDestroy(p_PkQPk, ierr)
    end do
    
    call VecDestroy(v_Diag, ierr)
    call MatDestroy(p_Q, ierr)

  end subroutine get_VU_SpD_DpS_v2
  
  subroutine get_VU_TDT(p_H, p_S, p_D, Top)
#include <petsc/finclude/petsc.h>
    use petsc  
    use petsc_mod
    use petsc_wrapper
    use globals
    use k_on_demand_mod
    use error_handler
    implicit none  
    
    Mat, target :: p_H, p_S, p_D
    integer :: top
    
    Mat, target :: p_PkQPk, p_TPkT, p_TPkQPkT, p_Q, p_S12
    Mat, pointer :: p_T    
    Vec :: v_Diag
    PetscScalar :: p_scal   
    integer :: ierr, ishell1
    PetscInt, pointer :: colmapS(:), colmapD(:), colmapQ(:)
    character(256) :: ss
    
    if (Top .eq. 1) then
      p_T => p_S
    else if (Top .eq. 2) then
      call petsc_get_densemat(p_S, p_S12, mattype_dense)
      p_T => p_S12
      call petsc_get_sqrt_mat(p_S, PETSC_NULL_MAT, p_S12)
    end if
    if (.not.l_dftu_sparse_projector) then
      call MatConvert(p_H, mattype_dense, MAT_INPLACE_MATRIX, p_H, ierr)
    end if
    
    call MatMatMatMult(p_T, p_D, p_T, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, &
   &  p_Q, ierr)
        
!~     p_scal = 2d0
!~     call MatScale(p_Q, p_scal, ierr)

     
    do ishell1 = 1, n_Ushells
    
      
    
      call MatMatMatMult(p_T, Uoperator(ishell1)%Pk, p_T, &
     &  MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, p_TPkT, ierr)
     
     
      call MatMatMatMult(Uoperator(ishell1)%Pk, p_Q, Uoperator(ishell1)%Pk, &
     &  MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, p_PkQPk, ierr)
     
      call MatMatMatMult(p_T, p_PkQPk, p_T, MAT_INITIAL_MATRIX, &
     &  PETSC_DEFAULT_REAL, p_TPkQPkT, ierr)
      
      call MatConjugate(p_TPkT, ierr)               
      p_scal = 0.5d0 * Uoperator(ishell1)%U
      if (l_dftu_sparse_projector) then
        call petsc_aXpY(p_H, p_TPkT, p_scal, PETSC_FALSE, .false.)
      else
        call MatAXPY(p_H, p_scal, p_TPkT, SUBSET_NONZERO_PATTERN, ierr)
      end if

      p_scal = -Uoperator(ishell1)%U
      call MatConjugate(p_TPkQPkT, ierr)
      if (l_dftu_sparse_projector) then
        call petsc_aXpY(p_H, p_TPkQPkT, p_scal, PETSC_FALSE, .false.)
      else
        call MatAXPY(p_H, p_scal, p_TPkQPkT, SUBSET_NONZERO_PATTERN, ierr)      
      end if      
      
      call MatDestroy(p_TPkT, ierr)
      call MatDestroy(p_PkQPk, ierr)
      call MatDestroy(p_TPkQPkT, ierr)
    end do
    
    call MatDestroy(p_Q, ierr)
    
    if (Top .eq. 1) then
      nullify(p_T)
    else if (Top .eq. 2) then
      call MatDestroy(p_T, ierr)
    end if

  end subroutine get_VU_TDT
  
  subroutine get_VU_SD_DS(p_H, p_S, p_D)
#include <petsc/finclude/petsc.h>
    use petsc  
    use petsc_mod
    use petsc_wrapper
    use globals
    use k_on_demand_mod
    use error_handler
    implicit none  
    
    Mat, target :: p_H, p_S, p_D
    integer :: top
    
    Mat :: p_PkQPk, p_SPk, p_PkS, p_PkSP,  p_Q, p_tmp, p_SPkQPk, p_PkQPkS
    Vec :: v_Diag
    PetscScalar :: p_scal   
    integer :: ierr, ishell1
    PetscInt, pointer :: colmapS(:), colmapD(:), colmapQ(:)
    character(256) :: ss

    if (.not.l_dftu_sparse_projector) then
      call MatConvert(p_H, mattype_dense, MAT_INPLACE_MATRIX, p_H, ierr)
    end if
    
    call MatMatMult(p_S, p_D, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, &
   &  p_Q, ierr)
    call MatMatMult(p_D, p_S, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, &
   &  p_tmp, ierr)
    call MatAXPY(p_Q, p_one, p_tmp, SUBSET_NONZERO_PATTERN, ierr)
    call MatDestroy(p_tmp, ierr)
    
    p_scal = 0.5d0
    call MatScale(p_Q, p_scal, ierr)

     
    do ishell1 = 1, n_Ushells
            
      call MatMatMult(p_S, Uoperator(ishell1)%Pk, &
     &  MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, p_SPk, ierr)     
      call MatMatMult(Uoperator(ishell1)%Pk, p_S, &
     &  MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, p_PkS, ierr)
                   
      call MatConjugate(p_SPk, ierr)               
      call MatConjugate(p_PkS, ierr)               
      p_scal = 0.25d0 * Uoperator(ishell1)%U
      if (l_dftu_sparse_projector) then
        call petsc_aXpY(p_H, p_SPk, p_scal, PETSC_FALSE, .false.)
        call petsc_aXpY(p_H, p_PkS, p_scal, PETSC_FALSE, .false.)
      else
        call MatAXPY(p_H, p_scal, p_SPk, SUBSET_NONZERO_PATTERN, ierr)
        call MatAXPY(p_H, p_scal, p_PkS, SUBSET_NONZERO_PATTERN, ierr)            
      end if

      call MatDestroy(p_SPk, ierr)
      call MatDestroy(p_PkS, ierr)

          
       call MatMatMatMult(Uoperator(ishell1)%Pk, p_Q, Uoperator(ishell1)%Pk, &
      &  MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, p_PkQPk, ierr)
     
      
       call MatMatMult(p_S, p_PkQPk, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, &
      &  p_SPkQPk, ierr)
       call MatMatMult(p_PkQPk, p_S, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, &
      &  p_PkQPkS, ierr)
      
      call MatConjugate(p_SPkQPk, ierr)
      call MatConjugate(p_PkQPkS, ierr)      
      p_scal = -0.5d0 * Uoperator(ishell1)%U
      if (l_dftu_sparse_projector) then
        call petsc_aXpY(p_H, p_SPkQPk, p_scal, PETSC_FALSE, .false.)
        call petsc_aXpY(p_H, p_PkQPkS, p_scal, PETSC_FALSE, .false.)      
      else
        call MatAXPY(p_H, p_scal, p_SPkQPk, SUBSET_NONZERO_PATTERN, ierr)
        call MatAXPY(p_H, p_scal, p_PkQPkS, SUBSET_NONZERO_PATTERN, ierr)
      end if
      
    
      call MatDestroy(p_PkQPk, ierr)
      call MatDestroy(p_SPkQPk, ierr)
      call MatDestroy(p_PkQPkS, ierr)
    
    end do
    
    call MatDestroy(p_Q, ierr)

  end subroutine get_VU_SD_DS
  
  
  subroutine get_VU_Tr_symD(p_H, p_S, p_D)
#include <petsc/finclude/petsc.h>
    use petsc  
    use petsc_mod
    use petsc_wrapper
    use globals
    use k_on_demand_mod
    use error_handler
    implicit none  
    
    Mat, target :: p_H, p_S, p_D
    integer :: top
    
    Mat :: p_SD, p_DS, p_SPk, p_PkS, p_tmp1
    Vec :: v_Diag
    PetscScalar :: p_scal   
    integer :: ierr, ishell1
    PetscInt, pointer :: colmapS(:), colmapD(:), colmapQ(:)
    character(256) :: ss

    if (.not.l_dftu_sparse_projector) then
      call MatConvert(p_H, mattype_dense, MAT_INPLACE_MATRIX, p_H, ierr)
    end if
    
    call MatMatMult(p_S, p_D, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, &
   &  p_SD, ierr)
    call MatMatMult(p_D, p_S, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, &
   &  p_DS, ierr)


     
    do ishell1 = 1, n_Ushells
            
      call MatMatMult(p_S, Uoperator(ishell1)%Pk, &
     &  MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, p_SPk, ierr)     
      call MatMatMult(Uoperator(ishell1)%Pk, p_S, &
     &  MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, p_PkS, ierr)
          
          
      
      call MatMatMult(p_SPk, p_DS, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, &
     &  p_tmp1, ierr)
      call MatConjugate(p_tmp1, ierr)    
      p_scal = -0.25d0 * Uoperator(ishell1)%U
      if (l_dftu_sparse_projector) then
        call petsc_aXpY(p_H, p_tmp1, p_scal, PETSC_FALSE, .false.)
      else
        call MatAXPY(p_H, p_scal, p_tmp1, SUBSET_NONZERO_PATTERN, ierr)
      end if
     call MatDestroy(p_tmp1, ierr)
      
      call MatMatMult(p_PkS, p_DS, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, &
     &  p_tmp1, ierr)
      call MatConjugate(p_tmp1, ierr)      
      p_scal = -0.25d0 * Uoperator(ishell1)%U
      if (l_dftu_sparse_projector) then
        call petsc_aXpY(p_H, p_tmp1, p_scal, PETSC_FALSE, .false.)
      else
        call MatAXPY(p_H, p_scal, p_tmp1, SUBSET_NONZERO_PATTERN, ierr)
      end if
     call MatDestroy(p_tmp1, ierr)
      
      call MatMatMult(p_SD, p_SPk, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, &
     &  p_tmp1, ierr)
      call MatConjugate(p_tmp1, ierr)      
      p_scal = -0.25d0 * Uoperator(ishell1)%U
      if (l_dftu_sparse_projector) then
        call petsc_aXpY(p_H, p_tmp1, p_scal, PETSC_FALSE, .false.)
      else
        call MatAXPY(p_H, p_scal, p_tmp1, SUBSET_NONZERO_PATTERN, ierr)
      end if
      call MatDestroy(p_tmp1, ierr)
      
      call MatMatMult(p_SD, p_PkS, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, &
     &  p_tmp1, ierr)
      call MatConjugate(p_tmp1, ierr)      
      p_scal = -0.25d0 * Uoperator(ishell1)%U
      if (l_dftu_sparse_projector) then
        call petsc_aXpY(p_H, p_tmp1, p_scal, PETSC_FALSE, .false.)
      else
        call MatAXPY(p_H, p_scal, p_tmp1, SUBSET_NONZERO_PATTERN, ierr)
      end if
     call MatDestroy(p_tmp1, ierr)
     
     
!~ diagonal part                   
      call MatConjugate(p_SPk, ierr)               
      call MatConjugate(p_PkS, ierr)               
      p_scal = 0.25d0 * Uoperator(ishell1)%U
      if (l_dftu_sparse_projector) then
        call petsc_aXpY(p_H, p_SPk, p_scal, PETSC_FALSE, .false.)
        call petsc_aXpY(p_H, p_PkS, p_scal, PETSC_FALSE, .false.)
      else
        call MatAXPY(p_H, p_scal, p_SPk, SUBSET_NONZERO_PATTERN, ierr)
        call MatAXPY(p_H, p_scal, p_PkS, SUBSET_NONZERO_PATTERN, ierr)            
      end if

      call MatDestroy(p_SPk, ierr)
      call MatDestroy(p_PkS, ierr)     
      
    
    
    end do
    
    call MatDestroy(p_SD, ierr)
    call MatDestroy(p_DS, ierr)

  end subroutine get_VU_Tr_symD
  
  subroutine get_VU_Tr_symD_v2(p_H, p_S, p_D)
#include <petsc/finclude/petsc.h>
    use petsc  
    use petsc_mod
    use petsc_wrapper
    use globals
    use k_on_demand_mod
    use error_handler
    implicit none  
    
    Mat, target :: p_H, p_S, p_D
    integer :: top
    
    Mat :: p_SD, p_DS, p_SPk, p_PkS, p_tmp1
    Vec :: v_Diag
    PetscScalar :: p_scal   
    integer :: ierr, ishell1
    PetscInt, pointer :: colmapS(:), colmapD(:), colmapQ(:)
    character(256) :: ss

    if (.not.l_dftu_sparse_projector) then
      call MatConvert(p_H, mattype_dense, MAT_INPLACE_MATRIX, p_H, ierr)
    end if

     
    do ishell1 = 1, n_Ushells
            
      call MatMatMult(p_S, Uoperator(ishell1)%Pk, &
     &  MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, p_SPk, ierr)     
      call MatMatMult(Uoperator(ishell1)%Pk, p_S, &
     &  MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, p_PkS, ierr)
          
!~ offdiagonal part              
      call MatMatMatMult(p_SPk, p_D, p_SPk, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, &
     &  p_tmp1, ierr)     
      call MatConjugate(p_tmp1, ierr)    
      p_scal = -0.5d0 * Uoperator(ishell1)%U
      if (l_dftu_sparse_projector) then
        call petsc_aXpY(p_H, p_tmp1, p_scal, PETSC_FALSE, .false.)
      else
        call MatAXPY(p_H, p_scal, p_tmp1, SUBSET_NONZERO_PATTERN, ierr)
      end if
     call MatDestroy(p_tmp1, ierr)
      
      call MatMatMatMult(p_PkS, p_D, p_PkS, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, &
     &  p_tmp1, ierr)     
      call MatConjugate(p_tmp1, ierr)    
      p_scal = -0.5d0 * Uoperator(ishell1)%U
      if (l_dftu_sparse_projector) then
        call petsc_aXpY(p_H, p_tmp1, p_scal, PETSC_FALSE, .false.)
      else
        call MatAXPY(p_H, p_scal, p_tmp1, SUBSET_NONZERO_PATTERN, ierr)
      end if
     call MatDestroy(p_tmp1, ierr)

!~ "diagonal" part                   
      call MatConjugate(p_SPk, ierr)               
      call MatConjugate(p_PkS, ierr)               
      p_scal = 0.25d0 * Uoperator(ishell1)%U
      if (l_dftu_sparse_projector) then
        call petsc_aXpY(p_H, p_SPk, p_scal, PETSC_FALSE, .false.)
        call petsc_aXpY(p_H, p_PkS, p_scal, PETSC_FALSE, .false.)
      else
        call MatAXPY(p_H, p_scal, p_SPk, SUBSET_NONZERO_PATTERN, ierr)
        call MatAXPY(p_H, p_scal, p_PkS, SUBSET_NONZERO_PATTERN, ierr)            
      end if

      call MatDestroy(p_SPk, ierr)
      call MatDestroy(p_PkS, ierr)     
      
    
    
    end do
    
    call MatDestroy(p_SD, ierr)
    call MatDestroy(p_DS, ierr)

  end subroutine get_VU_Tr_symD_v2
  
  subroutine get_VU_Tr_PkDPk_bar(p_H, p_S, p_D)
#include <petsc/finclude/petsc.h>
    use petsc  
    use petsc_mod
    use petsc_wrapper
    use globals
    use k_on_demand_mod
    use error_handler
    implicit none  
    
    Mat, target :: p_H, p_S, p_D
    integer :: top
    
    Mat :: p_PkSPkDPkSPk,p_PkSPk ,p_tmp1
    PetscScalar :: p_scal   
    integer :: ierr, ishell1
    PetscInt, pointer :: colmapS(:), colmapD(:), colmapQ(:)
    character(256) :: ss

    if (.not.l_dftu_sparse_projector) then
      call MatConvert(p_H, mattype_dense, MAT_INPLACE_MATRIX, p_H, ierr)
    end if

     
    do ishell1 = 1, n_Ushells
            
      call MatMatMatMult(Uoperator(ishell1)%Pk, p_S, Uoperator(ishell1)%Pk, &
     &  MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, p_PkSPk, ierr) 
     
      call MatMatMatMult(p_PkSPk, p_D, p_PkSPk, &
     &  MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, p_PkSPkDPkSPk, ierr) 
    
      call MatConjugate(p_PkSPk, ierr)    
      p_scal = 0.5d0 * Uoperator(ishell1)%U
      if (l_dftu_sparse_projector) then
        call petsc_aXpY(p_H, p_PkSPk, p_scal, PETSC_FALSE, .false.)
      else
        call MatAXPY(p_H, p_scal, p_PkSPk, SUBSET_NONZERO_PATTERN, ierr)
      end if
      
    
      call MatConjugate(p_PkSPkDPkSPk, ierr)    
      p_scal = -Uoperator(ishell1)%U
      if (l_dftu_sparse_projector) then
        call petsc_aXpY(p_H, p_PkSPkDPkSPk, p_scal, PETSC_FALSE, .false.)
      else
        call MatAXPY(p_H, p_scal, p_PkSPkDPkSPk, SUBSET_NONZERO_PATTERN, ierr)
      end if
    

      call MatDestroy(p_PkSPk, ierr)
      call MatDestroy(p_PkSPkDPkSPk, ierr)     
            
    end do


  end subroutine get_VU_Tr_PkDPk_bar
  


  
  
  
end module dftu_mod
