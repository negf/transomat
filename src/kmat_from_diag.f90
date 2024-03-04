module kmat_from_diag
#include <petsc/finclude/petsc.h>
  implicit none
  integer :: Ntot_diag
  
  contains
  
    subroutine get_kmat_from_diag()
#include <petsc/finclude/petsc.h>
      use slepceps
      use petsc    
      use petsc_mod
      use petsc_wrapper
      use k_on_demand_mod
      use kinds
      use globals
      use ft_mod
      use slepc_mod
      use misc
      use error_handler
      use dftu_mod
            
      implicit none
      
      Mat, target :: p_mat_coeff, p_tmp1, p_tmp2, A, p_Dmat_full
      Mat, pointer :: p_Dmat_for_dftU
      Mat, allocatable :: p_mat_tmp(:)
      Vec, allocatable :: p_vec_ew(:)
      Vec, pointer :: p_vec_stride(:)
      Vec :: p_vec_ew_all
      integer :: ierr, ik, nrow, ncol, n, nloc, nloc1, nloc2, i, j, N_ew_total, m, k, nz_proc, &
        ii(1), jj(1), i1, i2, i3, iunit, ncols_nz, ncols_nz_max, i_max_ew, N_diag, nn, iamk, iew
      integer, allocatable :: i_rows(:), i_rows_shift(:), idx(:), displs(:), counts(:), &
        cols_nz(:,:), i_max_ewk(:), i_row_nz(:)
      PetscScalar :: z_ew_n, c_jn, fac, p11, d11, nek, ne
      PetscScalar, pointer :: p_x1(:),p_x2(:), pp_coeff(:,:), local_Dmat_full(:, :)
      PetscReal, allocatable :: ew(:)
      PetscScalar, allocatable :: coeff_batch(:), coeff_batch_j(:), cols_buf(:), local_Dmat(:), local_Dmat2(:)
      PetscReal :: N_elec, N_elec_target, Ef, E1, E2, N_left, fermi_fac, ew_n, norm, ntrace, &
     &  norm01, norm02, norm03 
      integer(8) :: counti, count_rate, countf, ipos, counter(3,256)
      character(strln) :: file_ev, outstr, testfile, file_dmat_full, istr(3)
      logical :: l_exist

      counter = 0

      call MatGetSize(p_h00_ik_cc(1), nrow, ncol, ierr)
      call MatGetOwnershipRange(p_h00_ik_cc(1), nloc1, nloc2, ierr)
      nloc = nloc2 - nloc1
      
      allocate(i_rows(nloc), i_rows_shift(nloc), coeff_batch(nloc), &
        coeff_batch_j(ncol), displs(nprocs), counts(nprocs), stat = ierr)
      
      do i = 1, nloc
        i_rows(i) = nloc1 + i - 1
      end do

      counts(inode+1) = nloc
      displs(inode+1) = nloc1
      
    
      call MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, counts, 1, MPI_INTEGER, PETSC_COMM_WORLD, ierr)
      call MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, displs, 1, MPI_INTEGER, PETSC_COMM_WORLD, ierr)
      
      N_ew_total = nk_c * nrow
      allocate(p_vec_ew(nk_c), ew(N_ew_total), stat = ierr)

      allocate(p_mat_tmp(1), stat = ierr)
    
      call MatDuplicate(p_h00_ik_cc(1), MAT_DO_NOT_COPY_VALUES, p_mat_tmp(1), ierr)
!~       call MatSetOption(p_mat_tmp(1), MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE, ierr)
      call MatSetOption(p_mat_tmp(1), MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE, ierr)
      call MatSetOption(p_mat_tmp(1), MAT_IGNORE_OFF_PROC_ENTRIES, PETSC_FALSE, ierr)
      call MatSetOption(p_mat_tmp(1), MAT_NO_OFF_PROC_ZERO_ROWS, PETSC_FALSE, ierr)
      
      call petsc_get_densemat(p_mat_tmp(1), p_mat_coeff, mattype_dense)
!~       if (l_dftu) call petsc_get_densemat(p_mat_tmp(1), p_Dmat_full, mattype_dense)
      

      N_elec_target = n_electrons_ecc(spin_channel)
      write(pstr_out,fmt='(A,i8,3e24.12)') "N_elec_target ", spin_channel, &
     &  N_elec_target, n_electrons_ecc 
      call petsc_print_master()      
      
!~       N_diag = ceiling(N_elec_target * 1.4d0)
!~       N_diag = min(N_diag, nrow)
      N_diag = nrow
      Ntot_diag = N_diag * nk_c
            
      do ik = 1, nk_c       
        iik = ik
        call system_clock(counti, count_rate)
        write(pstr_out,fmt='(A,i8,3e24.12,i4)') "# k: ", nk_c-ik, kp_c(1:3,ik), wkp_c(ik)        
        call petsc_print_master(.false.)


        
        if (l_k_on_demand) call k_on_demand(ik, 3)
        
!~         if (hartree_shift .ne. 0d0) then
!~           write(0,*) "hartree_shift ", hartree_shift
!~           fac = -hartree_shift          
!~           call MatAXPY(p_h00k_cc(ik), fac , p_s00k_cc(ik), SAME_NONZERO_PATTERN, ierr)        
!~         end if

!~         if (l_dftu) then
        
!~           file_dmat_full = "D_"//trim(spin_str)//"_full_"//trim(int2str(ik))//".dat"
!~           inquire(file = trim(file_dmat_full), exist=l_exist)
!~           if ((l_exist) .and. (dftu_projector .ge. 2) .and. &
!~          &  (.not.l_dftu_sparse_projector)) then
!~             call petsc_mat_direct_load(p_Dmat_full, file_dmat_full, ierr)
!~             p_Dmat_for_dftU => p_Dmat_full
!~           else
!~             p_Dmat_for_dftU => p_k00k_cc(ik)
!~           end if
!~           call MatNorm(p_h00k_cc(ik), NORM_FROBENIUS, norm01, ierr)
!~           if (dftu_projector .eq. 1) then
!~             call get_VU_SpD_DpS_v2(p_h00k_cc(ik), p_s00k_cc(ik), p_Dmat_for_dftU)
!~           else if (dftu_projector .eq. 2) then
!~             call get_VU_TDT(p_h00k_cc(ik), p_s00k_cc(ik), p_Dmat_for_dftU, 2)
!~           else if (dftu_projector .eq. 3) then
!~             call get_VU_TDT(p_h00k_cc(ik), p_s00k_cc(ik), p_Dmat_for_dftU, 1)
!~           else if (dftu_projector .eq. 4) then
!~             call get_VU_SD_DS(p_h00k_cc(ik), p_s00k_cc(ik), p_Dmat_for_dftU)
!~           else if (dftu_projector .eq. 5) then
!~             call get_VU_Tr_symD_v2(p_h00k_cc(ik), p_s00k_cc(ik), p_Dmat_for_dftU)
!~           else if (dftu_projector .eq. 6) then
!~             call get_VU_Tr_PkDPk_bar(p_h00k_cc(ik), p_s00k_cc(ik), p_Dmat_for_dftU)
!~           else 
!~             write (errormsg, fmt='(A,i8)') "unknown DFTU projector ", dftu_projector
!~             call error()
!~           end if
!~           call MatHermitianTranspose(p_h00_ik_cc(1), MAT_INITIAL_MATRIX, A, ierr)
!~           call MatAXPY(A,p_minus1,p_h00_ik_cc(1),SAME_NONZERO_PATTERN, ierr)
!~           call MatNorm(A, NORM_FROBENIUS, norm, ierr)
!~           call MatNorm(p_h00_ik_cc(1), NORM_FROBENIUS, norm03, ierr)
!~           call MatDestroy(A,ierr)
!~           write (pstr_out, fmt='(A,3e18.8,l)') "Norm(K(k)-K(k)')", norm,norm01,norm03,l_exist
!~           call petsc_print_master()
!~         end if
          
        
        call MatCreateVecs(p_h00k_cc(ik), p_vec_ew(ik), PETSC_NULL_VEC, ierr) 
        call VecShift(p_vec_ew(ik), dcmplx(huge(1d0), 1d-64), ierr)
        call MatZeroEntries(p_mat_coeff, ierr)
        counter(3,6) = 0
        call system_clock(counter(1,6), count_rate)
        call diag_mat2(p_h00k_cc(ik), p_s00k_cc(ik), p_vec_ew(ik), p_mat_coeff, N_diag, eps_which_in=EPS_SMALLEST_REAL)      
        call system_clock(counter(2,6))          
        counter(3,6) = counter(3,6) + (counter(2,6) - counter(1,6))                
!~         call diag_mat2(p_h00k_cc(ik), p_s00k_cc(ik), p_vec_ew(ik), PETSC_NULL_MAT, N_diag, eps_which_in=EPS_SMALLEST_REAL)      
!~         call diag_mat(p_h00k_cc(ik), p_s00k_cc(ik), p_vec_ew(ik), PETSC_NULL_MAT, nmu_ecc, eps_which_in=EPS_SMALLEST_REAL)      

        write(file_ev,fmt='(A,i0,A)') "ev_mat_", ik,".dat"
        call petsc_mat_direct_save(p_mat_coeff, file_ev, ierr)
        
        call system_clock(countf)
        write (pstr_out, fmt='(X,2e24.12)') real(countf - counti, 8)/real(count_rate, 8), &
          real(counter(3,6), 8)/real(count_rate, 8)
        call petsc_print_master()
      end do
      
      write (pstr_out, fmt='(A)') repeat("-",80)
      call petsc_print_master()
      
      n = 0
      if (inode .eq. 0) n = N_ew_total
      call VecCreateMPI(PETSC_COMM_WORLD, n, N_ew_total, p_vec_ew_all, ierr)
      
      call VecGetArrayF90(p_vec_ew_all, p_x2,ierr)
      
      do ik = 1, nk_c
        
        call VecGetArrayF90(p_vec_ew(ik), p_x1, ierr)
        i_rows_shift = i_rows + nrow * (ik-1)
        call VecSetValues(p_vec_ew_all, nloc, i_rows_shift, p_x1, INSERT_VALUES, ierr)
        call VecRestoreArrayF90(p_vec_ew(ik), p_x1,ierr)
   
      end do
      call VecAssemblyBegin(p_vec_ew_all, ierr)
      call VecAssemblyEnd(p_vec_ew_all, ierr)
      
      if (maxval(abs(aimag(p_x2))).ge.1d-14) then
        write (pstr_out, fmt='(A, 2e24.12)') "warning max(abs(imag(ew))) ", &
          maxval(abs(aimag(p_x2)))
        call petsc_print_master()
      end if
      
      allocate(idx(N_ew_total), stat = ierr)
      
      if (inode.eq.0) then
        
        ew=real(p_x2, dp)
          
        allocate(i_max_ewk(nk_c), stat = ierr)
          
        do i = 1, N_ew_total
          idx(i) = i - 1
        end do
    
        call PetscSortRealWithPermutation(N_ew_total, ew, idx, ierr)
        open (newunit = j, file = "eigall_"//trim(spin_str)//".dat", action="write", &
       &  status = "replace")
       
        nn = N_ew_total / nk_c ! this should be nrow (ncol)
      
        do i = 1, N_ew_total
          i1 = ceiling((idx(i)+0.5d0)/nn)
          write(unit=j,fmt='(2i8,3e24.12,i8)') i, idx(i), p_x2(idx(i)+1), &          
         & real(wkp_c(i1), 8) / real(nktot, 8), Ntot_diag
        end do
        close(j)
        if (l_diag_fix_ef) then
          Ef = (ef_l + ef_r) * 0.5d0
        else
          E2 = ew(idx(Ntot_diag)+1)
          E1 = ew(1)
          Ef = ew(idx(Ntot_diag)*1) !(E2 + E1) * 0.5d0
          k = 0
          
          do
            k=k+1
            N_elec = 0d0
  
            call Nelec_min_fct(Ef, N_left, N_elec_target, ew, wkp_c, idx, temperature_el)          
  
            if (abs(N_left).le.1d-14) exit
            
            if (N_left.gt.0d0) then            
              E1 = Ef
              Ef = (E2 + E1) * 0.5d0
            else            
              E2 = Ef
              Ef = (E2 + E1) * 0.5d0        
            end if
            if (k.ge.1000) exit
          end do      
        end if
        
        call Nelec_min_fct(Ef, N_left, 0d0, ew, wkp_l, idx, temperature_el, i_max_ewk)  
        i_max_ew = maxval(i_max_ewk)
        i_max_ew = ceiling(i_max_ew * 1.1d0)
      end if
                              
      call VecRestoreArrayF90(p_vec_ew_all, p_x2,ierr)      
                
      call MPI_Bcast(Ef, 1, MPI_DOUBLE_PRECISION, 0, PETSC_COMM_WORLD, ierr)
      call MPI_Bcast(N_left, 1, MPI_DOUBLE_PRECISION, 0, PETSC_COMM_WORLD, ierr)
      call MPI_Bcast(i_max_ew, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
      call MPI_Bcast(idx, N_ew_total, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
      
        
      write (pstr_out, fmt='(A, 2e24.12, i8)') "Fermi energy ", Ef, -N_left, i_max_ew
      call petsc_print_master()
      write (pstr_out, fmt='(A)') repeat("-",80)
      call petsc_print_master()

      
      do i2 = -ncell_c(2), ncell_c(2)
        do i1 = -ncell_c(1), ncell_c(1)
          do i3 = -ncell_c(3), ncell_c(3)
            call MatZeroEntries(p_dmatxy(i1, i2, i3), ierr)          
          end do
        end do
      end do
      
      call system_clock(counti, count_rate)
      ncols_nz_max = 0
      do i = 1, nloc        
        call MatGetRow(p_mat_tmp(1), i_rows(i), ncols_nz, PETSC_NULL_INTEGER, &
          PETSC_NULL_SCALAR, ierr)
        ncols_nz_max = max(ncols_nz_max, ncols_nz)
        call MatRestoreRow(p_mat_tmp(1), i_rows(i), ncols_nz, PETSC_NULL_INTEGER, &
          PETSC_NULL_SCALAR, ierr)
      end do

      allocate(cols_nz(ncols_nz_max, nloc), cols_buf(ncols_nz_max), &
        i_row_nz(nloc), stat = ierr)

      do i = 1, nloc        
        call MatGetRow(p_mat_tmp(1), i_rows(i), i_row_nz(i), cols_nz(:,i), &
          PETSC_NULL_SCALAR, ierr)
        call MatRestoreRow(p_mat_tmp(1), i_rows(i), ncols_nz, PETSC_NULL_INTEGER, &
          PETSC_NULL_SCALAR, ierr)
      end do

      call system_clock(countf)
      write (pstr_out, fmt='(A,i8,e24.12)') "ncols_nz_max ", ncols_nz_max, &
        real(countf - counti, 8)/real(count_rate, 8)
      call petsc_print_master()  
      write (pstr_out, fmt='(A)') repeat("-",80)
      call petsc_print_master()
    
      
      Ntrace = 0d0
      call MatDenseGetArrayF90(p_mat_coeff, pp_coeff, ierr) 
!~       if (l_dftu) then        
!~         call MatDenseGetArrayF90(p_Dmat_full, local_Dmat_full, ierr)         
!~       end if
      
      if (inode.eq.0) open(newunit=iunit,file="eigs_"//trim(spin_str)//".dat", &
     &  action="write",status="replace")
      
        call MatSetOption(p_mat_tmp(1), MAT_NO_OFF_PROC_ENTRIES, PETSC_TRUE, ierr)
!~         call MatSetOption(p_mat_tmp(1), MAT_IGNORE_OFF_PROC_ENTRIES, PETSC_TRUE, ierr)
        call MatSetOption(p_mat_tmp(1), MAT_NEW_NONZERO_LOCATION_ERR , PETSC_TRUE, ierr)      

      nz_proc = sum(i_row_nz(1:nloc))
      allocate(local_Dmat(nz_proc))
      ne = 0d0
      do ik = 1, nk_c
        iik = ik
        counter = 0d0
        call system_clock(counti, count_rate)
        write(pstr_out,fmt='(A,i8,3e24.12,i4)') "# k: ", nk_c-ik, kp_c(1:3,ik), wkp_c(ik)
        if (inode.eq.0) write(iunit,fmt='(A)') trim(pstr_out)
        call petsc_print_master(.false.)
        
!~         if (l_k_on_demand) call k_on_demand(ik, 3)
        
        call MatZeroEntries(p_mat_coeff, ierr)
!~         call diag_mat2(p_h00k_cc(ik), p_s00k_cc(ik), p_vec_ew(ik), p_mat_coeff, i_max_ew, eps_which_in=EPS_SMALLEST_REAL)      
!~         call diag_mat(p_h00k_cc(ik), p_s00k_cc(ik), p_vec_ew(ik), p_mat_coeff, nmu_ecc, eps_which_in=EPS_SMALLEST_REAL)      

        write(file_ev,fmt='(A,i0,A)') "ev_mat_", ik,".dat"
        call system_clock(counter(1,5), count_rate)
        call petsc_mat_direct_load(p_mat_coeff, file_ev, ierr)
        call system_clock(counter(2,5))          
        counter(3,5) = counter(3,5) + (counter(2,5) - counter(1,5))

        local_Dmat = 0d0
!~         local_Dmat2 = 0d0
!~         if (l_dftu) local_Dmat_full = 0d0
        
        call MatZeroEntries(p_mat_tmp(1), ierr) ! zero as we add values          
        p11 = 0d0
        iew = 0
        
        do nn = 1, N_ew_total
          n = idx(nn) + 1 
          
          iamk = ceiling(1d0/real(nrow,dp)*(n)-1d0/(real(nrow,dp)*2d0))
          
          if (ik.ne.iamk) cycle
          
          n = n - (nrow * (ik - 1))
          
          iew = iew + 1
          call petsc_vec_getvalue(n - 1, z_ew_n, p_vec_ew(ik), .true.)             
          ew_n = z_ew_n ! EW should all be real.
          fermi_fac = fermi_function(ew_n, Ef, temperature_el)   
          if (inode.eq.0) write(iunit,fmt='(i8,2e24.12,E24.12E4)') iew, z_ew_n, fermi_fac
          if (fermi_fac .le. 1d-13) cycle
          call system_clock(counter(1,1), count_rate)
          coeff_batch_j = 0d0
          call MPI_Allgatherv(pp_coeff(1:nloc,n), nloc, MPI_DOUBLE_COMPLEX, coeff_batch_j, &
            counts, displs, MPI_DOUBLE_COMPLEX, PETSC_COMM_WORLD, ierr)
          call system_clock(counter(2,1))          
          counter(3,1) = counter(3,1) + (counter(2,1) - counter(1,1))
          
!~           i1 = 0
!~           call system_clock(counter(1,2), count_rate)
!~           do i = 1, nloc            
!~             fac = fermi_fac * pp_coeff(i,n)
!~             do j = 1, i_row_nz(i)
!~               i1 = i1 + 1
!~               k = cols_nz(j,i) + 1
!~               c_jn = dconjg(coeff_batch_j(k)) * fac
!~               local_Dmat(i1) = local_Dmat(i1) + c_jn
!~             end do                        
!~           end do
!~           call system_clock(counter(2,2))
!~           counter(3,2) = counter(3,2) + (counter(2,2) - counter(1,2))                        
          
          i1 = 0
          call system_clock(counter(1,2), count_rate)
          do i = 1, nloc            
            fac = fermi_fac * pp_coeff(i,n)
            i2 = 1
            do j = 1, ncol
              if ((cols_nz(i2,i) + 1).eq.j) then                
                i1 = i1 + 1
                c_jn = dconjg(coeff_batch_j(j)) * fac
                local_Dmat(i1) = local_Dmat(i1) + c_jn
!~                 if (l_dftu) local_Dmat_full(i, j) = local_Dmat_full(i, j) + c_jn
                i2 = i2 + 1                                                
!~               else if ((l_dftu) .and. (dftu_projector.ge.2) .and. &
!~              &  (.not.l_dftu_sparse_projector)) then
!~                 c_jn = dconjg(coeff_batch_j(j)) * fac
!~                 local_Dmat_full(i, j) = local_Dmat_full(i, j) + c_jn
              end if
            end do                        
          end do
          call system_clock(counter(2,2))
          counter(3,2) = counter(3,2) + (counter(2,2) - counter(1,2))                        
              
        end do  ! nn
        
!~         if ((l_dftu) .and. (dftu_projector.ge.2) .and. (.not.l_dftu_sparse_projector)) then
!~           file_dmat_full = "D_"//trim(spin_str)//"_full_"//trim(int2str(ik))//".dat"
!~           call petsc_mat_direct_save(p_Dmat_full, file_dmat_full, ierr)          
!~         end if
        
!~         write(0,fmt='(i8,A,3e24.12)') inode, "check Dmats ",maxval(abs(local_Dmat2-local_Dmat)), &
!~        &  maxval(abs(local_Dmat2-local_Dmat)),maxval(abs(local_Dmat2))
       
        call system_clock(counter(1,3), count_rate)
        i1 = 1
        do i = 1, nloc
          i2 = i1 + i_row_nz(i) - 1
          call MatSetValues(p_mat_tmp(1), 1, i_rows(i), i_row_nz(i), cols_nz(1:i_row_nz(i),i), &
            local_Dmat(i1:i2), INSERT_VALUES, ierr)
          i1 = i2 + 1
        end do
        call system_clock(counter(2,3))
        counter(3,3) = counter(3,3) +  (counter(2,3) - counter(1,3))        
            
        call system_clock(counter(1,4), count_rate)
        call MatAssemblyBegin(p_mat_tmp(1), MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyEnd(p_mat_tmp(1), MAT_FINAL_ASSEMBLY, ierr)         
        call system_clock(counter(2,4))
        counter(3,4) = counter(3,4) + (counter(2,4) - counter(1,4))

                
!~         call test_ne2(p_mat_tmp, nek)        
!~         ne = ne + nek
        
        call fourier3d_back_add(p_mat_tmp, p_dmatxy(:, :, :), kp_c(1:3, ik:ik), &
          wkp_c(ik:ik), 1, dlat_c, ncell_c(1), ncell_c(2), ncell_c(3), 0, 0, PETSC_FALSE)
        
        call system_clock(countf)
        write (pstr_out, fmt='(X,e24.12)') real(countf - counti, 8)/real(count_rate, 8)
        call petsc_print_master()
        write (pstr_out, fmt='(X,e24.12)') real(counter(3,1), 8)/real(count_rate, 8)
        call petsc_print_master()
        write (pstr_out, fmt='(X,e24.12)') real(counter(3,2), 8)/real(count_rate, 8)
        call petsc_print_master()
        write (pstr_out, fmt='(X,e24.12)') real(counter(3,3), 8)/real(count_rate, 8)
        call petsc_print_master()
        write (pstr_out, fmt='(X,e24.12)') real(counter(3,4), 8)/real(count_rate, 8)
        call petsc_print_master()
        write (pstr_out, fmt='(X,e24.12)') real(counter(3,5), 8)/real(count_rate, 8)
        
        call petsc_print_master()
        
      end do 
      
!~       if (l_dftu) then
!~         call MatDenseRestoreArrayF90(p_Dmat_full, local_Dmat_full, ierr)         
!~         call MatDestroy(p_Dmat_full, ierr)
!~       end if
      
!~       ne = ne / nktot
!~       write (pstr_out, fmt='(A,2e24.12)') "Nktot ",ne
!~       call petsc_print_master()
!~       write (pstr_out, fmt=*) "wkp_c ",wkp_c
!~       call petsc_print_master()
      call MatSetOption(p_mat_tmp(1), MAT_NO_OFF_PROC_ENTRIES, PETSC_FALSE, ierr)
      if (inode.eq.0) close(iunit)
      call MatDenseRestoreArrayF90(p_mat_coeff,pp_coeff, ierr)    
            
      fac = 1d0/(nktot)
      ne = 0d0
      do i1 = -ncell_c(1), ncell_c(1)
        do i2 = -ncell_c(2), ncell_c(2)
          do i3 = -ncell_c(3), ncell_c(3)
            call MatScale(p_dmatxy(i1, i2, i3), fac, ierr)
            call MatMatMult(p_dmatxy(i1,i2,i3), p_s00_diag(-i1,-i2,-i3)%p, MAT_INITIAL_MATRIX, &
          &  PETSC_DEFAULT_REAL, p_tmp1, ierr)
            call MatGetTrace(p_tmp1, nek, ierr)
            ne = ne + nek  
            call MatDestroy(p_tmp1, ierr)
          end do
        end do
      end do
      write (pstr_out, fmt='(A,2e24.12)') "Nrtot ", ne
      call petsc_print_master()       
      
  
      do ik = 1, nk_c
        call VecDestroy(p_vec_ew(ik), ierr)
      end do    
  
      deallocate(p_vec_ew, stat = ierr)

      call VecDestroy(p_vec_ew_all, ierr)
      
      call MatDestroy(p_mat_coeff, ierr)
      call MatDestroy(p_mat_tmp(1), ierr)
      
!~ save Fermi energy to sys_info.dat using MPI/IO as we can change a text stream easily
      call MPI_FILE_OPEN(MPI_COMM_WORLD, "sys_info.dat", MPI_MODE_WRONLY + MPI_MODE_CREATE, &
        MPI_INFO_NULL, iunit, ierr)
        
      if (l_ionode) then
        ipos = 3*3*24 + 3 + 3 * 22 + 1 + 22 + 1
        write(outstr,fmt='(e24.12, X, l, A)') Ef, .true., NEW_LINE('A')
        call MPI_FILE_WRITE_AT(iunit, ipos, trim(outstr), len(trim(outstr)), &
          MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)                        
      end if
      
      
      call MPI_FILE_CLOSE(iunit, ierr)
        
    end subroutine get_kmat_from_diag
  
    function fermi_function(x,Ef,T)
      use globals, only : kb
      use kinds
      implicit none
      
      real(dp) :: fermi_function
      real(dp), intent(in) :: x, Ef, T
      
      real(dp) :: beta, xx
      
      beta = 1d0/(T * kb) 
      xx = (x-Ef)*beta     
      fermi_function = 1d0 / ( exp(xx) + 1d0)
    
    end function fermi_function  

    subroutine Nelec_min_fct(x, f, N, ew, wk, idx, T, i_max_ewk)
#include <petsc/finclude/petsc.h>
      use petsc
      use globals
      implicit none
      
      PetscReal :: x, f, N, T, ew(:)
      integer :: wk(:), idx(:)
      integer, optional :: i_max_ewk(:)
      
      integer :: N_ew_total, i, ik, nrow      
      PetscReal :: N_elec, fac
      
      N_ew_total = size(ew)
      nrow  = N_ew_total / size(wk)
      
      N_elec = 0d0
      if (present(i_max_ewk)) i_max_ewk = 0
      do i = 1, Ntot_diag ! N_ew_total
!~ this should be save as long as the roundoff error is < 1.0/nrow*0.5d0      
        ik = ceiling(1d0/real(nrow,dp)*(idx(i)+1)-1d0/(real(nrow,dp)*2d0))        
        fac = fermi_function(ew(idx(i)+1),x,T)
        if ((fac.ge.1d-13).and.(present(i_max_ewk))) &
          i_max_ewk(ik) = i_max_ewk(ik) + 1
        N_elec = N_elec + fac*wk(ik)
      end do          
      N_elec = N_elec / nktot
      f = N - N_elec
    end subroutine Nelec_min_fct    
    

end module kmat_from_diag
