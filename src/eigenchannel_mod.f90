module eigenchannel_mod
  implicit none
  
  contains
  
    subroutine calc_eigenchannel_wf(t_i)
#include <slepc/finclude/slepceps.h>
      use slepceps
      use petscmat
      use petsc_mod
      use petsc_wrapper
      use globals
      use integrator_mod
      use slepc_mod
      
      Vec :: t_i
      
      Mat :: p_Ax, p_Dx, p_GammaX_full, p_rhs, p_Gr, p_Ga, preGamma, p_tmp1, p_tmp2, &
     &  p_Cx, p_trans
      Mat, pointer :: p_gr_inv
      Vec :: v_lambda_x, ti_x
      PetscScalar :: p_val, p_sum
      integer :: ierr, ii, nev_calc
      PetscBool :: l_ishermitian
      PetscReal :: eps_hermitian
      PetscScalar, pointer :: ev_pointer(:, :)
      character(256) :: outfile
      
      p_gr_inv => p_invGr ! should be still avaiable from previous call but now we need the whole matrix G
      
      call petsc_get_a_with_b_c(p_Gr, p_gr_inv, p_gr_inv, mattype_dense)
      
      solver_mode = 1 ! for now assume solver mode 1
      if (solver_mode .eq. 2) then
        call petsc_get_a_with_b_c(p_rhs, p_gr_inv, p_gr_inv, mattype_dense, .false.)
        call petsc_one_mat(p_rhs, 0, nmu_c -1) ! prepare RHS nmu_c,nmu_r ( .., 0 0 .. 1 0, 0 0 .. 0 1 )   C style 1=0
        call petsc_call_solver(p_gr_inv, p_rhs, p_Gr, matsolvertype_cc, solver_mode)
      else if ((solver_mode .eq. 1).or.(solver_mode .eq. 3).or.(solver_mode .eq. 4)) then
        call petsc_get_a_with_b_c(p_rhs, p_gr_inv, p_gr_inv, mattype_sparse)
        call petsc_one_mat(p_rhs,  0, nmu_c -1,2) ! prepare RHS nmu_c,nmu_r ( .., 0 0 .. 1 0, 0 0 .. 0 1 )   C style 1=0
        call petsc_call_solver(p_gr_inv, p_rhs, p_Gr, matsolvertype_cc, solver_mode)
      end if
  
      call MatDestroy(p_rhs, ierr)              
      call MatHermitianTranspose(p_Gr, MAT_INITIAL_MATRIX, p_Ga, ierr)
        
!~ first for Ar=Gr*GammaR*Ga

      write(outfile, fmt='(i8)') iik
      outfile = "S_"//trim(adjustl(outfile))//".dat"
      call dump_nonzero_structure(p_s00k_cc(iik), outfile, PETSC_COMM_WORLD, 0, 1) 
      write(outfile, fmt='(i8)') iik
      outfile = "Gr_"//trim(adjustl(outfile))//".dat"
      call dump_nonzero_structure(p_Gr, outfile, PETSC_COMM_WORLD, 0, 1) 
      write(outfile, fmt='(i8)') iik
      outfile = "Ga_"//trim(adjustl(outfile))//".dat"
      call dump_nonzero_structure(p_Ga, outfile, PETSC_COMM_WORLD, 0, 1) 
      write(outfile, fmt='(i8)') iik
      outfile = "GammaL_"//trim(adjustl(outfile))//".dat"
      call dump_nonzero_structure(p_GammaL, outfile, PETSC_COMM_WORLD, 0, 1) 
      write(outfile, fmt='(i8)') iik
      outfile = "GammaR_"//trim(adjustl(outfile))//".dat"
      call dump_nonzero_structure(p_GammaR, outfile, PETSC_COMM_WORLD, 0, 1) 


      call petsc_get_a_with_b_c(p_GammaX_full, p_gr_inv, p_gr_inv, mattype_sparse)
      call petsc_get_alloc_preG(p_GammaX_full, preGamma, nmu_l, 0)
      call MatDestroy(p_GammaX_full, ierr)
      call petsc_get_a_with_b_c(p_GammaX_full, p_gr_inv, p_gr_inv, mattype_sparse)
      call MatPreallocatorPreallocate(preGamma, PETSC_TRUE, p_GammaX_full, ierr)
      call MatDestroy(preGamma, ierr)

       
      call petsc_add_sub_B_to_A(p_GammaL, p_GammaX_full, 0, 0, p_one,&
     &  INSERT_VALUES, PETSC_FALSE)
     
       
      call MatMatMatMult(p_Gr, p_GammaX_full, p_Ga, MAT_INITIAL_MATRIX, &
     &  PETSC_DEFAULT_REAL, p_tmp1, ierr)
     
   

      call MatMatMatMult(p_s00k_cc(iik), p_tmp1, p_s00k_cc(iik), MAT_INITIAL_MATRIX, &
     &  PETSC_DEFAULT_REAL, p_Ax, ierr)

        
     
      call MatDestroy(p_GammaX_full, ierr)
      call MatDestroy(p_tmp1, ierr)


      

      call petsc_get_a_with_b_c(p_Dx, p_Ax, p_Ax, mattype_dense)
      call MatCreateVecs(p_Ax, v_lambda_x, PETSC_NULL_VEC, ierr)
      call diag_mat2(p_Ax, p_s00k_cc(iik), v_lambda_x, p_Dx, nmu_c, 1d-12, &
     &  eps_which_in = EPS_LARGEST_REAL, nev_calc = nev_calc)

      write (pstr_out, fmt='(A, i8)') "nev_calc ", nev_calc
      call petsc_print_master()              
          
      call MatDenseGetArrayF90(p_Dx, ev_pointer, ierr)
      do ii = 0, nmu_c - 1       
        call petsc_vec_getvalue(ii, p_val, v_lambda_x)        
        if (real(p_val, 8) .ge. 1d-32) then
!~           write (pstr_out, fmt='(A, i8, 2e24.12)') "EW ", ii, p_val
!~           call petsc_print_master()         
          p_val = dsqrt(real(p_val, 8) / (2d0 * pi))
        else
          p_val = 0d0
        end if
        ev_pointer(:, ii + 1) = ev_pointer(:, ii + 1) * p_val
      end do
      call MatDenseRestoreArrayF90(p_Dx, ev_pointer, ierr)

       
      call petsc_get_a_with_b_c(p_GammaX_full, p_gr_inv, p_gr_inv, mattype_sparse)
      call petsc_get_alloc_preG(p_GammaX_full, preGamma, 0, nmu_r)
      call MatDestroy(p_GammaX_full, ierr)
      call petsc_get_a_with_b_c(p_GammaX_full, p_gr_inv, p_gr_inv, mattype_sparse)
      call MatPreallocatorPreallocate(preGamma, PETSC_TRUE, p_GammaX_full, ierr)
      call MatDestroy(preGamma, ierr)
    
      call petsc_add_sub_B_to_A(p_GammaR, p_GammaX_full, nmu_c - nmu_r, nmu_c - nmu_r, &
     &  p_one, INSERT_VALUES, PETSC_FALSE)
     

     
      call MatMatMult(p_GammaX_full, p_Dx, MAT_INITIAL_MATRIX, &
     &  PETSC_DEFAULT_REAL, p_tmp1, ierr)
      call MatDestroy(p_GammaX_full, ierr)
      call MatHermitianTranspose(p_Dx, MAT_INPLACE_MATRIX, p_Dx, ierr)
      call MatMatMult(p_Dx, p_tmp1, MAT_INITIAL_MATRIX, &
     &  PETSC_DEFAULT_REAL, p_GammaX_full, ierr)
      call MatDestroy(p_tmp1, ierr)
      
      p_val = 2d0 * pi
      call MatScale(p_GammaX_full, p_val, ierr)
      
      call petsc_get_a_with_b_c(p_Cx, p_GammaX_full, p_GammaX_full, mattype_dense)
      call MatCreateVecs(p_Cx, ti_x, PETSC_NULL_VEC, ierr)
      
      
      
      call diag_mat2(p_GammaX_full, PETSC_NULL_MAT, ti_x , p_Cx, nmu_c, 1d-12, &
     &  eps_which_in = EPS_LARGEST_REAL, nev_calc = nev_calc)
      write (pstr_out, fmt='(A, i8)') "nev_calc ", nev_calc
      call petsc_print_master()              
   
       do ii = 0, nmu_c - 1       
        call petsc_vec_getvalue(ii, p_val, ti_x)        
        if (real(p_val, 8) .ge. 1d-8) then
          write (pstr_out, fmt='(A, i8, 4e24.12)') "ti ", ii, p_val
          call petsc_print_master()         
        end if
      end do


      call VecDestroy(v_lambda_x, ierr)        
      call VecDestroy(ti_x, ierr)        
      call MatDestroy(p_Ax, ierr)        
      call MatDestroy(p_Dx, ierr)        
      call MatDestroy(p_Gr, ierr)        
      call MatDestroy(p_Ga, ierr)        
      call MatDestroy(p_GammaX_full, ierr)        
      
        
    end subroutine calc_eigenchannel_wf
    
end module eigenchannel_mod
