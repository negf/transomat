module fix_mat_mod

  implicit none
  
  contains
  
    subroutine fix_mat(p_c00, p_l00, p_r00, p_c10, p_l10, p_c01,  p_r01)
#include <petsc/finclude/petsc.h>
      use slepceps
      use petscmat
      use petsc_mod
      use petsc_wrapper
      use kinds
      use globals
      use misc
      use error_handler

      implicit none
    
      Mat, allocatable :: p_c00(:, :), p_l00(:, :), p_r00(:, :)
      Mat, allocatable, optional :: p_c10(:, :), p_l10(:, :), p_c01(:, :), p_r01(:, :)
      integer :: j1, j2, i1, i2, ierr
      PetscScalar :: p_scal
    
      p_scal = p_one
      if (l_use_sigma_l) then
        pstr_out = "add electrode L to C"; call petsc_print_master()
      
      
      
        do i2 = -ncell_l(2), ncell_l(2)
          do i1 = -ncell_l(1), ncell_l(1)
            write(pstr_out, fmt = '(2i8)') i1,i2
            call petsc_print_master()
            j1 = 0
            j2 = j1 + nmu_l

            call petsc_add_sub_B_to_A(p_l00(i1, i2), p_c00(i1, i2), j1, j1, &
           &  p_scal, INSERT_VALUES, PETSC_FALSE)

                              
            if ( (k_mat_mode .eq. 2) .and. (diag_dim .eq. 3) ) then
              if (.not.present(p_c10)) then
                write (errormsg, fmt='(A,i8)') "diag + fix needs p_h10_l  ", 666
                call error()
              end if
              call MatZeroEntries(p_h10_cc(i1,i2),ierr)
              call petsc_add_sub_B_to_A(p_l10(i1, i2),p_c10(i1,i2), 0, &
             &  nmu_ecc - nmu_l, p_scal, INSERT_VALUES, PETSC_TRUE)
            end if
            
          end do
        end do
      end if

      if (l_use_sigma_r) then
        pstr_out = "add electrode R to C"; call petsc_print_master()
      
        do i2 = -ncell_r(2), ncell_r(2)
          do i1 = -ncell_r(1), ncell_r(1)
          
            j1 = nmu_c - nmu_r
            j2 = j1 - nmu_r          
            
            call petsc_add_sub_B_to_A(p_r00(i1, i2), p_c00(i1, i2), j1, j1, &
           & p_scal, INSERT_VALUES, PETSC_FALSE)

          
            if ((k_mat_mode .eq. 2).and.(diag_dim .eq. 3)) then
              if (.not.present(p_c01)) then
                write (errormsg, fmt='(A,i8)') "diag + fix needs p_h01_r  ", 666
                call error()
              end if            
              call MatZeroEntries(p_h01_cc(i1,i2),ierr)
              call petsc_add_sub_B_to_A(p_r01(i1, i2), p_c01(i1,i2), nmu_ecc - nmu_r, &
             & 0, p_scal, INSERT_VALUES, PETSC_TRUE)
            end if
            
          end do
        end do
      end if
      
    end subroutine fix_mat
    
end module fix_mat_mod
