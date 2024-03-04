module j_dump_mod

contains

  subroutine dump_j(p_out, outdir, outfile_in)
#include <petsc/finclude/petsc.h>
    use petscmat
    use petsc_mod
    use petsc_wrapper
    use kinds
    use globals
    use write_conquest_dump

    implicit none

    Mat, allocatable :: p_out(:,:,:)
    character(*) :: outfile_in, outdir

    integer :: i1, i2, ierr
    
    Mat :: p_mat1
    character(256) :: outfile
        
!~       outfile = adjustl(trim(outfile)//"_matrix2.i00.p000000")
!~       call write_conquest2(outfile, xyz_ecc, dlat_c, nat_ecc, tomat, inao_c, neigh_c,&
!~       nlist_c, ndim_c, atoms_c, p_out)    

!~ we can only calculate the matrices defined on the scattering region (SR), to be precise 
!~ for R_i,j,k=0 (assuming PBC in i,j) for the representation on the real space grid we also 
!~ need the contribution from R_i,j_k=-1,+1, for the density matrix we subsitute here the
!~ bulk values from the electrodes. For j_c and p_nl this is not possible as we do not know
!~ them. So we make the following assumption that the extended scattering region extSR 
!~ (effective SR used in the calculation) contains enough parts of the electrodes
!~ such that LL|extSR|RR=LL|Cll3-Cll2-Cll1-Cls-SR-CsR-Crr1-Crr2-Crr3|RR with LL=Cll3=Cll2 and
!~ RR=Crr3=Crr2 hold. And we use the coupling between Cll3-Cll2 (Crr2-Crr3) as the coupling 
!~ of the extSR to the electrodes LL-Cll3=Cll3-Cll2 (Crr3-RR=Crr2-Crr3). 
!~
!~ Moreover, as Conquest does construct the neigboring units cells directly from 000
!~ So we assume PBC and substitute once Cll3->Crr3 and once Crr3->Cll3 and write out both
!~ seperately and stich everything togther in Conquest.

    do i2 = -ncell_l(2), ncell_l(2)
      do i1 = -ncell_l(1), ncell_l(1)
      
        call MatZeroEntries(p_k10_l(i1 ,i2), ierr)
        call petsc_split_matrix(p_out(i1, i2, 0), p_mat1, &
          nmu_l  + 1 , nmu_l + nmu_l, 1, nmu_l, MAT_INITIAL_MATRIX)
        call petsc_aXpY(p_k10_l(i1, i2), p_mat1, p_one, PETSC_FALSE, .false.)
        call MatDestroy(p_mat1, ierr)
        
        call MatZeroEntries(p_k00_l(i1 ,i2), ierr)
        call petsc_split_matrix(p_out(i1, i2, 0), p_mat1, &
         1+nmu_l, nmu_l+nmu_l, 1+nmu_l, nmu_l+nmu_l,  MAT_INITIAL_MATRIX)
        call petsc_aXpY(p_k00_l(i1, i2), p_mat1, p_one, PETSC_FALSE, .false.)
        call MatDestroy(p_mat1, ierr)
        
        call MatZeroEntries(p_k01_l(i1 ,i2), ierr)
        call petsc_split_matrix(p_out(i1, i2, 0), p_mat1, &
          1, nmu_l, nmu_l  + 1 , nmu_l + nmu_l, MAT_INITIAL_MATRIX)
        call petsc_aXpY(p_k01_l(i1, i2), p_mat1, p_one, PETSC_FALSE, .false.)
        call MatDestroy(p_mat1, ierr)
        
      end do
    end do
    
    do i2 = -ncell_r(2), ncell_r(2)
      do i1 = -ncell_r(1), ncell_r(1)
      
        call MatZeroEntries(p_k01_r(i1, i2), ierr)
        call petsc_split_matrix(p_out(i1, i2, 0), p_mat1, &
         nmu_ecc - nmu_r - nmu_r + 1, nmu_ecc - nmu_r, nmu_ecc - nmu_r + 1 , nmu_ecc, &
          MAT_INITIAL_MATRIX)            
        call petsc_aXpY(p_k01_r(i1, i2), p_mat1, p_one, PETSC_FALSE, .false.) 
        call MatDestroy(p_mat1, ierr)   
      
        call MatZeroEntries(p_k00_r(i1, i2), ierr)
        call petsc_split_matrix(p_out(i1, i2, 0), p_mat1, &
         nmu_ecc - nmu_r - nmu_r + 1, nmu_ecc - nmu_r, nmu_ecc - nmu_r - nmu_r + 1 , nmu_ecc - nmu_r , &
          MAT_INITIAL_MATRIX)            
        call petsc_aXpY(p_k00_r(i1, i2), p_mat1, p_one, PETSC_FALSE, .false.) 
        call MatDestroy(p_mat1, ierr)   
        
        call MatZeroEntries(p_k10_r(i1, i2), ierr)
        call petsc_split_matrix(p_out(i1, i2, 0), p_mat1, &
          nmu_ecc - nmu_r + 1 , nmu_ecc, nmu_ecc - nmu_r - nmu_r + 1, nmu_ecc - nmu_r, &
          MAT_INITIAL_MATRIX)            
        call petsc_aXpY(p_k10_r(i1, i2), p_mat1, p_one, PETSC_FALSE, .false.) 
        call MatDestroy(p_mat1, ierr)   
                          
      end do
    end do

    do i2 = -ncell_l(2), ncell_l(2)
      do i1 = -ncell_l(1), ncell_l(1)
        call MatZeroEntries(p_out(i1, i2, -1), ierr)
        call MatZeroEntries(p_out(i1, i2,  1), ierr)
        call petsc_add_sub_B_to_A(p_k10_l(i1, i2),p_out(i1, i2, -1), 0, nmu_ecc - nmu_l, p_one, INSERT_VALUES, PETSC_FALSE)
        call petsc_add_sub_B_to_A(p_k00_l(i1, i2),p_out(i1, i2, 0), nmu_ecc - nmu_l, nmu_ecc - nmu_l, p_one, INSERT_VALUES, PETSC_FALSE)   
        call petsc_add_sub_B_to_A(p_k00_l(i1, i2),p_out(i1, i2, 0), 0, 0, p_one, INSERT_VALUES, PETSC_FALSE)   
        call petsc_add_sub_B_to_A(p_k01_l(i1, i2),p_out(i1, i2,  1), nmu_ecc - nmu_l, 0, p_one, INSERT_VALUES, PETSC_FALSE)
      end do
    end do
    
    outfile = adjustl(trim(outfile_in)//"_l_matrix2.i00.p000000")
    call write_conquest2(trim(outfile), outdir, xyz_ecc, dlat_c, nat_ecc, tomat, tomatr, &
   & inao_c, neigh_c, nlist_c, ndim_c, atoms_c, p_out)

    do i2 = -ncell_r(2), ncell_r(2)
      do i1 = -ncell_r(1), ncell_r(1)
        call MatZeroEntries(p_out(i1, i2, -1), ierr)
        call MatZeroEntries(p_out(i1, i2,  1), ierr)
        call petsc_add_sub_B_to_A(p_k01_r(i1, i2), p_out(i1, i2, 1), nmu_ecc - nmu_r, 0, p_one, INSERT_VALUES, PETSC_FALSE)
        call petsc_add_sub_B_to_A(p_k00_r(i1, i2), p_out(i1, i2, 0), 0, 0, p_one, INSERT_VALUES, PETSC_FALSE)   
        call petsc_add_sub_B_to_A(p_k00_r(i1, i2), p_out(i1, i2, 0), nmu_ecc - nmu_r, nmu_ecc - nmu_r, p_one, INSERT_VALUES, PETSC_FALSE)   
        call petsc_add_sub_B_to_A(p_k10_r(i1, i2), p_out(i1, i2, -1), 0, nmu_ecc - nmu_r, p_one, INSERT_VALUES, PETSC_FALSE)                        
      end do
    end do
  
    outfile = adjustl(trim(outfile_in)//"_r_matrix2.i00.p000000")
    call write_conquest2(trim(outfile), outdir, xyz_ecc, dlat_c, nat_ecc, tomat, tomatr, &
   &  inao_c, neigh_c, nlist_c, ndim_c, atoms_c, p_out)      
    
  end subroutine dump_j
  
end module j_dump_mod
