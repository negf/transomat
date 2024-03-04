module analysis
  implicit none
  
  contains
  
  subroutine bond_order()
#include <petsc/finclude/petsc.h>
    use petsc  
    use petsc_mod
    use globals
    use k_on_demand_mod
    implicit none
    
    Mat :: p_DS
    integer :: ik, iat1, iat2, imu1, imu2, ipos, jpos, ierr, iunit, iam, &
      nao_max, nn, n1, n2
    PetscScalar :: p_Bmunu, p_Bnumu, Bk_ij
    PetscScalar, allocatable :: B_ij(:,:), p_BBmunu(:,:), p_BBnumu(:,:)
    integer(8) :: counti, count_rate, countf
    
    nao_max = maxval(imu_ecc)
    allocate(B_ij(nat_ecc,nat_ecc), p_BBmunu(nao_max, nao_max), &
      p_BBnumu(nao_max, nao_max), stat = ierr)
    B_ij = 0d0
    do ik = nk_c, 1, -1
      call k_on_demand(ik, 2)
      call MatMatMult(p_k00k_cc(ik), p_s00k_cc(ik),  MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, p_DS, ierr)      
      write (pstr_out, fmt='(A, 2i8)') "k-point, atom ", ik, nat_ecc ; call petsc_print_master(.false.)
      call system_clock(counti, count_rate)
      do iat1 = 1, nat_ecc
        write (pstr_out, fmt='(A)') repeat(ACHAR(8), 8); call petsc_print_master(.false.)            
        write (pstr_out, fmt='(i8)') nat_ecc - iat1; call petsc_print_master(.false.)
        do iat2 = 1, nat_ecc
          Bk_ij = 0d0
          p_BBmunu = 0d0
          p_BBnumu = 0d0
          do imu1 = 1, imu_ecc(iat1)
            ipos = sum(imu_ecc(1:iat1-1)) + imu1 - 1
            do imu2 = 1, imu_ecc(iat2)
              jpos = sum(imu_ecc(1:iat2-1)) + imu2 - 1
              call petsc_mat_getvalue(ipos, jpos, p_BBmunu(imu1, imu2), p_DS, 1,&
                PETSC_COMM_WORLD, iam, .false.)
              call petsc_mat_getvalue(jpos, ipos, p_BBnumu(imu2, imu1), p_DS, 1,&
                PETSC_COMM_WORLD, iam, .false.)
            end do           
          end do
          n1 = imu_ecc(iat1)
          n2 = imu_ecc(iat2)
          nn = n1 * n2
          call MPI_AllReduce(MPI_IN_PLACE, p_BBmunu(1:n1,1:n2), nn, MPI_DOUBLE_COMPLEX, &
            MPI_SUM, PETSC_COMM_WORLD, ierr)  
          call MPI_AllReduce(MPI_IN_PLACE, p_BBnumu(1:n2,1:n1), nn, MPI_DOUBLE_COMPLEX, &
            MPI_SUM, PETSC_COMM_WORLD, ierr)  
          do imu1 = 1, imu_ecc(iat1)
            do imu2 = 1, imu_ecc(iat2)
              Bk_ij = Bk_ij + p_BBmunu(imu1, imu2)*p_BBnumu(imu2, imu1)
            end do           
          end do
          
          
          B_ij(iat1, iat2) = B_ij(iat1, iat2) + Bk_ij
          if (wkp_l(ik).eq.2) B_ij(iat1, iat2) = B_ij(iat1, iat2) + conjg(Bk_ij) 
        end do         
      end do
      call MatDestroy(p_DS, ierr)
      call system_clock(countf)
      write (pstr_out, fmt='(X,e24.12)') real(countf - counti, 8)/real(count_rate, 8); call petsc_print_master()
    end do
    
    B_ij = B_ij / (nktot)
    
    if (l_ionode) then
      open(newunit = iunit, file = "bond_order.dat", action = "write", status = "replace")
      do iat1 = 1, nat_ecc
        do iat2 = 1, nat_ecc
          write(unit = iunit, fmt = '(E24.12E4)', advance="no") real(B_ij(iat1, iat2),dp)
        end do
        write(iunit,fmt=*)
      end do      
      close(iunit)
    end if
    
    call MPI_BARRIER(PETSC_COMM_WORLD, ierr)
    
  end subroutine bond_order
  
end module analysis
