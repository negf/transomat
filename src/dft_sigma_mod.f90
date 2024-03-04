      !! DFT+sigma scf
!if (dftsigma) then
!  ndfts=jmu_dftsigma-imu_dftsigma+1
!  if (inode.eq.0) write(6,*) "imu_dftsigma,jmu_dftsigma",imu_dftsigma,jmu_dftsigma
!  call petsc_split_matrix(p_h00_cc(0,0),p_tmp1,imu_dftsigma,jmu_dftsigma,imu_dftsigma,jmu_dftsigma)
!  call petsc_split_matrix(p_s00_cc(0,0),p_tmp2,imu_dftsigma,jmu_dftsigma,imu_dftsigma,jmu_dftsigma)
!  call petsc_split_matrix(p_k00_cc(0,0),p_tmp3,imu_dftsigma,jmu_dftsigma,imu_dftsigma,jmu_dftsigma)
!  call MatGetType(p_tmp1,sdummy,ierr)
!  call MatGetInfo(p_tmp1, MAT_LOCAL, info,p_ierr)
!  call MatGetOwnershipRange(p_tmp1,nl1,nl2,ierr)
!  write(sdummy,fmt='(A,5i8)') "p_tmp1 "//trim(sdummy),inode,int(info(mat_info_nz_allocated)),int(info(mat_info_nz_used)),nl1,nl2
!  call petsc_print(sdummy)
!
!  call MatGetOwnershipRange(p_tmp1,nl1,nl2,ierr)
!  call MatCreateDense(PETSC_COMM_WORLD,nl2-nl1,PETSC_DECIDE,ndfts,ndfts,PETSC_NULL_SCALAR ,p_ev,ierr)
!
!
!
!  call MatCreateVecs(p_ev,p_ew1,PETSC_NULL_vec,ierr)
!  call MatCreateVecs(p_ev,p_ew2,PETSC_NULL_vec,ierr)
!
!~    !     call diag_mat(p_tmp1,p_tmp2,p_ew1,p_ev,ndfts,ierr)
!
!
!  p_scal=d_virt
!  call MatAXPY(p_tmp1,p_scal,p_tmp2, DIFFERENT_NONZERO_PATTERN,ierr)  ! h=h+d_virt*s
!
!
!
!  call MatMatMatMult(p_tmp2,p_tmp3,p_tmp2,MAT_INITIAL_MATRIX,PETSC_DEFAULT_REAL,p_tmp4,ierr) !=s*k*s
!
!
!
!  p_scal=d_occ-d_virt
!
!  call MatAXPY(p_tmp1,p_scal,p_tmp4,DIFFERENT_NONZERO_PATTERN,ierr) !h=h+(d_occ-d_virt)*s*k*s
!
!
!
!~    !     call diag_mat(p_tmp1,p_tmp2,p_ew2,p_ev,ndfts,ierr)
!
!  call MatDestroy(p_tmp4,ierr)
!
!  call petsc_split_matrix(p_h00_cc(0,0),p_tmp4,imu_dftsigma,jmu_dftsigma,imu_dftsigma,jmu_dftsigma)
!
!
!
!
!  allocate(cols(ndfts),p_vals(ndfts),stat=ierr)
!  if (ierr.ne.0) then
!    write(0,*) "allocation error init_sys ",ierr
!    stop
!  end if
!
!  call MatGetOwnershipRange(p_tmp4,nl1,nl2,ierr)
!
!  do imu1=nl1,nl2-1
!    call MatGetRow(p_tmp4,imu1,ncols,cols,p_vals,p_ierr)
!    do imu2=1,ncols
!      idxm=imu1
!      idym=cols(imu2)
!      call MatGetValues(p_tmp1,1,idxm,1,idym,p_vals(1),ierr)
!      p_scal=p_vals(1)
!      idxm=imu1+imu_dftsigma-1
!      idym=cols(imu2)+imu_dftsigma-1
!      call MatSetValue(p_h00_cc(0,0),idxm,idym,p_scal,INSERT_VALUES,ierr)
!    end do
!    call MatRestoreRow(p_tmp4,imu1,ncols,cols,p_vals,p_ierr)
!  end do
!
!  call MatAssemblyBegin(p_h00_cc(0,0),MAT_FINAL_ASSEMBLY,ierr)
!  call MatAssemblyEnd(p_h00_cc(0,0),MAT_FINAL_ASSEMBLY,ierr)
!
!
!  call MatDestroy(p_tmp1,ierr)
!  call MatDestroy(p_tmp2,ierr)
!
!  call petsc_split_matrix(p_h00_cc(0,0),p_tmp1,imu_dftsigma,jmu_dftsigma,imu_dftsigma,jmu_dftsigma)
!  call petsc_split_matrix(p_s00_cc(0,0),p_tmp2,imu_dftsigma,jmu_dftsigma,imu_dftsigma,jmu_dftsigma)
!~    !     call diag_mat(p_tmp1,p_tmp2,p_ew1,p_ev,ndfts,ierr)
!
!  call MatDestroy(p_tmp3,ierr)
!  call MatDestroy(p_tmp4,ierr)
!  call MatDestroy(p_ev,ierr)
!  call VecDestroy(p_ew1,ierr)
!  call VecDestroy(p_ew2,ierr)
!
!  deallocate(cols)
!
!end if

