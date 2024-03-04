module petsc_mod
#include <petsc/finclude/petsc.h>
  use petscmat
  use kinds
  implicit none

! some useful constant

  PetscScalar :: p_minus1 = -1d0, p_zione = dcmplx(0d0, 1d0), p_one = 1d0, p_zero = 0d0

! real space
  Mat, allocatable, target ::     p_h00_cc(:, :), p_s00_cc(:, :), p_k00_cc(:, :), &
    p_h01_cc(:, :), p_s01_cc(:, :), p_h10_cc(:, :), p_s10_cc(:, :), &
    p_h00_l(:, :), p_s00_l(:, :), p_k00_l(:, :), &
    p_h01_l(:, :), p_s01_l(:, :), p_k01_l(:, :), &
    p_h10_l(:, :), p_s10_l(:, :), p_k10_l(:, :), &
    p_h00_r(:, :), p_s00_r(:, :), p_k00_r(:, :), &
    p_h01_r(:, :), p_s01_r(:, :), p_k01_r(:, :), &
    p_h10_r(:, :), p_s10_r(:, :), p_k10_r(:, :), &
    p_vnl00_cc(:, :)

! real space for diag
!~ fortran does not have array of pointer only pointer to array,
!~ AFAIK at the moment this is the only workaround
  type Mat_pointer_array
    Mat, pointer :: p
  end type Mat_pointer_array
  type(Mat_pointer_array), allocatable :: p_h00_diag(:,:,:), p_s00_diag(:,:,:), p_k00_diag(:,:,:)

! reciprocal space
  Mat, pointer ::     p_h00k_cc(:), p_s00k_cc(:), p_k00k_cc(:), &
    p_h00k_l(:), p_s00k_l(:), p_k00k_l(:), &
    p_h01k_l(:), p_s01k_l(:), p_k01k_l(:), &
    p_h10k_l(:), p_s10k_l(:), p_k10k_l(:), &
    p_h00k_r(:), p_s00k_r(:), p_k00k_r(:), &
    p_h01k_r(:), p_s01k_r(:), p_k01k_r(:), &
    p_h10k_r(:), p_s10k_r(:), p_k10k_r(:), &
    p_vnl00k_cc(:)

! electrodes
  Mat :: p_grrr, p_gllr, p_sigmalr, p_sigmarr, p_gammal, p_gammar

! ep coupling in real space
  Mat, allocatable :: p_dh_cc(:, :, :, :)

! ep coupling in reciprocal space
  Mat, allocatable :: p_dhk_cc(:, :, :)

! phonon dynamical matrix, EVs, EWs
  Mat, pointer :: p_Kphonon00_cc(:, :), p_Kphonon00k_cc(:), p_phonon_EV(:), p_ep_lambda_k(:, :)
  Mat, target :: p_invsqrt_mass, p_mass_matrix
  Vec, allocatable :: p_phonon_EW(:)
  Mat :: p_tmp_ep_mat

! density matrix and related auxiliary matrices
! use sparse matrices following Conquest structure,
! for current density one maybe have to use the full matrices as we have to contraction with S^-1 -> S^-1*D*S^-1
  Mat, target, allocatable :: p_dmatxy(:, :, :), p_dmatxy_l(:, :), p_dmatxy_r(:, :)
!~   Mat :: p_d_tmp1,p_d_tmp2
! auxiliary matrix holding the greens function on CC, for simplicity they are dens for now.
! eventually only the LL and RR part would be dens blocks and we could do the inversion only for selected entries using MUMPS.
  Mat, target ::     p_tmpcc1, p_tmpcc2, p_trace, p_invGr
  Mat, pointer :: p_p_tmp
  Mat :: p_preG, p_G 

! targets for k space matricies when using k-on-demand
  Mat, target ::  p_h00_ik_cc(1), p_s00_ik_cc(1), p_h00_ik_l(1), p_s00_ik_l(1), p_h10_ik_l(1), &
    p_s10_ik_l(1), p_h01_ik_l(1), p_s01_ik_l(1), p_h00_ik_r(1), p_s00_ik_r(1), p_h10_ik_r(1), &
    p_s10_ik_r(1), p_h01_ik_r(1), p_s01_ik_r(1), p_k00_ik_cc(1), p_vnl00_ik_cc(1)

  Mat, allocatable :: p_mat_d(:,:,:)

  MatReuse :: matuse

  character(strln) :: pstr_out

contains

  subroutine petsc_get_alloc_preG(A, preG, nl_block, nr_block, l_transpose_in)
#include <petsc/finclude/petsc.h>
    use petscmat
    use globals, only: l_ionode, inode
    use error_handler
    implicit none

    integer, intent(in) :: nl_block, nr_block
    logical , optional :: l_transpose_in
    Mat :: A, preG

    Mat :: preA
    integer :: ierr, irow, nrow, ncol, nzcols, nl1, nl2, nlc1, nlc2, ii(1), jj
    integer, allocatable :: cols(:)
    logical :: l_transpose, ll
    PetscScalar :: pp(1)

    l_transpose = .false.
    if (present(l_transpose_in)) l_transpose = l_transpose_in

    call MatGetSize(A, nrow, ncol, ierr)
    call MatGetOwnershipRange(A, nl1, nl2, ierr)
    call MatGetOwnershipRangeColumn(A, nlc1, nlc2, ierr)

    call MatCreate(PETSC_COMM_WORLD, preG, ierr)
    if (l_transpose) then
      call MatSetSizes(preG, PETSC_DECIDE, PETSC_DECIDE, ncol, nrow,  ierr)
    else
      call MatSetSizes(preG, nl2 - nl1, nlc2 - nlc1, nrow, ncol, ierr)
    end if
    call MatSetType(preG, MATPREALLOCATOR, ierr)
    call MatSetUp(preG, ierr)

!~     call MatXAIJSetPreallocation(preG,1,1,PETSC_NULL_INTEGER,1,PETSC_NULL_INTEGER,ierr) ! workaround to use MatGetOwnerShipRange on b to allocate specific blocks -> petsc_alloc_mat_block

    allocate (cols(ncol), stat=ierr)
    if (ierr .ne. 0) then
      write (errormsg, *) "allocation error ", ierr
      call error()
    end if

    do irow = nl1, nl2 - 1
      call MatGetRow(A, irow, nzcols, cols, PETSC_NULL_SCALAR, ierr)
      ii = irow
      if (l_transpose) then
        call MatSetValues(preG, nzcols, cols(1:nzcols), 1, ii, PETSC_NULL_SCALAR, &
          INSERT_VALUES, ierr)
      else
        call MatSetValues(preG, 1, ii, nzcols, cols(1:nzcols), p_zero, &
          INSERT_VALUES, ierr)
      end if
      call MatRestoreRow(A, irow, nzcols, cols, PETSC_NULL_SCALAR, ierr)
    end do

!~ add dense blocks to G due to dense selfenergy blocks (inverse of surface self-engery)
!~ left
    do irow = nl1, nl2 - 1     
     
      ii = irow
      if (irow.gt.nl_block - 1) cycle

      do jj = 1, nl_block        
        cols(jj) = jj - 1
      end do        

      if (l_transpose) then
        call MatSetValues(preG, nl_block, cols(1:nl_block), 1, ii, PETSC_NULL_SCALAR, &
          INSERT_VALUES, ierr)
      else
        call MatSetValues(preG, 1, ii, nl_block, cols(1:nl_block), PETSC_NULL_SCALAR, &
          INSERT_VALUES, ierr)
      end if
      
    end do
!~ right
    do irow = nl1, nl2 - 1     
     
      ii = irow
      if (irow.lt.(nrow - nr_block + 1) - 1) cycle
      
      do jj = 1, nrow - nr_block + 1
        cols(jj) = (nrow - nr_block + jj) - 1
      end do      
      if (l_transpose) then
        call MatSetValues(preG, nr_block, cols(1:nr_block), 1, ii, PETSC_NULL_SCALAR, &
          INSERT_VALUES, ierr)
      else
        call MatSetValues(preG, 1, ii, nr_block, cols(1:nr_block), PETSC_NULL_SCALAR, &
          INSERT_VALUES, ierr)
      end if
      
    end do


    
!~ It seems off-proc entries cause suddenly a  seg-fault when using preallocator.
!~ Not sure if this was an intentional update in PETSC or a bug. Eitherway
!~ above the corresponding entries are set at the correct processes.


!~     if (l_ionode) then

!~       if (nl_block.ge.1) then
!~         do irow = 1, nl_block
!~           cols(irow) = irow - 1
!~         end do
!~         do irow = 0, nl_block - 1
!~           ii = irow
!~           call MatSetValues(preG, 1, ii, nl_block, cols(1:nl_block), PETSC_NULL_SCALAR, INSERT_VALUES, ierr)
!~         end do
!~       end if

!~       if (nr_block.ge.1) then
!~         do irow = 1, nrow - nr_block + 1
!~           cols(irow) = (nrow - nr_block + irow) - 1
!~         end do
!~         do irow = (nrow - nr_block + 1) - 1, nrow - 1
!~           ii = irow
!~           call MatSetValues(preG, 1, ii, nr_block, cols(1:nr_block), PETSC_NULL_SCALAR, INSERT_VALUES, ierr)
!~         end do
!~       end if

!~     end if

    call MatAssemblyBegin(preG, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(preG, MAT_FINAL_ASSEMBLY, ierr)
    
  end subroutine petsc_get_alloc_preG

! get a single value from a vector and broadcast it to all other processes
! takes C numbering i.e. start at 0.
  subroutine petsc_vec_getvalue(ipos, p_val, p_vec, l_bcast_result)
    use globals, only: inode
    use error_handler
    implicit none

    PetscScalar :: p_val
    integer, intent(in) :: ipos
    Vec, intent(in) :: p_vec
    logical, optional :: l_bcast_result

    PetscScalar :: pp(1)
    integer :: iroot, nl1, nl2, ii(1), ierr, nrow
    logical :: l_bcast

    call VecGetSize(p_vec, nrow, ierr)
    if ((ipos .lt. 0) .or. (ipos .gt. (nrow - 1))) then
      write (errormsg, fmt='(A,2i8)') "ipos out of bounds in petsc_vec_getvalue", &
        ipos, nrow
      call error()
    end if

    l_bcast = .true.
    if (present(l_bcast_result)) l_bcast = l_bcast_result

    ii = ipos

    call VecGetOwnershipRange(p_vec, nl1, nl2, ierr)
    iroot = -1
    if ((ipos .ge. nl1) .and. (ipos .le. nl2 - 1)) then
      call VecGetValues(p_vec, 1, ii, pp, ierr)
      p_val = pp(1)
      iroot = inode
    end if

    if (l_bcast) then
      call MPI_AllReduce(MPI_IN_PLACE, iroot, 1, MPI_INTEGER, MPI_MAX, PETSC_COMM_WORLD, ierr)  ! is this actually a clever way to distribute iroot? looks stupid.
      call MPI_Bcast(p_val, 1, MPI_DOUBLE_COMPLEX, iroot, PETSC_COMM_WORLD, ierr)
    end if

  end subroutine petsc_vec_getvalue

  subroutine petsc_mat_getvalue(ipos, jpos, p_val, p_mat, commsel, icomm, iroot_out,&
    lcast_in)
    use globals, only: inode, inode_sub
    use error_handler
    implicit none

    PetscScalar :: p_val
    integer, intent(in) :: ipos, jpos, commsel, icomm
    integer, optional :: iroot_out
    logical, optional :: lcast_in
    Mat, intent(in) :: p_mat

    PetscScalar :: pp(1)
    integer :: iroot, nl1, nl2, ii(1), jj(1), ierr, nlc1, nlc2, &
      nrow, ncol
    logical :: lcast

    lcast = .true.
    if (present(lcast_in)) lcast = lcast_in

    call MatGetSize(p_mat, nrow, ncol, ierr)
    if ((ipos .lt. 0) .or. (ipos .gt. (nrow - 1))) then
      write (errormsg, fmt='(A,2i8)') "ipos out of bounds in petsc_mat_getvalue", &
        ipos, nrow
      call error()
    end if
    if ((jpos .lt. 0) .or. (jpos .gt. (ncol - 1))) then
      write (errormsg, fmt='(A,2i8)') "jpos out of bounds in petsc_mat_getvalue", &
        jpos, ncol
      call error()
    end if

    ii = ipos
    jj = jpos

    call MatGetOwnershipRange(p_mat, nl1, nl2, ierr)
    call MatGetOwnershipRangeColumn(p_mat, nlc1, nlc2, ierr)

    iroot = -1
    p_val = 0d0
    if (present(iroot_out)) iroot_out = iroot
    if ((ipos .ge. nl1) .and. (ipos .le. nl2 - 1)) then
      if (commsel .eq. 1) then
        iroot = inode
      else if (commsel .eq. 2) then
        iroot = inode_sub
      end if
      if (present(iroot_out)) iroot_out = inode
      p_val = 0d0
      call MatGetValues(p_mat, 1, ii, 1, jj, pp, ierr)
      p_val = pp(1)
    end if

    if (.not. lcast) return
    call MPI_AllReduce(MPI_IN_PLACE, iroot, 1, MPI_INTEGER, MPI_MAX, icomm, ierr)  ! is this actually a clever way to distribute iroot? looks stupid.
    call MPI_Bcast(p_val, 1, MPI_DOUBLE_COMPLEX, iroot, icomm, ierr)

  end subroutine petsc_mat_getvalue

  subroutine rows_to_procs(nrow_pproc, nrow)

    use globals, only: nprocs, inode
    implicit none

    integer :: nrow_pproc(:, :)
    integer :: nrow

    integer :: iproc, nrow_local

    nrow_local = nrow / nprocs
    do iproc = 0, nprocs - 1
      nrow_pproc(1, iproc + 1) = nrow_local * iproc + 1
      nrow_pproc(2, iproc + 1) = nrow_pproc(1, iproc + 1) + nrow_local - 1
    end do
    nrow_pproc(2, nprocs) = nrow

  end subroutine rows_to_procs

  subroutine nzs_to_procs(nrow_pproc, nzrow)
    use globals, only: nprocs
    implicit none

    integer, allocatable :: nrow_pproc(:, :), nzrow(:)

    integer :: ipos, i, rsum_low, overhead, rsum_high, nz_average_pproc, nrow

    nrow = size(nzrow)
    nz_average_pproc = max(500, nint(real(sum(nzrow), 8)/real(nprocs, 8)))

    ipos = 1
    nrow_pproc = 0
    nrow_pproc(1, 1) = 1
    rsum_low = 0
    i = 0
    overhead = 0
    do
      i = i + 1
      rsum_low = rsum_low + nzrow(i)
      rsum_high = rsum_low
      if (i + 1 .le. nrow) rsum_high = rsum_high + nzrow(i + 1)
      !~         if (inode.eq.0) write(6,*) i,ipos,rsum_low,rsum_high
      if ((rsum_low .ge. nz_average_pproc) .or. (rsum_high .ge. nz_average_pproc) .or. (i .eq. nrow)) then
        if (abs(nz_average_pproc - rsum_high) .lt. (abs(nz_average_pproc - rsum_low))) then
          nrow_pproc(2, ipos) = i + 1
          overhead = (-nz_average_pproc + rsum_high)
          i = i + 1
        else
          nrow_pproc(2, ipos) = i
          overhead = (-nz_average_pproc + rsum_low)
        end if
        rsum_low = overhead
        if ((i .eq. nrow) .or. (ipos .eq. nprocs)) exit
        ipos = ipos + 1

        nrow_pproc(1, ipos) = i + 1
      end if

    end do
    if (nrow_pproc(2, ipos) .ne. nrow) nrow_pproc(2, ipos) = nrow

  end subroutine nzs_to_procs

  subroutine petsc_print_master(lnewline, comm_in)
#include <petsc/finclude/petsc.h>
    use petscmat
    implicit none
    integer :: ierr
    MPI_Comm, optional  :: comm_in
    logical, optional :: lnewline
    logical :: lnewline2
    MPI_Comm :: comm
    lnewline2 = .true.
    if (present(lnewline)) lnewline2 = lnewline
    comm = PETSC_COMM_WORLD
    if (present(comm_in)) comm = comm_in
    if (lnewline2) then
      call PetscPrintf(comm, trim(pstr_out)//New_line('A'), ierr)
    else if (.not. lnewline2) then
      call PetscPrintf(comm, trim(pstr_out), ierr)
    end if
    pstr_out = ""
  end subroutine petsc_print_master

!~     subroutine petsc_one_mat(A,i1,i2)
! i1 and i2 in C ordering not fortran, i.e. 0...n-1
  subroutine petsc_one_mat(A, i1, i2,itype_in)
#include <petsc/finclude/petscmat.h>
    use globals, only: inode
    implicit none

    Mat :: A
    integer :: i1, i2
    integer, optional :: itype_in

    integer :: nl1, nl2, ierr, imu, idxm(1), idym(1), nrow, ncol, nlc1, nlc2
    PetscScalar :: v(1)
    Mat :: p_preA
    integer :: itype

    itype=1
    if (present(itype_in)) itype=itype_in


    call MatGetSize(A, nrow, ncol, ierr)
    call MatGetOwnershipRange(A, nl1, nl2, ierr)
    call MatGetOwnershipRangeColumn(A, nlc1, nlc2, ierr)

    v = p_one

    if (itype.eq.2) then
      call MatCreate(PETSC_COMM_WORLD, p_preA, ierr)
      call MatSetSizes(p_preA, nl2 - nl1, nlc2 - nlc1, nrow, ncol, ierr)
      call MatSetType(p_preA, MATPREALLOCATOR, ierr)
      call MatSetUp(p_preA, ierr)

      do imu = nl1, nl2 - 1
        if ((imu .ge. i1) .and. (imu .le. i2)) then
          idym = imu - i1
          idxm = imu
          call MatSetValues(p_preA, 1, idxm, 1, idym, v, INSERT_VALUES, ierr)
        end if
      end do
      call MatAssemblyBegin(p_preA, MAT_FINAL_ASSEMBLY, ierr)
      call MatAssemblyEnd(p_preA, MAT_FINAL_ASSEMBLY, ierr)

      call MatPreallocatorPreallocate(p_preA, PETSC_TRUE, A, ierr)
      call MatDestroy(p_preA,ierr)
    end if

    do imu = nl1, nl2 - 1
      if ((imu .ge. i1) .and. (imu .le. i2)) then
        idym = imu - i1
        idxm = imu
        call MatSetValues(A, 1, idxm, 1, idym, v, INSERT_VALUES, ierr)
      end if
    end do
    call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr)

  end subroutine petsc_one_mat

  subroutine petsc_alloc_mat_block(A, nrow1, nrow2, ncol1, ncol2)
#include <petsc/finclude/petscmat.h>
    use globals, only: inode
    use error_handler
    implicit none

    Mat :: A
    MatType ::  mattype_in
    integer :: nrow1, nrow2, ncol1, ncol2

    integer :: nrow_global, ncol_global, nrow_loc1, nrow_loc2, ncol_loc1, ncol_loc2, ierr, nrow_loc, ncol_loc
    integer :: irow, i, dnz, onz, d1, d2
    integer, allocatable :: d_nnz(:), o_nnz(:)

    call MatGetOwnershipRange(A, nrow_loc1, nrow_loc2, ierr)
    nrow_loc2 = nrow_loc2 - 1
    call MatGetOwnershipRangeColumn(A, ncol_loc1, ncol_loc2, ierr)
    ncol_loc2 = ncol_loc2 - 1
    call MatGetSize(A, nrow_global, ncol_global, ierr)

    if ((nrow_global .lt. nrow2) .or. (ncol_global .lt. ncol2)) then
      if (inode .eq. 0) then
        write (errormsg, fmt='(A,4i8)') "alloc block larger than global bounds",&
         &ncol_global, ncol2, nrow_global, nrow2
        call error()
      end if
    end if

    nrow_loc = nrow_loc2 - nrow_loc1 + 1
    allocate (d_nnz(nrow_loc), o_nnz(nrow_loc))
    d_nnz = 0
    o_nnz = 0
    d1 = 0
    d2 = 0

    if ((ncol1 .le. ncol_loc2) .and. (ncol1 .ge. ncol_loc1)) then
      d1 = ncol1
    else if ((ncol2 .gt. ncol_loc2) .and. (ncol1 .gt. ncol_loc1)) then
      d1 = 0
    else if (ncol1 .lt. ncol_loc1) then
      d1 = ncol_loc1
    end if

    if ((ncol2 .le. ncol_loc2) .and. (ncol2 .ge. ncol_loc1)) then
      d2 = ncol2
    else if (ncol2 .gt. ncol_loc2) then
      d2 = ncol_loc2
    end if

    dnz = max(d2 - d1 + 1, 0)
    onz = ncol2 - ncol1 + 1 - dnz

!~       write(0,fmt='(9i8,A,6i8)') inode,dnz,onz,max(nrow1,nrow_loc1),min(nrow2,nrow_loc2),nrow1,nrow_loc1,nrow2,nrow_loc2," XX ",d1,d2,ncol1,ncol2,ncol_loc1,ncol_loc2

    irow = max(nrow1, nrow_loc1) - nrow_loc1
    do i = max(nrow1, nrow_loc1), min(nrow2, nrow_loc2)
      irow = irow + 1
      if (dnz .ge. 1) then
        d_nnz(irow) = dnz
      end if
      if (onz .ge. 1) then
        o_nnz(irow) = onz
      end if
    end do

!~       write(0,*) inode,"d_nnz ",d_nnz
    call MatMPIAIJSetPreallocation(A, -1, d_nnz, -1, o_nnz, ierr)
!~       call MatMPIAIJSetPreallocation(A,1,PETSC_NULL_INTEGER,0,PETSC_NULL_INTEGER,ierr)
    call MatAssemblyBegin(a, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(a, MAT_FINAL_ASSEMBLY, ierr)
  end subroutine petsc_alloc_mat_block

!~     subroutine petsc_get_a_with_b_c(a,b,c,mattype_in)
!~ create Matrix A with number of rows from B and number of columns from C
!~ with the same local structure which is important for matrix products.
  subroutine petsc_get_a_with_b_c(a, b, c, mattype_in, l_prealloc_in)
#include <petsc/finclude/petsc.h>
    use petscmat
    use globals, only: nprocs
    implicit none

    Mat :: a, b, c
    character(*) ::  mattype_in
    logical, optional :: l_prealloc_in

    integer :: nlrow1_b, nlrow2_b, nlcol1_b, nlcol2_b,&
   &           nlrow1_c, nlrow2_c, nlcol1_c, nlcol2_c,&
   &           ngrow_b, ngcol_b, ngrow_c, ngcol_c, ierr,&
   &           ii(1)
    logical :: l_prealloc
    
    l_prealloc = .true.
    if (present(l_prealloc_in))  l_prealloc =  l_prealloc_in

    call MatGetSize(b, ngrow_b, ngcol_b, ierr)
    call MatGetSize(c, ngrow_c, ngcol_c, ierr)
    call MatGetOwnershipRange(b, nlrow1_b, nlrow2_b, ierr)
    call MatGetOwnershipRangeColumn(b, nlcol1_b, nlcol2_b, ierr)
    call MatGetOwnershipRange(c, nlrow1_c, nlrow2_c, ierr)
    call MatGetOwnershipRangeColumn(c, nlcol1_c, nlcol2_c, ierr)
    call MatCreate(PETSC_COMM_WORLD, a, ierr)
    call MatSetType(a, mattype_in, ierr)

    call MatSetSizes(a, nlrow2_b - nlrow1_b, nlcol2_c - nlcol1_c, ngrow_b, ngcol_c, ierr)

    if ((MatType_in .eq. MATELEMENTAL) .or. (MatType_in .eq. MATMPIDENSE)) then
      call MatSetUp(a, ierr)
    else
      if (l_prealloc) then
        if (nprocs .eq. -1) then
          call MatSeqAIJSetPreallocation(a, 0, PETSC_NULL_INTEGER, ierr) ! workaround to use MatGetOwnerShipRange on b to allocate specific blocks -> petsc_alloc_mat_block
        else
          call MatMPIAIJSetPreallocation(a, 0, PETSC_NULL_INTEGER, 0, PETSC_NULL_INTEGER, ierr) ! workaround to use MatGetOwnerShipRange on b to allocate specific blocks -> petsc_alloc_mat_block
        end if
      end if
    end if

    if (l_prealloc) then
      call MatAssemblyBegin(a, MAT_FINAL_ASSEMBLY, ierr)
      call MatAssemblyEnd(a, MAT_FINAL_ASSEMBLY, ierr)
    end if

  end subroutine petsc_get_a_with_b_c

!~     subroutine petsc_get_b_from_a_with_col(a,b,ncol,mattype_in)
!~ this can cause problems if it generate matricies with different local col structure
!~ better use petsc_get_a_with_b_c
  subroutine petsc_get_b_from_a_with_col(a, b, ncol, mattype_in)
#include <petsc/finclude/petsc.h>
    use petscmat
    implicit none

    Mat :: a, b
    MatType ::  mattype_in
    integer :: ncol

    integer :: nl1, nl2, ierr, ngr, ngc, nlc1, nlc2

    call MatGetSize(a, ngr, ngc, ierr)

    call MatGetOwnershipRange(a, nl1, nl2, ierr)
    call MatGetOwnershipRangeColumn(a, nlc1, nlc2, ierr)
    call MatCreate(PETSC_COMM_WORLD, b, ierr)
    call MatSetType(b, mattype_in, ierr)
    if (MatType_in .eq. MATELEMENTAL) then
      call MatSetSizes(b, PETSC_DECIDE, PETSC_DECIDE, ngr, ncol, ierr)
    else
      call MatSetSizes(b, nl2 - nl1, PETSC_DECIDE, ngr, ncol, ierr)
    end if
    if ((MatType_in .eq. MATELEMENTAL) .or. (MatType_in .eq. MATMPIDENSE)) then
      call MatSetUp(b, ierr)
    else
      call MatMPIAIJSetPreallocation(b, 0, PETSC_NULL_INTEGER, 0, PETSC_NULL_INTEGER, ierr) ! workaround to use MatGetOwnerShipRange on b to allocate specific blocks -> petsc_alloc_mat_block
    end if

    call MatAssemblyBegin(b, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(b, MAT_FINAL_ASSEMBLY, ierr)

  end subroutine petsc_get_b_from_a_with_col

  subroutine petsc_add_sub_B_to_A(B, A, nroff, ncoff, fac, addv, add_new)
#include <petsc/finclude/petscmat.h>
    use globals, only: inode
    use error_handler
    implicit none

    Mat :: A, B
    PetscScalar :: fac
    InsertMode :: addv
    PetscBool :: add_new

    integer :: nroff, ncoff

    integer :: nrow, ncol, nrowB, ncolB, ierr, nrowA, ncolA
    integer :: nl1, nl2, i

    PetscScalar, allocatable :: buf_row(:)
    PetscInt :: ncols, rows(1)
    PetscInt, allocatable :: cols(:)
    PetscBool :: loffproc_save
    PetscErrorCode :: p_ierr

!~       call petsc_mat_info(A,"A ",ierr)
!~       call petsc_mat_info(B,"B ",ierr)

    call MatGetSize(B, nrowB, ncolB, ierr)
    call MatGetSize(A, nrowA, ncolA, ierr)

    if (((nrowB + nroff) .gt. nrowA) .or. ((ncolB + ncoff) .gt. ncolA)) then
      write (errormsg, fmt='(A,6i8)') "sub size B + off .gt. A",&
      &nrowB, nroff, nrowA, ncolB, ncoff, ncolA
      call error()
    end if

    call MatSetOption(A, MAT_NEW_NONZERO_LOCATIONS, add_new, ierr)
    call MatSetOption(A, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE.and.(.not.add_new), ierr)
    call MatSetOption(A, MAT_IGNORE_OFF_PROC_ENTRIES, PETSC_FALSE, ierr)

    allocate (buf_row(ncolB), cols(ncolB), stat=ierr)

    call MatGetOwnershipRange(B, nl1, nl2, ierr)
    do i = nl1, nl2 - 1
      call MatGetRow(B, i, ncols, cols, buf_row, p_ierr)
      if (fac .ne. 1d0) buf_row = buf_row*fac
      rows = i + nroff ! C not fortran index
      cols = cols + ncoff
      call MatSetValues(A, 1, rows, ncols, cols, buf_row(1:ncols), addv, ierr)
      if (ierr .ne. 0) then
        write (errormsg, fmt='(A,i8)') "MatSetValues error ", ierr          
        call error()   
      end if
      call MatRestoreRow(B, i, ncols, cols, buf_row, p_ierr)

    end do

    call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr)
    if (ierr .ne. 0) then
      write (errormsg, fmt='(A,i8)') "MatAssemblyBegin error ", ierr          
      call error()   
    end if    
    call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr)
    if (ierr .ne. 0) then
      write (errormsg, fmt='(A,i8)') "MatAssemblyEnd error ", ierr          
      call error()   
    end if    

    call MatSetOption(A, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE, ierr)
    call MatSetOption(A, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE, ierr)

  end subroutine petsc_add_sub_B_to_A

! Y=a*X+Y
! possible to restrict to nonzero structure of Y or not
  subroutine petsc_aXpY(Y, X, a, laddnew, lt, lignoreoffprog_in)
#include <petsc/finclude/petscmat.h>
    use globals, only: inode
    use error_handler
    implicit none

    Mat :: Y, X
    PetscBool :: laddnew
    PetscBool, optional :: lignoreoffprog_in
    logical :: lt

    PetscScalar :: a

    integer :: nroff, ncoff

    integer :: nrow, ncol, nrowY, ncolY, ierr, nrowX, ncolX
    integer :: nl1Y, nl2Y, nl1X, nl2X, i

    PetscScalar, allocatable :: buf_row(:)
    PetscInt :: ncols, rows(1)
    PetscInt, allocatable :: cols(:)
    PetscBool :: laddnew_save, loffproc_save, lignoreoffprog
    PetscErrorCode :: p_ierr

    lignoreoffprog = PETSC_FALSE
    if (present(lignoreoffprog_in)) lignoreoffprog =  lignoreoffprog_in

    call MatGetSize(Y, nrowY, ncolY, ierr)
    call MatGetSize(X, nrowX, ncolX, ierr)

    if (((nrowY .ne. nrowX) .or. (ncolY .ne. ncolX)).and.(.not.lt)) then
      write (errormsg, fmt='(A,4i16)') "petsc_aXpY: ((nrowY.ne.nrowX).or.(ncolY.ne.ncolX))",&
      &nrowY, nrowX, ncolY, ncolX
      call error()
    else if (((ncolY .ne. nrowX) .or. (nrowY .ne. ncolX)).and.(lt)) then
      write (errormsg, fmt='(A,4i16)') "petsc_aXpY transpose: ((ncolY.ne.nrowX).or.(nrowY.ne.ncolX))",&
      &ncolY, nrowX, nrowY, ncolX
      call error()
    end if

    call MatGetOwnershipRange(X, nl1X, nl2X, ierr)
    call MatGetOwnershipRange(Y, nl1Y, nl2Y, ierr)

!~       if ((nl1X.ne.nl1Y).and.(nl2X.ne.nl2Y)) then
!~         write(0,*) "petsc_aXpY: local layout differes"
!~         write(0,*) inode,nl1x,nl2x,nl1y,nl2y
!~         stop
!~       end if

    call MatGetOption(Y, MAT_NEW_NONZERO_LOCATIONS, laddnew_save, ierr)
    call MatSetOption(Y, MAT_NEW_NONZERO_LOCATIONS, laddnew, ierr)
    call MatSetOption(Y, MAT_IGNORE_OFF_PROC_ENTRIES, lignoreoffprog, ierr)

    allocate (buf_row(ncolX), cols(ncolX), stat=ierr)

    do i = nl1X, nl2X - 1
      call MatGetRow(X, i, ncols, cols, buf_row, p_ierr)
      buf_row = buf_row*a
      rows = i
      if (lt) then
        call MatSetValues(Y, ncols, cols, 1, rows, buf_row, ADD_VALUES, ierr)
      else
        call MatSetValues(Y, 1, rows, ncols, cols, buf_row, ADD_VALUES, ierr)
      end if
      call MatRestoreRow(X, i, ncols, cols, buf_row, p_ierr)
    end do
    call MatAssemblyBegin(Y, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(Y, MAT_FINAL_ASSEMBLY, ierr)
    call MatSetOption(Y, MAT_NEW_NONZERO_LOCATIONS, laddnew_save, ierr)
  end subroutine petsc_aXpY

  subroutine petsc_get_densemat(a, b, mattype_in)
#include <petsc/finclude/petsc.h>
    use petscmat
    implicit none

    Mat :: a, b
    MatType ::  mattype_in

    integer :: nl1, nl2, ierr, ngr, ngc, nlc1, nlc2, nloc, mloc

    call MatGetSize(a, ngr, ngc, ierr)

    call MatGetOwnershipRange(a, nl1, nl2, ierr)
    call MatGetOwnershipRangeColumn(a, nlc1, nlc2, ierr)
    call MatCreate(PETSC_COMM_WORLD, b, ierr)
    call MatSetType(b, mattype_in, ierr)
!~       call MatCreateDense(PETSC_COMM_WORLD,nl2-nl1,PETSC_DECIDE,ngr,ngc,PETSC_NULL_SCALAR ,b,ierr)
    if (MatType_in .eq. MATSCALAPACK) then
      nloc = PETSC_DECIDE
      mloc = PETSC_DECIDE
      call PetscSplitOwnershipEqual(PETSC_COMM_WORLD, nloc, ngr, ierr)
      call PetscSplitOwnershipEqual(PETSC_COMM_WORLD, mloc, ngc, ierr)
      call MatSetSizes(b, nloc, mloc, ngr, ngc, ierr)  !so we have to use this which now breaks AYPX. -> update the corresponding call somehow
      call MatSetFromOptions(b, ierr)
    else if (MatType_in .eq. MATELEMENTAL) then
!~         call MatSetSizes(b,nl2-nl1,nlc2-nlc1,ngr,ngc,ierr) ! this used to work and now it doen't anymore.
      call MatSetSizes(b, PETSC_DECIDE, PETSC_DECIDE, ngr, ngc, ierr)  !so we have to use this which now breaks AYPX. -> update the corresponding call somehow
    else
      call MatSetSizes(b, nl2 - nl1, nlc2 - nlc1, ngr, ngc, ierr)
    end if

    if ((MatType_in .eq. MATELEMENTAL) .or. (MatType_in .eq. MATMPIDENSE) &
      .or. (MatType_in .eq. MATSCALAPACK)) then
      call MatSetUp(b, ierr)
    end if

    call MatAssemblyBegin(b, MAT_FINAL_ASSEMBLY, ierr) !is this call necessary ?
    call MatAssemblyEnd(b, MAT_FINAL_ASSEMBLY, ierr)
  end subroutine petsc_get_densemat

!~ read binary matricies in something like the PETSc format but without ordering of the nz positions in the column vector
!  therefore PETSc MatLoad does not accept this. Moreover (in contrast to MatLoad) the files are stored NOT in big endian.

  subroutine petsc_mat_load(A, infile, initmat, initnrowproc, provide_nrow_pproc, &
    mattype, cols_loc, nzrow_loc, nzrow_out, nnz_out)
#include <petsc/finclude/petsc.h>
    use petscmat
    use kinds
    use globals, only: inode, nprocs, l_ionode
    use error_handler
    implicit none

!~      type(MPI_Request):: request

    Mat, intent(inout) :: A
    logical, intent(in) :: initmat, initnrowproc
    integer, allocatable, optional :: provide_nrow_pproc(:, :), cols_loc(:), &
      nzrow_loc(:), nzrow_out(:)
    integer, optional :: nnz_out

    character(*), intent(in) ::infile, mattype

    character(strln), parameter :: endian(2)=(/ "LITTLE_ENDIAN", "BIG_ENDIAN" /)
    PetscInt :: nz, nrow, ncol, nz_average_pproc, nrow_loc, nzdiag_loc, nzoff_loc, nnz
    PetscInt, allocatable :: nzrow(:), d_nnz(:), o_nnz(:), nzcol(:), ncsr(:)
    PetscScalar, allocatable :: buf_row(:)
    PetscScalar :: zz(1)
    PetscInt :: ii(1), jj(1)
    integer :: ipos, iunit, iclassid(1), info, fh, i, rsum_low, rsum_high, overhead, iproc, j, icol, nn, ierr
    integer, allocatable :: nrow_pproc(:, :), nz_pproc(:)
    integer  :: istatus(MPI_STATUS_SIZE)
    integer(8)  :: ioff, ioffheader
    integer :: ikind, iendian
    character(strln) :: sdummy, sprint


    iproc = inode + 1

    ikind = kind(ierr)

!~       if (inode.eq.0) write(6,fmt='(A)') "read file "//trim(infile)
    ! read file header (see PETSc MatLoad)
!~ determine endian by checking first header entry (of ikind), i.e. PETSC classid = 1211216
    if (l_ionode) then
!~       call MPI_FILE_OPEN(PETSC_COMM_WORLD, trim(infile), MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr)
!~       ioff = 0
!~       call MPI_FILE_READ_AT(fh, ioff, iclassid, 1, MPI_INTEGER, istatus, ierr)
!~       ioff = ioff + ikind
!~       call MPI_FILE_READ_AT(fh, ioff, nrow, 1, MPI_INTEGER, istatus, ierr)
!~       ioff = ioff + ikind
!~       call MPI_FILE_READ_AT(fh, ioff, ncol, 1, MPI_INTEGER, istatus, ierr)
!~       ioff = ioff + ikind
!~       call MPI_FILE_READ_AT(fh, ioff, nz, 1, MPI_INTEGER, istatus, ierr)
!~       ioff = ioff + ikind
      open(newunit = fh, file = trim(infile),access = "stream",status = "old",CONVERT = 'LITTLE_ENDIAN')
      ioff = 1
      read(unit = fh, pos = ioff)  iclassid
      iendian = 1 ! little endian
      if (iclassid(1) .ne. 1211216) then
        open(newunit = fh, file = trim(infile),access = "stream",status = "old",CONVERT = 'BIG_ENDIAN')
        read(unit = fh, pos = ioff)  iclassid
        if (iclassid(1) .ne. 1211216) then
          write (errormsg, *) "could not determine endian of PETSc files"
          call error()
        end if
        iendian = 2
      end if
      ioff = ioff + ikind
      read(unit = fh, pos = ioff) nrow
      ioff = ioff + ikind
      read(unit = fh, pos = ioff) ncol
      ioff = ioff + ikind
      read(unit = fh, pos = ioff) nz
      ioff = ioff + ikind

    end if

    call MPI_Bcast(iendian, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
    call MPI_Bcast(nrow, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
    call MPI_Bcast(ncol, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
    call MPI_Bcast(nz, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
    call MPI_Bcast(ioff, 1, MPI_INTEGER8, 0, PETSC_COMM_WORLD, ierr)
    ioffheader = ioff - 1

    call MPI_Barrier(PETSC_COMM_WORLD, ierr)
    nnz = nz
    if (present(nnz_out)) nnz_out = nnz

    allocate (nzrow(nrow), nrow_pproc(2, nprocs), nz_pproc(nprocs), stat=ierr)
    if (ierr .ne. 0) then
      write (errormsg, *) "allocation error in petsc_mat_load: &
     &nzrow,nrow_pproc,nz_pproc", ierr, nrow, nprocs
      call error()
    end if

    if (l_ionode) then
      nzrow = 0
!~       call MPI_FILE_READ_AT(fh, ioff, nzrow(1:nrow), nrow, MPI_INTEGER, istatus, ierr)
      read(unit = fh, pos = ioff) nzrow(1:nrow)
      close(fh)
    end if

    call MPI_Bcast(nzrow, nrow, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
!~       if (inode.eq.0) then
!~         write(0,*) inode," nzrow ",nrow,nzrow(1:nrow),sum(nzrow(1:nrow))
!~         do i=1,nrow
!~           write(0,fmt='(2i16)') i,nzrow(i)
!~         end do
!~       end if

    if (initnrowproc .and. present(provide_nrow_pproc)) then
      ! move this maybe to a seperate subroutine.

      ! work out how to distribute the rows on the processes. for now simple approach:
      ! start from 1 proc and try to get as close to the average number (below or above) of rows per process and
      ! add some overhead (difference to average number) to avoid stacking at the end. probably there are much more
      ! clever ways to do it.

      nz_average_pproc = max(500, nint(real(nz, 8)/real(nprocs, 8)))
      write (sdummy, fmt=*) nz_average_pproc
      call PetscPrintf(PETSC_COMM_WORLD, "nz_average_pproc="//trim(adjustl(sdummy))//NEW_LINE('A'), ierr)

      ipos = 1
      nrow_pproc = 0
      nrow_pproc(1, 1) = 1
      rsum_low = 0
      i = 0
      overhead = 0
      do
        i = i + 1
        rsum_low = rsum_low + nzrow(i)
        rsum_high = rsum_low
        if (i + 1 .le. nrow) rsum_high = rsum_high + nzrow(i + 1)
        !~          if (inode.eq.0) write(6,*) i,ipos,rsum_low,rsum_high
        if ((rsum_low .ge. nz_average_pproc) .or. (rsum_high .ge. nz_average_pproc) .or. (i .eq. nrow)) then
          if (abs(nz_average_pproc - rsum_high) .lt. (abs(nz_average_pproc - rsum_low))) then
            nrow_pproc(2, ipos) = i + 1
            overhead = (-nz_average_pproc + rsum_high)
            i = i + 1
          else
            nrow_pproc(2, ipos) = i
            overhead = (-nz_average_pproc + rsum_low)
          end if
          rsum_low = overhead
          if ((i .eq. nrow) .or. (ipos .eq. nprocs)) exit
          ipos = ipos + 1

          nrow_pproc(1, ipos) = i + 1
        end if

      end do
      if (nrow_pproc(2, ipos) .ne. nrow) nrow_pproc(2, ipos) = nrow
      provide_nrow_pproc = nrow_pproc
      if (present(nzrow_out)) then
        if (.not. allocated(nzrow_out)) allocate (nzrow_out(nrow))
        nzrow_out = nzrow
      end if
      return
    else if (initnrowproc .and. .not. present(provide_nrow_pproc)) then
      write (errormsg, fmt='(A)') "initnrowproc=.true. requiers provide_nrow_pproc"
      call error()
    else
      nrow_pproc = provide_nrow_pproc
    end if

    nn = 0
    do i = 1, nprocs
      if ((nrow_pproc(2, i) - nrow_pproc(1, i)) .eq. 0) cycle
      nn = nn + (nrow_pproc(2, i) - nrow_pproc(1, i) + 1)
      if (nn .ge. nrow) then
        nrow_pproc(2, i) = nrow
        nrow_pproc(:, i + 1:nprocs) = 0
        exit
      end if
    end do

    nz_pproc = 0
    do i = 1, nprocs
      if (nrow_pproc(1, i) .ge. 1) nz_pproc(i) = sum(nzrow(nrow_pproc(1, i):nrow_pproc(2, i)))
!~         write(6,fmt='(A,4i8)') "nrow_pproc ",i,nrow_pproc(1:2,i),nz_pproc(i)
    end do
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    allocate (nzcol(nz_pproc(iproc)), stat=ierr)
    if (ierr .ne. 0) then
      write (errormsg, fmt='(A,2i8)') "allocation error petsc_mat_load: nz_pproc ", ierr, nz_pproc(iproc)
      call error()
    end if

    nrow_loc = 0

    if (nrow_pproc(1, iproc) .ne. 0) nrow_loc = nrow_pproc(2, iproc) - nrow_pproc(1, iproc) + 1

    allocate (d_nnz(nrow_loc), o_nnz(nrow_loc), stat=ierr)
    if (ierr .ne. 0) then
      write (errormsg, *) "allocation error petsc_mat_load: d_nnz,o_nnz ", ierr, nrow_loc
      call error()
    end if
    d_nnz = 0d0
    o_nnz = 0d0

!~     call MPI_FILE_OPEN(PETSC_COMM_WORLD, trim(infile), MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr)
    open(newunit = fh, file = trim(infile), access = "stream", status = "old", &
      CONVERT = trim(endian(iendian)))
    if (nrow_loc .ge. 1) then ! only read if process holds rows
      nzcol = 0
      ioff = ioffheader + nrow*ikind + sum(nz_pproc(1:inode))*ikind + 1
      if (nz_pproc(iproc).ge.1) read(unit = fh, pos = ioff) nzcol(1:nz_pproc(iproc))
!~       call MPI_FILE_READ_AT(fh, ioff, nzcol(1:nz_pproc(iproc)), nz_pproc(iproc), MPI_INTEGER, istatus, ierr)

      nzcol = nzcol + 1 ! to fortran numbering

      icol = 0
      d_nnz = 0
      o_nnz = 0
      do i = nrow_pproc(1, iproc), nrow_pproc(2, iproc)
        nzdiag_loc = 0
        do j = 1, nzrow(i)
          icol = icol + 1
!~             if (i.eq.nrow_pproc(2,iproc)) write(0,fmt='(A,4i8)') "A ",iproc,i,icol,nzcol(icol)
          if ((nzcol(icol) .ge. nrow_pproc(1, iproc)) .and. (nzcol(icol) .le. nrow_pproc(2, iproc))) then
            nzdiag_loc = nzdiag_loc + 1
            d_nnz(i - nrow_pproc(1, iproc) + 1) = d_nnz(i - nrow_pproc(1, iproc) + 1) + 1
          else
            o_nnz(i - nrow_pproc(1, iproc) + 1) = o_nnz(i - nrow_pproc(1, iproc) + 1) + 1
          end if
        end do
      end do

    end if
    close(fh)
    call MPI_Barrier(PETSC_COMM_WORLD, ierr)
    if (initmat) then
      if (present(cols_loc) .and. present(nzrow_loc)) then
        nz = sum(nzrow(nrow_pproc(1, iproc):nrow_pproc(2, iproc)))
        allocate (cols_loc(nz), nzrow_loc(nrow_loc + 1), stat=ierr)
      end if

      sdummy = "file: "//trim(adjustl(infile))
      write (sprint, fmt=*) nrow*ncol
      sdummy = trim(sdummy)//" ndim: "//trim(adjustl(sprint))
      write (sprint, fmt=*) nnz
      sdummy = trim(sdummy)//" nnz: "//trim(adjustl(sprint))
      write (sprint, fmt='(f6.3)') real(nnz, 8)/real(nrow*ncol, 8)
      sdummy = trim(sdummy)//" ratio: "//trim(adjustl(sprint))//NEW_LINE('A')
      call PetscPrintf(PETSC_COMM_WORLD, trim(sdummy), ierr)

      ! this has to be called by ALL processes
      call MatCreate(PETSC_COMM_WORLD, A, ierr)

      if (trim(mattype) .eq. trim(MATELEMENTAL)) then
        call MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, nrow, ncol, ierr)
!~           call MatSetSizes(A,nrow_loc,nrow_loc,nrow,ncol,ierr)
        call MatSetType(A, mattype, ierr)
        call MatSetUp(A, ierr)
      else if (trim(mattype) .eq. trim(MATMPIDENSE)) then
!~           call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,nrow,ncol,ierr)
        call MatSetSizes(A, nrow_loc, nrow_loc, nrow, ncol, ierr)
        call MatSetType(A, mattype, ierr)
        call MatSetUp(A, ierr)
      else
        allocate (ncsr(nrow + 1), stat=ierr)
        if (ierr .ne. 0) then
          write (errormsg, fmt='(A,i8)') "allocation error", ierr
          call error()
        end if
        call MatSetSizes(A, nrow_loc, nrow_loc, nrow, ncol, ierr)
        call MatSetType(A, mattype, ierr)
!~           call MatSetUp(A,ierr)
        ncsr = 0
        if (nrow_loc .ge. 1) then
          do i = 2, nrow_loc
            ncsr(i) = ncsr(i - 1) + nzrow(nrow_pproc(1, iproc) - 1 + i - 1)
          end do
          ncsr(nrow_loc + 1) = ncsr(nrow_loc) + nzrow(nrow_pproc(1, iproc) - 1 + nrow_loc)
        end if
        nzcol = nzcol - 1
        allocate (buf_row(nz_pproc(iproc)))
        buf_row = 0d0

        if (nprocs .eq. -1) then
          call MatSeqAIJSetPreallocationCSR(A, ncsr(1:nrow_loc + 1), nzcol(1:nz_pproc(iproc)), buf_row(1:nz_pproc(iproc)), ierr)
        else
          call MatMPIAIJSetPreallocationCSR(A, ncsr(1:nrow_loc + 1), nzcol(1:nz_pproc(iproc)), buf_row(1:nz_pproc(iproc)), ierr)

!~              call MatMPIAIJSetPreallocation(A,1,d_nnz,1,o_nnz,ierr)
!~              call MatXAIJSetPreallocation(A,1,d_nnz,o_nnz,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,ierr)

        end if
        if (present(cols_loc) .and. present(nzrow_loc)) then
          nzrow_loc = ncsr
          cols_loc = nzcol
        end if
        nzcol = nzcol + 1
        deallocate (buf_row)

        call MatAssemblyBegin(a, MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyEnd(a, MAT_FINAL_ASSEMBLY, ierr)
!~           call MatSetOption(A,        MAT_USE_HASH_TABLE,PETSC_TRUE,ierr)
!~         call petsc_mat_info(A,"A ",ierr)
      end if

    end if

!~     call MatSetOption(A, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE, ierr)
!~     call MatSetOption(A,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_FALSE,ierr)
!~     call MatSetOption(A,MAT_NEW_NONZERO_ALLOCATION_ERR ,PETSC_FALSE,ierr)
!~       if (trim(mattype).ne.trim(MATELEMENTAL)) call MatSetOption(A,MAT_IGNORE_OFF_PROC_ENTRIES ,PETSC_TRUE,ierr) ! probably that should be ok
!~     call MatSetOption(A, MAT_IGNORE_OFF_PROC_ENTRIES, PETSC_TRUE, ierr) ! probably that should be ok

    call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr)

    if (iendian .eq. 2) then
      if (nnz .ge. 1) call petsc_mat_direct_load(A, trim(infile), ierr)
    else if (iendian .eq. 1) then
      call MPI_FILE_OPEN(PETSC_COMM_WORLD, trim(infile), MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr)

      if (nrow_loc .ge. 1) then ! only read if process holds rows
        allocate (buf_row(maxval(nzrow(nrow_pproc(1, iproc):nrow_pproc(2, iproc)))), stat=ierr)
        if (ierr .ne. 0) then
          write (errormsg, *) "allocation error ", ierr
          call error()
        end if

        icol = 0
        nzcol = nzcol - 1 ! back to C numbering

  !~         call MatScale(A,p_zero,ierr)

        do i = nrow_pproc(1, iproc), nrow_pproc(2, iproc)
          if (nzrow(i) .eq. 0) cycle
          ioff = ioffheader + nrow*ikind + sum(nz_pproc(1:nprocs))*ikind + sum(nzrow(1:i - 1))*16
          call MPI_FILE_READ_AT(fh, ioff, buf_row(1:nzrow(i)), nzrow(i), MPI_COMPLEX16, istatus, ierr)
          call MatSetValues(A, 1, i - 1, nzrow(i), nzcol(icol + 1:icol + nzrow(i)), buf_row(1:nzrow(i)), INSERT_VALUES, ierr)
          icol = icol + nzrow(i)
        end do
      end if

      ! this has to be called by ALL processes


  !~       if (trim(mattype).ne.trim(MATELEMENTAL)) call MatSetOption(A,MAT_IGNORE_OFF_PROC_ENTRIES ,PETSC_FALSE,ierr)
      call MPI_FILE_CLOSE(fh, ierr)
    end if

    call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr)

    call MatSetOption(A, MAT_IGNORE_OFF_PROC_ENTRIES, PETSC_TRUE, ierr) ! probably that should be ok

  end subroutine petsc_mat_load

  subroutine petsc_mat_info(A, what, ierr)
#include <petsc/finclude/petscmat.h>
    use kinds
    use globals, only: inode, inode_group
    implicit none

    Mat, intent(in) :: A
    PetscInt, intent(out) :: ierr
    PetscErrorCode :: p_ierr
    character(*), intent(in) :: what

    PetscInt :: nl1, nl2, m, n, nlc1, nlc2, icomm, nlocalrow, nlocalcol
    real(dp) :: info(MAT_INFO_SIZE), nzratio
    character(strln) :: stype, sdummy

    call MatGetType(A, stype, ierr)
    call MatGetInfo(A, MAT_LOCAL, info, p_ierr)
    call MatGetOwnershipRange(A, nl1, nl2, ierr)
    call MatGetOwnershipRangeColumn(A, nlc1, nlc2, ierr)
    call MatGetLocalSize(A, nlocalrow, nlocalcol, ierr)
    call MatGetSize(A, m, n, ierr)
    call PetscObjectGetComm(A, icomm, ierr)
    nzratio = info(mat_info_nz_used)/real((nlc2 - nlc1)*m, 8)

    write (sdummy, fmt='(A,i6,i8,11i10,e12.4)') trim(what)//" "//trim(stype), inode, inode_group, &
   &  int(info(mat_info_nz_allocated)), int(info(mat_info_nz_used)), nl1, nl2, nl2 - nl1, &
   &  nlc1, nlc2, nlocalrow, nlocalcol, m, n, nzratio

    call petsc_print(sdummy)

  end subroutine petsc_mat_info

  subroutine petsc_mat_nzused(A,nz)
#include <petsc/finclude/petscmat.h>
    use kinds
    use globals, only: inode
    implicit none

    Mat :: A
    integer :: nz

    PetscErrorCode :: p_ierr

    real(dp) :: info(MAT_INFO_SIZE)
    call MatGetInfo(A,  MAT_GLOBAL_SUM, info, p_ierr)

    nz = int(info(mat_info_nz_used))

  end subroutine petsc_mat_nzused

! this takes fortran limits, i.e. 1....n (not C 0....n-1)
  subroutine petsc_split_matrix(p_mat, p_sub, imu1, imu2, jmu1, jmu2, mat_reuse)
#include <petsc/finclude/petscmat.h>
    use kinds
    use globals
    use mathstuff, only: ininterval
    use error_handler

    implicit none

    Mat :: p_mat, p_sub
    integer :: imu1, imu2, jmu1, jmu2
    MatReuse :: mat_reuse

    IS :: isrow, iscol
    integer, allocatable :: irow(:), icol(:)
    PetscInt :: nl1, nl2, ierr
    PetscInt :: ione, itwo, indices(2)
    real(dp) :: info(MAT_INFO_SIZE)
    integer :: imin, n, i1, i2, ip, irow_min, irow_max, icol_min, icol_max, pirow, picol
    character(256) :: sinfo

    call MatGetOwnershipRange(p_mat, nl1, nl2, ierr)

    irow_min = 1
    irow_max = 0
    icol_min = 1
    icol_max = 0

    nl1 = nl1 + 1
    nl2 = nl2

!~           write(0,fmt='(A,5i8)') "PETSc info:",inode,int(info(mat_info_nz_allocated)),int(info(mat_info_nz_used)),nl1,nl2

    if (ininterval(nl1, nl2, imu1, imu2)) then
      irow_min = max(imu1, nl1)
      irow_max = min(imu2, nl2)
    end if
    if (ininterval(nl1, nl2, jmu1, jmu2)) then
      icol_min = max(jmu1, nl1)
      icol_max = min(jmu2, nl2)
    end if

!~       write(sinfo,fmt='(A,11i8)') trim(what),inode,imu1,imu2,jmu1,jmu2,irow_min,irow_max,icol_min,icol_max,nl1,nl2
!~       call petsc_print(sinfo)

    pirow = irow_max - irow_min + 1
    picol = icol_max - icol_min + 1 ! columns in the diagonal block of each inode

    allocate (irow(pirow), icol(picol), stat=ierr)
    if (ierr .ne. 0) then
      write (errormsg, fmt='(A,i8)') "allocation error ", ierr
      call error()
    end if
    i1 = 0
    do i2 = irow_min, irow_max
      i1 = i1 + 1
      irow(i1) = i2
    end do
    i1 = 0
    do i2 = icol_min, icol_max
      i1 = i1 + 1
      icol(i1) = i2
    end do
    irow = irow - 1
    icol = icol - 1

    call ISCreateGeneral(PETSC_COMM_WORLD, pirow, irow(1:pirow),                   &
           &                     PETSC_COPY_VALUES, isrow, ierr)

    call ISCreateGeneral(PETSC_COMM_WORLD, picol, icol(1:picol),                   &
           &                     PETSC_COPY_VALUES, iscol, ierr)

    call MatCreateSubMatrix(p_mat, isrow, iscol, mat_reuse, p_sub, ierr)

    call ISDestroy(isrow,ierr)
    call ISDestroy(iscol,ierr)

  end subroutine petsc_split_matrix

! this takes fortran limits, i.e. 1....n (not C 0....n-1)
  subroutine petsc_matrix_get_is(p_mat, imu1, imu2, jmu1, jmu2, isrow, iscol)
#include <petsc/finclude/petscmat.h>
    use kinds
    use globals
    use mathstuff, only: ininterval
    use error_handler

    implicit none

    Mat :: p_mat
    integer :: imu1, imu2, jmu1, jmu2

    IS :: isrow, iscol
    integer, allocatable :: irow(:), icol(:)
    PetscInt :: nl1, nl2, ierr
    PetscInt :: ione, itwo, indices(2)
    real(dp) :: info(MAT_INFO_SIZE)
    integer :: imin, n, i1, i2, ip, irow_min, irow_max, icol_min, icol_max, pirow, picol
    character(256) :: sinfo

    call MatGetOwnershipRange(p_mat, nl1, nl2, ierr)

    irow_min = 1
    irow_max = 0
    icol_min = 1
    icol_max = 0

    nl1 = nl1 + 1
    nl2 = nl2

!~           write(0,fmt='(A,5i8)') "PETSc info:",inode,int(info(mat_info_nz_allocated)),int(info(mat_info_nz_used)),nl1,nl2

    if (ininterval(nl1, nl2, imu1, imu2)) then
      irow_min = max(imu1, nl1)
      irow_max = min(imu2, nl2)
    end if
    if (ininterval(nl1, nl2, jmu1, jmu2)) then
      icol_min = max(jmu1, nl1)
      icol_max = min(jmu2, nl2)
    end if

!~       write(sinfo,fmt='(A,11i8)') trim(what),inode,imu1,imu2,jmu1,jmu2,irow_min,irow_max,icol_min,icol_max,nl1,nl2
!~       call petsc_print(sinfo)

    pirow = irow_max - irow_min + 1
    picol = icol_max - icol_min + 1 ! columns in the diagonal block of each inode

    allocate (irow(pirow), icol(picol), stat=ierr)
    if (ierr .ne. 0) then
      write (errormsg, fmt='(A,i8)') "allocation error ", ierr
      call error()
    end if
    i1 = 0
    do i2 = irow_min, irow_max
      i1 = i1 + 1
      irow(i1) = i2
    end do
    i1 = 0
    do i2 = icol_min, icol_max
      i1 = i1 + 1
      icol(i1) = i2
    end do
    irow = irow - 1
    icol = icol - 1

    call ISCreateGeneral(PETSC_COMM_WORLD, pirow, irow(1:pirow),                   &
           &                     PETSC_COPY_VALUES, isrow, ierr)

    call ISCreateGeneral(PETSC_COMM_WORLD, picol, icol(1:picol),                   &
           &                     PETSC_COPY_VALUES, iscol, ierr)

  end subroutine petsc_matrix_get_is

  subroutine petsc_split_matrix_virtual(p_mat, p_sub, imu1, imu2, jmu1, jmu2, what)
#include <petsc/finclude/petscmat.h>
    use kinds
    use globals
    use mathstuff, only: ininterval
    use error_handler

    implicit none

    Mat :: p_mat, p_sub
    integer :: imu1, imu2, jmu1, jmu2
    character(*) :: what

    IS :: isrow, iscol
    integer, allocatable :: irow(:), icol(:)
    PetscInt :: nl1, nl2, ierr
    PetscInt :: ione, itwo, indices(2)
    real(dp) :: info(MAT_INFO_SIZE)
    integer :: imin, n, i1, i2, ip, irow_min, irow_max, icol_min, icol_max, pirow, picol
    character(256) :: sinfo

    ione = 1
    itwo = 2

    call MatGetOwnershipRange(p_mat, nl1, nl2, ierr)

    irow_min = 1
    irow_max = 0
    icol_min = 1
    icol_max = 0

    nl1 = nl1 + 1
    nl2 = nl2

!~           write(0,fmt='(A,5i8)') "PETSc info:",inode,int(info(mat_info_nz_allocated)),int(info(mat_info_nz_used)),nl1,nl2

    if (ininterval(nl1, nl2, imu1, imu2)) then
      irow_min = max(imu1, nl1)
      irow_max = min(imu2, nl2)
    end if
    if (ininterval(nl1, nl2, jmu1, jmu2)) then
      icol_min = max(jmu1, nl1)
      icol_max = min(jmu2, nl2)
    end if

!~       write(sinfo,fmt='(A,11i8)') trim(what),inode,imu1,imu2,jmu1,jmu2,irow_min,irow_max,icol_min,icol_max,nl1,nl2
!~       call petsc_print(sinfo)

    pirow = irow_max - irow_min + 1
    picol = icol_max - icol_min + 1 ! columns in the diagonal block of each inode

    allocate (irow(pirow), icol(picol), stat=ierr)
    if (ierr .ne. 0) then
      write (errormsg, fmt='(A,i8)') "allocation error ", ierr
      call error()
    end if
    i1 = 0
    do i2 = irow_min, irow_max
      i1 = i1 + 1
      irow(i1) = i2
    end do
    i1 = 0
    do i2 = icol_min, icol_max
      i1 = i1 + 1
      icol(i1) = i2
    end do
    irow = irow - 1
    icol = icol - 1

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    call ISCreateGeneral(PETSC_COMM_WORLD, pirow, irow(1:pirow),                   &
           &                     PETSC_COPY_VALUES, isrow, ierr)

    call ISCreateGeneral(PETSC_COMM_WORLD, picol, icol(1:picol),                   &
           &                     PETSC_COPY_VALUES, iscol, ierr)

    call MatCreateSubMatrixVirtual(p_mat, isrow, iscol, p_sub, ierr)
!~       call MatGetInfo(p_sub, MAT_LOCAL, info,ierr)
!~       call MatGetOwnershipRange(p_sub,nl1,nl2,ierr)
!~       write(sinfo,fmt='(A,5i8)') trim(what),inode,int(info(mat_info_nz_allocated)),int(info(mat_info_nz_used)),nl1,nl2
    call petsc_print(sinfo)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

  end subroutine petsc_split_matrix_virtual

  subroutine petsc_split_matrix_seq(p_mat, p_sub, imu1, imu2, jmu1, jmu2, what)
#include <petsc/finclude/petscmat.h>
    use kinds
    use globals
    use mathstuff, only: ininterval
    use error_handler

    implicit none

    Mat :: p_mat, p_sub
    integer :: imu1, imu2, jmu1, jmu2
    character(*) :: what

    IS :: isrow, iscol
    integer, allocatable :: irow(:), icol(:)
    PetscInt :: nl1, nl2, ierr
    PetscInt :: ione, itwo, indices(2)
    real(dp) :: info(MAT_INFO_SIZE)
    integer :: imin, n, i1, i2, ip, irow_min, irow_max, icol_min, icol_max, pirow, picol
    character(256) :: sinfo

    ione = 1
    itwo = 2

    call MatGetOwnershipRange(p_mat, nl1, nl2, ierr)

    irow_min = 1
    irow_max = 0
    icol_min = 1
    icol_max = 0

    nl1 = nl1 + 1
    nl2 = nl2

!~           write(0,fmt='(A,5i8)') "PETSc info:",inode,int(info(mat_info_nz_allocated)),int(info(mat_info_nz_used)),nl1,nl2

    if (ininterval(nl1, nl2, imu1, imu2)) then
      irow_min = max(imu1, nl1)
      irow_max = min(imu2, nl2)
    end if
    if (ininterval(nl1, nl2, jmu1, jmu2)) then
      icol_min = max(jmu1, nl1)
      icol_max = min(jmu2, nl2)
    end if

!~       write(sinfo,fmt='(A,11i8)') trim(what),inode,imu1,imu2,jmu1,jmu2,irow_min,irow_max,icol_min,icol_max,nl1,nl2
!~       call petsc_print(sinfo)

    pirow = irow_max - irow_min + 1
    picol = icol_max - icol_min + 1 ! columns in the diagonal block of each inode

    allocate (irow(pirow), icol(picol), stat=ierr)
    if (ierr .ne. 0) then
      write (errormsg, fmt='(A,i8)') "allocation error ", ierr
      call error()
    end if
    i1 = 0
    do i2 = irow_min, irow_max
      i1 = i1 + 1
      irow(i1) = i2
    end do
    i1 = 0
    do i2 = icol_min, icol_max
      i1 = i1 + 1
      icol(i1) = i2
    end do
    irow = irow - 1
    icol = icol - 1

    call mpi_barrier(PETSC_COMM_WORLD, ierr)

    call ISCreateGeneral(PETSC_COMM_WORLD, pirow, irow(1:pirow),                   &
           &                     PETSC_COPY_VALUES, isrow, ierr)

    call ISCreateGeneral(PETSC_COMM_WORLD, picol, icol(1:picol),                   &
           &                     PETSC_COPY_VALUES, iscol, ierr)

!~       call MatMPIAIJGetLocalMat(p_mat, MAT_INITIAL_MATRIX,p_sub,ierr)
!~       call MatMPISGetLocalMatCondensed(p_mat, MAT_INITIAL_MATRIX,isrow,iscol,p_sub,ierr)
!~       call MatCreateSubMatrixVirtual(p_mat,isrow,iscol,p_sub,ierr)
!~       call MatGetInfo(p_sub, MAT_LOCAL, info,ierr)
!~       call MatGetOwnershipRange(p_sub,nl1,nl2,ierr)
!~       write(sinfo,fmt='(A,5i8)') trim(what),inode,int(info(mat_info_nz_allocated)),int(info(mat_info_nz_used)),nl1,nl2
!~       call petsc_print(sinfo)
    call mpi_barrier(PETSC_COMM_WORLD, ierr)

  end subroutine petsc_split_matrix_seq

  subroutine petsc_print(sprint)
#include <petsc/finclude/petscmat.h>
    implicit none
    character(*) :: sprint
    integer :: ierr
    call PetscSynchronizedPrintf(PETSC_COMM_WORLD, trim(sprint)//NEW_LINE('A'), ierr)
    call PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT, ierr)

  end subroutine petsc_print

  subroutine petsc_cleanup()
#include <petsc/finclude/petscmat.h>
    use petscmat
    use globals
    implicit none

    integer :: ierr, ik, i1, i2

    do i2 = -ncell_c(2), ncell_c(2)
      do i1 = -ncell_c(1), ncell_c(1)

        if (l_use_sigma_l) then
          call MatDestroy(p_k00_l(i1, i2), ierr)
          call MatDestroy(p_k10_l(i1, i2), ierr)
          call MatDestroy(p_k01_l(i1, i2), ierr)
        end if
        if (l_use_sigma_r) then
          call MatDestroy(p_k00_r(i1, i2), ierr)
          call MatDestroy(p_k10_r(i1, i2), ierr)
          call MatDestroy(p_k01_r(i1, i2), ierr)
        end if

        if (k_mat_mode .eq. 2) then
          call MatDestroy(p_h01_cc(i1, i2), ierr)
          call MatDestroy(p_s01_cc(i1, i2), ierr)
          call MatDestroy(p_h10_cc(i1, i2), ierr)
          call MatDestroy(p_s10_cc(i1, i2), ierr)
        end if


        call MatDestroy(p_dmatxy(i1, i2, -1), ierr)
        call MatDestroy(p_dmatxy(i1, i2, 0), ierr)
        call MatDestroy(p_dmatxy(i1, i2, 1), ierr)
      end do
    end do

    if (l_ep) then
      do i1 = 1, nk_c
        do i2 = 1, n_ep_modes_max
          call MatDestroy(p_ep_lambda_k(i2, i1), ierr)
        end do
        call VecDestroy(p_phonon_EW(i1), ierr)
      end do
      call MatDestroy(p_tmp_ep_mat, ierr)
    end if

!~     call MatDestroy(p_preG, ierr)

    call MatDestroy(p_h00_ik_cc(1), ierr)
    call MatDestroy(p_s00_ik_cc(1), ierr)

    if (k_mat_mode .eq. 1) then
      call MatDestroy(p_invGr, ierr)
      if (l_use_sigma_l) then
        call MatDestroy(p_h00_ik_l(1), ierr)
        call MatDestroy(p_s00_ik_l(1), ierr)
        call MatDestroy(p_h10_ik_l(1), ierr)
        call MatDestroy(p_s10_ik_l(1), ierr)
        call MatDestroy(p_h01_ik_l(1), ierr)
        call MatDestroy(p_s01_ik_l(1), ierr)
      end if
      if (l_use_sigma_r) then
        call MatDestroy(p_h00_ik_r(1), ierr)
        call MatDestroy(p_s00_ik_r(1), ierr)
        call MatDestroy(p_h10_ik_r(1), ierr)
        call MatDestroy(p_s10_ik_r(1), ierr)
        call MatDestroy(p_h01_ik_r(1), ierr)
        call MatDestroy(p_s01_ik_r(1), ierr)
      end if
    end if

    do i2 = -ncell_c(2), ncell_c(2)
      do i1 = -ncell_c(1), ncell_c(1)
        call MatDestroy(p_h00_cc(i1, i2), ierr)
        call MatDestroy(p_s00_cc(i1, i2), ierr)

        if (l_use_sigma_l) then
          call MatDestroy(p_h00_l(i1, i2), ierr)
          call MatDestroy(p_s00_l(i1, i2), ierr)
          call MatDestroy(p_h10_l(i1, i2), ierr)
          call MatDestroy(p_s10_l(i1, i2), ierr)
          call MatDestroy(p_h01_l(i1, i2), ierr)
          call MatDestroy(p_s01_l(i1, i2), ierr)
        end if

        if (l_use_sigma_r) then
          call MatDestroy(p_h00_r(i1, i2), ierr)
          call MatDestroy(p_s00_r(i1, i2), ierr)
          call MatDestroy(p_h10_r(i1, i2), ierr)
          call MatDestroy(p_s10_r(i1, i2), ierr)
          call MatDestroy(p_h01_r(i1, i2), ierr)
          call MatDestroy(p_s01_r(i1, i2), ierr)
        end if
      end do
    end do

    if (l_dftu) then
      call MatDestroy(p_k00_ik_cc(1), ierr)
    end if


  end subroutine petsc_cleanup


!~ distribute the rows among the subcomms, such that they remain on the corresponding
!~ process. That is we split
!~ 1. 11 12 13 14  sub1.1
!~ 1. 11 12 13 14  sub1.2
!~ 2. 11 12 13 14  sub2.1
!~ 3. 11 12 13 14  sub2.2

  subroutine petsc_rows_to_subcomms(A, row_range_subcomm, nrows_subcomm)
#include <petsc/finclude/petsc.h>
    use petsc
    use globals, only: inode, nprocs, l_ionode, inode_sub, psubcomm, inode_group, &
      ngroups, nodes_group
    implicit none

    Mat :: A
    integer :: row_range_subcomm(2), nrows_subcomm

    integer :: nrow_loc1, nrow_loc2, ierr


    call MatGetOwnershipRange(A, nrow_loc1, nrow_loc2, ierr)
    nrows_subcomm = nrow_loc2 - nrow_loc1
    row_range_subcomm(1) = nrow_loc1
    row_range_subcomm(2) = nrow_loc2 - 1


  end subroutine petsc_rows_to_subcomms

  subroutine petsc_cols_to_subcomms(ncols, ncols_subcomm, cols_offsets, cols_count)
    use globals, only: inode, nprocs, l_ionode, inode_sub, psubcomm, inode_group, ngroups, nodes_group
    implicit none
    integer :: ncols, ncols_subcomm(2), nup, nlow, ilow, iup, i, i1, i2, ierr
    integer :: cols_offsets(:), cols_count(:)

    nup = ceiling(real(ncols)/real(ngroups))
    nlow = nup - 1
    iup = ncols - ngroups*(nup - 1)
    ilow = ngroups - iup
!~    if (inode .eq. 0) write (0, *) nup, nlow, iup, ilow

    ncols_subcomm = -1
    if (inode_group .le. iup) then
      ncols_subcomm(1) = nup*(inode_group - 1)
      ncols_subcomm(2) = ncols_subcomm(1) + nup - 1
    else
      ncols_subcomm(1) = iup*nup + nlow*(inode_group - iup - 1)
      ncols_subcomm(2) = ncols_subcomm(1) + nlow - 1
    end if

    i1 = 0
    do i = 1, ngroups
      i1 = i1 + 1
      i2 = i1 + nodes_group(i) - 1

      call flush (6)
      if (i .le. iup) then
        cols_count(i1:i2) = nup
      else
        cols_count(i1:i2) = nlow
      end if
      call MPI_Barrier(PETSC_COMM_WORLD, ierr)
      i1 = i2
    end do
    cols_offsets = -1
    cols_offsets(1:nodes_group(1)) = 0
    i1 = nodes_group(1)
    do i = 2, ngroups
      i1 = i1 + 1
      i2 = i1 + nodes_group(i) - 1
      cols_offsets(i1:i2) = cols_offsets(i1 - 1) + cols_count(i1 - 1)
      i1 = i2
    end do
!~     do i=1,iup
!~       cols_offsets(i)=(i-1)*nup
!~       cols_count(i)=nup
!~     end do

!~     do i=1,ilow
!~       cols_offsets(i)=iup*nup+(i-1)*nlow
!~       cols_count(i)=nlow
!~     end do
!~    do i1=1,nprocs
!~      if (i1.eq.inode+1) then
!~         do i=1,ngroups
!~          write(0,fmt='(2i8,A,i8,A,5i8)') inode,inode_group," cols_offsets ",&
!~          cols_offsets(i1)," cols_count ",cols_count(inode+1),nodes_group(inode_group),ncols_subcomm,ncols_subcomm(2)-ncols_subcomm(1)+1
!~         end do
!~      end if
!~      call MPI_Barrier(PETSC_COMM_WORLD,ierr)
!~    end do
  end subroutine petsc_cols_to_subcomms
! this should only be used for DENSE and ELEMENTAL because of MatSetUp
  subroutine petsc_create_mat_on_subcomms(A, psubcomm, nrows, ncols, mattype_in)
#include <petsc/finclude/petscmat.h>
    use petscmat
    use globals, only: inode_sub
    use error_handler
    implicit none

    Mat :: A
    integer :: psubcomm, ncols, nrows
    MatType :: mattype_in

    integer :: ierr

    if (.not. ((mattype_in .eq. trim(MATMPIDENSE)) .or. (mattype_in .eq. trim(MATELEMENTAL)))) then
      write (errormsg, fmt='(A,4A)') "ERROR: only MATMPIDENSE and MATELEMENTAL are supported, not ",&
      &trim(mattype_in)
      call error()
    end if

    call MatCreate(psubcomm, A, ierr)
    call MatSetType(A, mattype_in, ierr)
    call MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, nrows, ncols, ierr)
!~     if (inode_sub.eq.0) then
!~       call MatSetSizes(A,PETSC_DECIDE,536,nrows,ncols,ierr)
!~     else
!~       call MatSetSizes(A,PETSC_DECIDE,0,nrows,ncols,ierr)
!~     end if
    call MatSetUp(A, ierr)
    call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
    call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)

  end subroutine petsc_create_mat_on_subcomms

  subroutine petsc_mat_direct_load(p_A, file_name, ioerr)

    implicit none

    Mat :: p_A
    character(*) :: file_name
    integer :: ioerr

    PetscViewer :: v_file
    integer :: ierr

    call PetscViewerCreate(PETSC_COMM_WORLD, v_file, ierr)
    call PetscViewerPushFormat(v_file, PETSC_VIEWER_NATIVE, ierr)
    call PetscViewerSetType(v_file, PETSCVIEWERBINARY, ierr)
    call PetscViewerFileSetMode(v_file, FILE_MODE_READ, ierr)
    call PetscViewerFileSetName(v_file, trim(file_name), ierr)
    call MatLoad(p_A, v_file, ioerr)
    call PetscViewerDestroy(v_file, ierr)
    
  end subroutine petsc_mat_direct_load

  subroutine petsc_mat_direct_save(p_A, file_name, ioerr, viewer_type_in, comm_in)

    implicit none

    Mat :: p_A
    character(*) :: file_name
    integer :: ioerr
    character(*), optional :: viewer_type_in
    integer, optional :: comm_in

    PetscViewer :: v_file
    integer :: ierr
    PetscViewerType, :: viewer_type
    integer :: comm

    viewer_type = PETSCVIEWERBINARY
    comm = PETSC_COMM_WORLD
    if (present(viewer_type_in)) viewer_type = trim(viewer_type_in)
    if (present(comm_in)) comm = comm_in

    call PetscViewerCreate(comm, v_file, ierr)
    call PetscViewerPushFormat(v_file, PETSC_VIEWER_NATIVE, ierr)
    call PetscViewerSetType(v_file, viewer_type, ierr)
    call PetscViewerFileSetMode(v_file, FILE_MODE_WRITE, ierr)
    call PetscViewerFileSetName(v_file, trim(file_name), ierr)
    call MatView(p_A, v_file, ioerr)
    call PetscViewerDestroy(v_file, ierr)

  end subroutine petsc_mat_direct_save

  subroutine dump_nonzero_structure(p_A, filename, comm, ionode, dump_what_in)

    use globals, only: inode, nprocs, l_ionode
    use error_handler

    implicit none
    Mat :: p_A
    integer :: comm, ionode
    integer, optional :: dump_what_in
    character(*) :: filename
    integer, allocatable ::cols(:)
    integer :: nzcol, ierr, nrow, ncol, iunit, i, irow, j, k, nrow_loc1, nrow_loc2
    PetscScalar, allocatable :: prow(:)
    integer :: dump_what

    dump_what = 0
    if (present(dump_what_in)) dump_what = dump_what_in

    call MatGetSize(p_A, nrow, ncol, ierr)
    call MatGetOwnershipRange(p_A, nrow_loc1, nrow_loc2, ierr)
    allocate (prow(ncol))

    if (inode.eq.ionode) then
      open (newunit=iunit, file=trim(filename), action="write", status="replace")
      close (iunit)
    end if

    allocate (cols(ncol), stat=ierr)
    if (ierr .ne. 0) then
      write (errormsg, fmt='(A,i8)') "allocation failure ", ierr
      call error()
    end if
    
    call MPI_Barrier(comm, ierr)

    do i = 0, nprocs - 1

      if (i .eq. inode) then
        open (newunit=iunit, file=trim(filename), action="write", status="old", position="append")
        do irow = 0, nrow - 1
          cols = -1
          if ((irow .ge. nrow_loc1) .and. (irow .lt. nrow_loc2)) then
            if (dump_what .eq. 1) then
              call MatGetRow(p_A, irow, nzcol, cols, prow, ierr)
            else
              call MatGetRow(p_A, irow, nzcol, cols, PETSC_NULL_SCALAR, ierr)
            end if
            k = 1
            do j = 0, ncol - 1
              if (cols(k) .eq. j) then
                if (dump_what .eq. 1) then
                  write (iunit, fmt='(2e24.12)', advance="no") prow(k)
                else
                  write (iunit, fmt='(i2)', advance="no") 1
                end if
                k = k + 1
              else
                if (dump_what .eq. 1) then
                  write (iunit, fmt='(2e24.12)', advance="no") cmplx(0d0,0d0)
                else
                  write (iunit, fmt='(i2)', advance="no") 0
                end if
              end if
            end do
            write (iunit, fmt=*)
            if (dump_what .eq. 1) then
              call MatRestoreRow(p_A, irow, nzcol, cols, prow, ierr)
            else
              call MatRestoreRow(p_A, irow, nzcol, cols, PETSC_NULL_SCALAR, ierr)
            end if
          end if
        end do
        close (iunit)
      end if

      call MPI_Barrier(comm, ierr)

    end do

  end subroutine dump_nonzero_structure

  function petsc_get_nnz(p_A, what)

    implicit none

      Mat :: p_A
      MatInfoType :: what

      integer :: petsc_get_nnz

      real(C_double) :: info(MAT_INFO_SIZE)
      integer :: ierr

      call MatGetInfo(p_A, what, info,ierr)
      petsc_get_nnz=int(info(mat_info_nz_used))

  end function petsc_get_nnz
  
  subroutine petsc_A_px_B(A, B, C, fac, l_conjA, l_conjB, insadd)
  use petsc
  use error_handler
  implicit none

  Mat :: A, B, C
  InsertMode :: insadd
  PetscScalar :: fac
  logical :: l_conjA, l_conjB
  PetscErrorCode :: ierr
  PetscInt :: i, j, ncolsA, ncolsB, rowStart, rowEnd, n, m
  PetscInt, allocatable :: colsA(:), colsB(:)
  PetscScalar, allocatable :: valsA(:), valsB(:), valsC(:)
  PetscScalar :: valC
  logical :: l_nzsame

  ! Get the ownership range of the matrices
  call MatGetSize(C, n, m, ierr)
  call MatGetOwnershipRange(A, rowStart, rowEnd, ierr)
  allocate(valsC(m), valsB(m), valsA(m), colsA(m), colsB(m), stat = ierr)
  ! Perform pointwise multiplication
  do i = rowStart, rowEnd - 1
  
    call MatGetRow(A, i, ncolsA, colsA, PETSC_NULL_SCALAR, ierr)
    
    call MatGetValues(A, 1, [i], ncolsA, colsA(1:ncolsA), valsA(1:ncolsA), ierr)
    call MatGetValues(B, 1, [i], ncolsA, colsA(1:ncolsA), valsB(1:ncolsA), ierr)
    
    if (l_conjA) valsA(1:ncolsA) = dconjg(valsA(1:ncolsA))
    if (l_conjB) valsB(1:ncolsA) = dconjg(valsB(1:ncolsA))
    valsC(1:ncolsA) = fac * (valsA(1:ncolsA) * valsB(1:ncolsA))
    
    call MatSetValues(C, 1, [i], ncolsA, colsA(1:ncolsA), valsC(1:ncolsA), insadd, ierr)
    
    call MatRestoreRow(A, i, ncolsA, colsA, valsA, ierr)

  end do

  ! Assemble the result matrix
  call MatAssemblyBegin(C, MAT_FINAL_ASSEMBLY, ierr)
  call MatAssemblyEnd(C, MAT_FINAL_ASSEMBLY, ierr)  
  deallocate(valsC)
    
  end subroutine petsc_A_px_B
end module petsc_mod
