module scf_mod
  implicit none
contains

  subroutine scf()
#include <petsc/finclude/petsc.h>
    use petscmat
    use petsc_mod
    use petsc_wrapper
    use kinds
    use globals
    use integrator_mod
    use ft_mod
    use kmat_mod
    use timing
    use k_on_demand_mod
    use j_dump_mod
    use error_handler


    character(1),parameter :: lr_str(2) = [ "L", "R" ]
    real(dp) :: deq

    real(dp) :: err_eq_int
    complex(dp) :: z, znorm, zdos
    integer :: ierr, ik, ifp, nmu_save, imu1, imu2, iat1, iat2, intne, i1, i2, ii(1), jj(1)
    integer(8) :: counti, count_rate, countf

    real(dp) :: de, ee, counterstart, counterend
    character(strln) :: subfile, strk
    

    PetscErrorCode :: p_ierr
    Mat :: f, a, b, p_mat1, p_mat2, p_mat3, p_Sm1, p_GVSm1
    Mat, allocatable :: p_d_tmp1(:), p_d_tmp2(:), p_pnl(:,:,:), p_Dneq_dc(:,:), p_d_dc_k(:)
    PetscScalar :: fac, v(1)
    PetscScalar, allocatable :: dc_weights_k(:), dc_weights(:, :, :)
    PetscReal :: norm
    PetscViewer :: viewer
!~       call PetscLogDefaultBegin(ierr)
    write (pstr_out, fmt='(A)') "start SCF... "; call petsc_print_master()

!~     vl = 0d0
!~     vr = vb
    
    mul = ef_l + vl
    mur = ef_r + vr
    if (lr_select .eq. -1) then
      if (mul .lt. mur) then
        mu_fermi_pol = mul
        select_lr = 1
      else
        mu_fermi_pol = mur
        select_lr = 2
      end if    
    else
      if (lr_select .eq. 1) then
        mu_fermi_pol = mul
        select_lr = 1
      else if (lr_select .eq. 2) then
        mu_fermi_pol = mur
        select_lr = 2
      else 
        write (errormsg, fmt='(A,i8)') "lr_select must be 1 or 2 ", lr_select
        call error()        
      end if    
    end if

    
    write (pstr_out, fmt='(A,i3,6e16.8)') "lr, vb, vl, vr, mul, mur mu_fermi_pol", select_lr, vb, vl ,vr, mul, mur, mu_fermi_pol
    call petsc_print_master()

    nsim_inv = nmu_c

    nint_order = nint_order_el
    l_output_progress = .true.
    ldouble_contour = ldouble_contour.and.(vb.ne.0d0)

    allocate(p_dmat_error(-ncell_c(1):ncell_c(1), -ncell_c(2):ncell_c(2)) , &
   &  p_d_tmp1(1), p_d_tmp2(1))
!~     if (ldouble_contour) allocate(p_Dneq_dc(1))
    do i1 = -ncell_c(1), ncell_c(1)
      do i2 = -ncell_c(2), ncell_c(2)
        call MatDuplicate(p_dmatxy(i1, i2, 0), MAT_SHARE_NONZERO_PATTERN, &
       &  p_dmat_error(i1, i2), ierr)
      end do
    end do
! initialize necessary matricies
    if (.not. l_k_on_demand) then
      call petsc_get_densemat(p_h00k_r(1), p_grrr, mattype_surf)     ! elemental
      call petsc_get_densemat(p_h00k_l(1), p_gllr, mattype_surf)     ! elemental
      call petsc_get_densemat(p_h00k_r(1), p_sigmarr, mattype_dense) ! dense
      call petsc_get_densemat(p_h00k_r(1), p_gammar, mattype_dense)  ! dense
      call petsc_get_densemat(p_h00k_l(1), p_sigmalr, mattype_dense) ! dense
      call petsc_get_densemat(p_h00k_l(1), p_gammal, mattype_dense)  ! dense

      call MatDuplicate(p_h00k_cc(1), MAT_SHARE_NONZERO_PATTERN, p_tmpcc1, ierr)
      call MatDuplicate(p_h00k_cc(1), MAT_SHARE_NONZERO_PATTERN, p_d_tmp1(1), ierr)
      call MatDuplicate(p_h00k_cc(1), MAT_SHARE_NONZERO_PATTERN, p_d_dc_k(1), ierr)
    else if (l_k_on_demand) then
    
      if (l_use_sigma_l) then
        call petsc_get_densemat(p_h00_ik_l(1), p_gllr, mattype_surf)     ! elemental
        call petsc_get_densemat(p_h00_ik_l(1), p_sigmalr, mattype_dense) ! dense
        call petsc_get_densemat(p_h00_ik_l(1), p_gammal, mattype_dense)  ! dense
      end if
      if (l_use_sigma_r) then
        call petsc_get_densemat(p_h00_ik_r(1), p_grrr, mattype_surf)     ! elemental
        call petsc_get_densemat(p_h00_ik_r(1), p_sigmarr, mattype_dense) ! dense
        call petsc_get_densemat(p_h00_ik_r(1), p_gammar, mattype_dense)  ! dense
      end if
      
      if (calc_current_density) then
      
        allocate(p_pnl(-ncell_c(1):ncell_c(1), -ncell_c(2):ncell_c(2), -1:1))
      
!~ addional matrix to hold non-local potential contribution to current density      
        do i2 = -ncell_c(2), ncell_c(2)
          do i1 = -ncell_c(1), ncell_c(1)
            call MatDuplicate(p_dmatxy(i1, i2, -1), MAT_SHARE_NONZERO_PATTERN, p_pnl(i1, i2, -1), ierr)
            call MatDuplicate(p_dmatxy(i1, i2, 0), MAT_SHARE_NONZERO_PATTERN, p_pnl(i1, i2, 0), ierr)
            call MatDuplicate(p_dmatxy(i1, i2, 1), MAT_SHARE_NONZERO_PATTERN, p_pnl(i1, i2, 1), ierr)
            call MatZeroEntries(p_pnl(i1, i2, -1), ierr)            
            call MatZeroEntries(p_pnl(i1, i2, 0), ierr)            
            call MatZeroEntries(p_pnl(i1, i2, 1), ierr)            
          end do
        end do

!~ for p_nl we need Dnl=G<Vnl*invS while for Dnl we only need to keep the nonzero terms due to phi_i(r)*Dnl_ij*phi_j(r)
!~ only have nonzeros if phi_i(r)*phi_j(r) .ne. 0 due to the finite cutoff, the intermidate matricies G<, Vnl, invS
!~ must be dense otherwise we will miss elements when doing the matrix multipilcation. I think all solvers for the noneq part
!~ of G< should be able to handle mattype_dense correctly.
!~ for the conventional contribution j_c we could keep the sparsity of the density matrix (likewise H matrix) but as p_nl and
!~ j_c are calculated from G< we will do this in one run using dense temporaty matrices.
        call petsc_get_densemat(p_h00_ik_cc(1), p_tmpcc1, mattype_dense) 
        call petsc_get_densemat(p_h00_ik_cc(1), p_d_tmp1(1), mattype_dense)         
        if (ldouble_contour) then
          call petsc_get_densemat(p_h00_ik_cc(1), p_tmpcc2, mattype_dense) 
        end if        
      else
        call MatDuplicate(p_h00_ik_cc(1), MAT_SHARE_NONZERO_PATTERN, p_tmpcc1, ierr)
        call MatDuplicate(p_h00_ik_cc(1), MAT_SHARE_NONZERO_PATTERN, p_d_tmp1(1), ierr)        
        if (ldouble_contour) then
          allocate(p_d_dc_k(1))
          call petsc_get_densemat(p_h00_ik_cc(1), p_d_dc_k(1), mattype_dense)     
                
          call MatDuplicate(p_h00_ik_cc(1), MAT_SHARE_NONZERO_PATTERN, p_tmpcc2, ierr)      
          call MatDuplicate(p_h00_ik_cc(1), MAT_SHARE_NONZERO_PATTERN, p_d_tmp2(1), ierr)           
        end if                
      end if
    end if

! use dense matrix. for debugging only ----------------------------------
!~       call petsc_get_densemat(p_h00k_cc(1),p_tmpcc1,mattype_dense)
!~       call petsc_get_densemat(p_h00k_cc(1),p_tmpcc1,mattype_dense)
!~       call petsc_get_densemat(p_h00k_cc(1),p_tmpcc2,mattype_dense)
!~       call petsc_get_densemat(p_h00k_cc(1),p_d_tmp1(1),mattype_dense)
!~       do i2=-ncell_c(2),ncell_c(2)
!~         do i1=-ncell_c(1),ncell_c(1)
!~           call MatConvert(p_dmatxy(i1,i2,0),mattype_dense,MAT_INPLACE_MATRIX,p_dmatxy(i1,i2,0),ierr)
!~           call MatConvert(p_dmatxy(i1,i2,-1),mattype_dense,MAT_INPLACE_MATRIX,p_dmatxy(i1,i2,-1),ierr)
!~           call MatConvert(p_dmatxy(i1,i2,1),mattype_dense,MAT_INPLACE_MATRIX,p_dmatxy(i1,i2,1),ierr)
!~         end do
!~       end do
!------------------------------------------------------------------------

!!~       call dump_kmat(dmatxy)
!!~       return

    if (ldouble_contour) then
!~       call MatCreateVecs(p_h00_ik_cc(1), vec_dc_weights, PETSC_NULL_vec, ierr)
      allocate(dc_weights(nmu_c, -ncell_c(1):ncell_c(1), -ncell_c(2):ncell_c(2)), &
     &  dc_weights_k(nmu_c), stat = ierr)
      if (ierr .ne. 0) then
        write (errormsg, fmt='(A,i8)') "allocation error ", ierr
        call error()
      end if    
      dc_weights = 0d0  
      dc_weights_k = 0d0  
      intne = 2
    else
      intne = 1
    end if

    do i2 = -ncell_c(2), ncell_c(2)
      do i1 = -ncell_c(1), ncell_c(1)
        if (ncell_c(3) .ne. 0) call MatZeroEntries(p_dmatxy(i1, i2, -1), ierr)
        call MatZeroEntries(p_dmatxy(i1, i2, 0), ierr)
        if (ncell_c(3) .ne. 0) call MatZeroEntries(p_dmatxy(i1, i2, 1), ierr)
      end do
    end do

    
    call init_eq_int(mul, mur, select_lr)
    write (pstr_out, fmt='(A,10e16.8)') "Test int: ", 0d0, 5d-1, 1d0, contour_x(0d0), contour_x(5d-1), contour_x(1d0)
    call petsc_print_master()
    write (pstr_out, fmt='(A,6e16.8)') "low: ", fermi(contour_x(0d0) - mul, temperature_el), & 
      fermi(contour_x(0d0) - mur, temperature_el)
    call petsc_print_master()
    write (pstr_out, fmt='(A,6e16.8)') "low f_l-f_r: ", fermi(contour_x(0d0) - mul, temperature_el) - fermi(contour_x(0d0) - mur, temperature_el)
    call petsc_print_master()
    write (pstr_out, fmt='(A,6e16.8)') "mid: ", fermi(contour_x(5d-1) - mul, temperature_el), & 
      fermi(contour_x(5d-1) - mur, temperature_el)
    call petsc_print_master()
    write (pstr_out, fmt='(A,6e16.8)') "mid f_l-f_r: ", fermi(contour_x(5d-1) - mul, temperature_el) - fermi(contour_x(5d-1) - mur, temperature_el)
    call petsc_print_master()
    write (pstr_out, fmt='(A,6e16.8)') "up: ", fermi(contour_x(1d0) - mul, temperature_el), & 
      fermi(contour_x(1d0) - mur, temperature_el)
    call petsc_print_master()
    write (pstr_out, fmt='(A,6e16.8)') "up f_l-f_r: ", fermi(contour_x(1d0) - mul, temperature_el) - fermi(contour_x(1d0) - mur, temperature_el)
    call petsc_print_master()
    
    do ik = nk_c, 1, -1
      strk = int2str(ik)
      write (pstr_out, fmt='(A)') repeat("-", 80); call petsc_print_master()
      if (l_k_on_demand) call k_on_demand(ik, 1)
      timing_surf = 0d0
      timing_matsolve = 0d0
      call system_clock(counti, count_rate)
      iik = ik

!~      if (((mul .ne. mur) .and. (vb .ne. 0d0)).and.(.not.l_no_noneq)) then
     if ( ( (abs(mul -  mur).ge.1d-4) .or. (vb .ne. 0d0) ).and.(.not.l_no_noneq)) then
        write (pstr_out, fmt='(A,i8)') "D_noneq. left for k-point ", iik
        call petsc_print_master()
        call init_neq_int(select_lr)
        error_reim_select = 0
        call integrate(get_gr_cc_neq, el, eu, p_d_tmp1(1), eps_int, err_eq_int)
          
        call fourier_back_add(p_d_tmp1, p_dmatxy(:, :, 0), kp_c(1:2, ik:ik), wkp_c(ik:ik), 1, &
          dlat_c, ncell_c(1), ncell_c(2), 0, 0)                                    
      

        if (ldouble_contour) then
          call MatZeroEntries(p_d_dc_k(1), ierr)
          call get_alpha(p_d_tmp1(1), p_d_tmp2(1), dc_weights_k)
        
          call MatAXPY(p_d_dc_k(1), p_minus1, p_d_tmp1(1), SAME_NONZERO_PATTERN, ierr)
          call MatAXPY(p_d_dc_k(1), p_minus1, p_d_tmp2(1), SAME_NONZERO_PATTERN, ierr)
        end if
        write (pstr_out, fmt='(A)') repeat("-", 4); call petsc_print_master()
      end if

      write (pstr_out, fmt='(A,i8)') "D_eq. left for k-point ", iik; call petsc_print_master()
      write (pstr_out, fmt='(A)') repeat("-", 4); call petsc_print_master()

! semi circle -------------------------------------
      call init_eq_int(mu_fermi_pol, mu_fermi_pol, select_lr) !1=left, 2=right, 3=average left and right fermi functions
      contour_select = 1 ! semi circle

      write (pstr_out, fmt='(A,6e16.8)') "Semi circle "//lr_str(select_lr)//": ", 0d0, phi_low, contour(0d0), contour(phi_low)
      call petsc_print_master()
      
      error_reim_select = 2
      
      call integrate(get_gr_cc_eq, 0d0, phi_low, p_d_tmp1(1), eps_int_contour, err_eq_int)

      call fourier_back_add(p_d_tmp1, p_dmatxy(:, :, 0), kp_c(1:2, ik:ik), wkp_c(ik:ik), 1, &
        dlat_c, ncell_c(1), ncell_c(2), 2, 0)
!---------------------------------------

! along E axis and fermi poles --------------------------
        contour_select = 2 ! parallel to x        
        write (pstr_out, fmt='(A,6e16.8)') "Along E-axis "//lr_str(select_lr)//": ", 0d0, 1d0, contour_x(0d0), contour_x(1d0)
        call petsc_print_master()
        
        call integrate(get_gr_cc_eq, 0d0, 1d0, p_d_tmp1(1), eps_int_contour, err_eq_int)
  
  
        call fourier_back_add(p_d_tmp1, p_dmatxy(:, :, 0), kp_c(1:2, ik:ik), wkp_c(ik:ik), 1, &
          dlat_c, ncell_c(1), ncell_c(2), 2, 0)
  
        write (pstr_out, fmt='(A)') repeat("-", 4); call petsc_print_master()
  
        write (pstr_out, fmt='(A,2i8)') "Fermi poles "//lr_str(select_lr)//": ", n_fermipols, iik; call petsc_print_master()
        call get_fermipole_contribution(p_d_tmp1(1), mu_fermi_pol) ! fermi poles
        
        call fourier_back_add(p_d_tmp1, p_dmatxy(:, :, 0), kp_c(1:2, ik:ik), wkp_c(ik:ik), 1, &
          dlat_c, ncell_c(1), ncell_c(2), 2, 0)
          
        if (ldouble_contour) then
          write (pstr_out, fmt='(A)') repeat("-", 4); call petsc_print_master()
                    
          call MatImaginaryPart(p_d_tmp1(1), ierr)
          call MatAXPY(p_d_dc_k(1), p_minus1, p_d_tmp1(1), SAME_NONZERO_PATTERN, ierr)
         
          call init_eq_int(mul, mul, 1) ! 3=average left and right fermi functions
          
          write (pstr_out, fmt='(A,i8)') "Fermi poles L", iik; call petsc_print_master()
          call get_fermipole_contribution(p_d_tmp1(1), mul) ! fermi poles
          call MatImaginaryPart(p_d_tmp1(1), ierr)
          call  MatAXPY(p_d_dc_k(1), p_one, p_d_tmp1(1), SAME_NONZERO_PATTERN, ierr)
          
          call init_eq_int(mul, mur, 4) ! 4=for DC, 3=average left and right fermi functions
          contour_select = 2 ! parallel to x        
          write (pstr_out, fmt='(A,6e16.8)') "Along E-axis DC-LR: ", 0d0, 1d0, contour_x(0d0), contour_x(1d0)
          write (pstr_out, fmt='(A,6e16.8)') "f_l-f_r: ", &
         &  fermi(contour_x(1d0) - mul, temperature_el) - &
         &  fermi(contour_x(1d0) - mur, temperature_el), &
         &  fermi(contour_x(0d0) - mul, temperature_el) - &
         &  fermi(contour_x(0d0) - mur, temperature_el)
          call petsc_print_master()
          
          subfile = "subdiv_contour2_DC-LR_"//trim(strk)//".dat"
          ldouble_contour = .false. ! switch off to avoid unecessary evaluation of p_d_tmp2 during integration
          call integrate(get_gr_cc_eq, 0d0, 1d0, p_d_tmp1(1), eps_int_contour, err_eq_int)
          
          ldouble_contour = .true. 
          call MatImaginaryPart(p_d_tmp1(1), ierr)
          call MatAXPY(p_d_dc_k(1), p_one, p_d_tmp1(1), SAME_NONZERO_PATTERN, ierr)
          
          call scale_Dne_dc(p_d_dc_k(1), dc_weights_k)
          
          call fourier_back_add(p_d_dc_k, p_dmatxy(:, :, 0), kp_c(1:2, ik:ik), wkp_c(ik:ik), 1, &
          dlat_c, ncell_c(1), ncell_c(2), 0, 0)
          
          write (pstr_out, fmt='(A)') repeat("-", 4); call petsc_print_master()
    

      
        end if
  
  
        write (pstr_out, fmt='(A)') repeat("-", 4); call petsc_print_master()



      call system_clock(countf)
      write (pstr_out, fmt='(A,e24.12)') "time surf     ", timing_surf; call petsc_print_master()
      write (pstr_out, fmt='(A,e24.12)') "time matsolve ", timing_matsolve; call petsc_print_master()
      write (pstr_out, fmt='(A,e24.12)') "time          ", real(countf - counti, 8)/real(count_rate, 8); call petsc_print_master()
    end do ! ik_loop

    write (pstr_out, fmt='(A,e24.12)') "time MatMatSolve(AIJ)       ", timing_matmat_solve(1); call petsc_print_master()
    write (pstr_out, fmt='(A,e24.12)') "time MatMatSolve(Elemental) ", timing_matmat_solve(2); call petsc_print_master()

    fac = 1d0/(pi*nktot)
!~     dc_weights = dc_weights / nktot
    do i2 = -ncell_c(2), ncell_c(2)
      do i1 = -ncell_c(1), ncell_c(1)
!~         if (ldouble_contour) then
!~           call MatRealPart(p_Dneq_dc(i1, i2), ierr)
!~           call MatRealPart(p_dmatxy(i1, i2, 0), ierr)
!~           call get_mat_norm(p_Dneq_dc(i1, i2), norm)
!~           write (pstr_out, fmt='(A, e24.12)') "Dneq_dc  : ", norm; call petsc_print_master()
!~           call get_mat_norm(p_dmatxy(i1, i2, 0), norm)
!~           write (pstr_out, fmt='(A, e24.12)') "D        : ", norm; call petsc_print_master()
!~           call MatAXPY(p_dmatxy(i1, i2, 0), p_one, p_Dneq_dc(i1, i2), SAME_NONZERO_PATTERN, ierr)
!~           call get_mat_norm(p_dmatxy(i1, i2, 0), norm)
!~           write (pstr_out, fmt='(A, e24.12)') "D+Dneq_dc: ", norm; call petsc_print_master()        
!~           if (l_ionode) then
!~             do imu1=1, nmu_c
!~               do imu2=1, nmu_c
!~                 tmpm(imu1,imu2)=dc_weights(imu1,i1,i2)*dc_weights(imu2,i1,i2)
!~               end do
!~             end do
!~             write(0, fmt='(2i4, 2e24.12)') i1, i2,maxval(abs(real(tmpm))),&
!~            &  maxval(abs(aimag(tmpm)))
!~             do imu1 = 1, nmu_c
!~               write(0, fmt='(i8, 2i4, 4e24.12)') imu1, i1, i2, dc_weights(imu1, imu1, i1, i2), &
!~              &   zsqrt(dc_weights(imu1, imu1, i1, i2))
!~             end do
!~           end if
!~         end if
        call MatScale(p_dmatxy(i1, i2, 0), fac, ierr)
      end do
    end do

    ii = 0
    jj = 0
    v = 0d0
    if (inode .eq. 0) call MatGetValues(p_dmatxy(0, 0, 0), 1, ii, 1, jj, v, ierr)
    write (pstr_out, fmt='(2i8,2e24.12)') ii, jj, v; call petsc_print_master()
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    if (calc_current_density) then
    
!~       do i2 = -ncell_c(2), ncell_c(2)
!~         do i1 = -ncell_c(1), ncell_c(1)
!~           call get_mat_norm(p_dmatxy(i1, i2, 0), norm)
!~           write (pstr_out, fmt='(A, 2i4, e24.12)') "norm[(G<-G<^T)(i1,i2,0)]: ", i1, i2, norm
!~           call petsc_print_master()                
!~           call MatDuplicate(p_dmatxy(i1, i2, 0),  MAT_COPY_VALUES, p_mat3, ierr)
!~           call MatImaginaryPart(p_mat3, ierr)
!~           call get_mat_norm(p_mat3, norm)
!~           call MatDestroy(p_mat3, ierr)
!~           write (pstr_out, fmt='(A, e24.12)') "imag: ", norm; call petsc_print_master()
!~           call MatRealPart(p_dmatxy(i1, i2, 0), ierr)
!~           call get_mat_norm(p_dmatxy(i1, i2, 0), norm)
!~           write (pstr_out, fmt='(A, e24.12)') "Real: ", norm; call petsc_print_master()                    
!~       end do
!~     end do
    call dump_j(p_dmatxy, "j_c", ecc_dir)
    
!~     do i2 = -ncell_c(2), ncell_c(2)
!~       do i1 = -ncell_c(1), ncell_c(1)      
!~         call get_mat_norm(p_pnl(i1, i2, 0), norm)
!~         write (pstr_out, fmt='(A, 2i4, e24.12)') "norm[(G<VnlinvS+(G<VnlinvS)^(T*))(i1,i2,0)]: ", i1, i2 , norm
!~         call petsc_print_master()                
!~         call MatDuplicate(p_pnl(i1, i2, 0),  MAT_COPY_VALUES, p_mat3, ierr)
!~         call MatImaginaryPart(p_mat3, ierr)
!~         call get_mat_norm(p_mat3, norm)
!~         call MatDestroy(p_mat3, ierr)
!~         write (pstr_out, fmt='(A, e24.12)') "imag: ", norm; call petsc_print_master()
!~         call MatRealPart(p_pnl(i1, i2, 0), ierr)
!~         call get_mat_norm(p_pnl(i1, i2, 0), norm)
!~         write (pstr_out, fmt='(A, e24.12)') "Real: ", norm; call petsc_print_master()                  
!~       end do
!~     end do
    call dump_j(p_pnl, "p_nl", ecc_dir)
          
    else
      call dump_kmat(l_diag_fixK, "K_negf_"//trim(spin_str)//"_matrix2.i00.p000000", ecc_dir)
      if (l_reaktor) then
        if (l_use_sigma_l) then
          call dump_kmat_lr("K_negf_matrix2.i00.p000000", trim(elec_l_dir)//"_reaktor", "l")
        end if
        if (l_use_sigma_r) then
          call dump_kmat_lr("K_negf_matrix2.i00.p000000", trim(elec_r_dir)//"_reaktor", "r")      
        end if
      end if
    end if

!~       call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'profile.txt',viewer,ierr)
!~       call PetscViewerPushFormat(viewer,PETSC_VIEWER_DEFAULT,ierr)
!~       call PetscLogView(viewer,ierr)
!~       call PetscViewerDestroy(viewer,ierr)



    write (pstr_out, fmt='(A)') "clean up..."; call petsc_print_master()

    call MatDestroy(p_grrr, ierr)
    call MatDestroy(p_gllr, ierr)
    call MatDestroy(p_sigmarr, ierr)
    call MatDestroy(p_gammar, ierr)
    call MatDestroy(p_sigmalr, ierr)
    call MatDestroy(p_gammal, ierr)

    do i1 = -ncell_c(1), ncell_c(1)
      do i2 = -ncell_c(2), ncell_c(2)
        call MatDestroy(p_dmat_error(i1, i2), ierr)
      end do
    end do
    

    call MatDestroy(p_tmpcc1, ierr)
    call MatDestroy(p_d_tmp1(1), ierr)
    if (ldouble_contour) then
      call MatDestroy(p_tmpcc2, ierr)
      call MatDestroy(p_d_tmp2(1), ierr)
      call MatDestroy(p_d_dc_k(1), ierr)
    end if

  end subroutine scf

end module scf_mod
