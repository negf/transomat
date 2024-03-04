subroutine init_control(transpfile)

  use globals
  use control
  use error_handler

  implicit none

  type(cntrl_) :: cntrl
  integer :: io
  character(*) :: transpfile

  call init_cntrl(cntrl)

  open (newunit=iunit_control, file=trim(transpfile), action="read", &
    position="rewind", FORM='FORMATTED', status="old")

  call get_cntrl_key(cntrl%ecc_dir%kname, ecc_dir)
  call get_cntrl_key(cntrl%elec_l_dir%kname, elec_l_dir)
  call get_cntrl_key(cntrl%elec_r_dir%kname, elec_r_dir)
  
  call get_cntrl_key(cntrl%use_sigma_l%kname, l_use_sigma_l)
  call get_cntrl_key(cntrl%use_sigma_r%kname, l_use_sigma_r)
  ef_c_fix = -huge(ef_c_fix)
  call get_cntrl_key(cntrl%ef_c_fix%kname, ef_c_fix, l_iscritical = .false.)

  call get_cntrl_key(cntrl%trans_file%kname, trans_file)
  call get_cntrl_key(cntrl%eta_cc%kname, eta_cc)
  call get_cntrl_key(cntrl%eta_elec%kname, eta_elec)
  call get_cntrl_key(cntrl%eps_geo%kname, eps_geo)
  call get_cntrl_key(cntrl%eps_real%kname, eps_real)
  call get_cntrl_key(cntrl%conv_dez%kname, conv_dez)
  call get_cntrl_key(cntrl%kx_dez%kname, nkx_dez)
  call get_cntrl_key(cntrl%ky_dez%kname, nky_dez)
  call get_cntrl_key(cntrl%estart%kname, estart)
  call get_cntrl_key(cntrl%eend%kname, eend)
  call get_cntrl_key(cntrl%n_energy_steps%kname, n_energy_steps)
  call get_cntrl_key(cntrl%i_energy_start%kname, i_energy_start)
  call get_cntrl_key(cntrl%i_energy_end%kname, i_energy_end)
!~   call get_cntrl_key(cntrl%bias%kname, vb)
  call get_cntrl_key(cntrl%bias_l%kname, vl)
  call get_cntrl_key(cntrl%bias_r%kname, vr)

!~   vb = vb/eh
  vl = vl / eh
  vr = vr / eh
  vb = vr - vl
  
  call get_cntrl_key(cntrl%integrator%kname, integrator)
  call get_cntrl_key(cntrl%nint_order%kname, nint_order_el)
  call get_cntrl_key(cntrl%maxsub%kname, maxsub)
  
  ldouble_contour = .false.
  call get_cntrl_key(cntrl%double_contour%kname, ldouble_contour, l_iscritical = .false.)
  
  l_no_noneq = .false.
  call get_cntrl_key(cntrl%no_noneq%kname, l_no_noneq, l_iscritical = .false.)

  call get_cntrl_key(cntrl%epsfermi%kname, epsfermi)
  call get_cntrl_key(cntrl%delta_imag%kname, delta_imag)
  call get_cntrl_key(cntrl%eps_int%kname, eps_int)
  call get_cntrl_key(cntrl%eps_int_contour%kname, eps_int_contour)
  call get_cntrl_key(cntrl%elow%kname, elow)
  lr_select = -1
  call get_cntrl_key(cntrl%lr_select%kname, lr_select, l_iscritical = .false.)
  call get_cntrl_key(cntrl%temperature%kname, temperature_el)
  call get_cntrl_key(cntrl%oneshot%kname, oneshot)
  call get_cntrl_key(cntrl%currentdensity%kname, calc_current_density)
  call get_cntrl_key(cntrl%dftsigma%kname, dftsigma)

  if (dftsigma) then
    call get_cntrl_key(cntrl%d_occ%kname, d_occ)
    call get_cntrl_key(cntrl%d_virt%kname, d_virt)
    call get_cntrl_key(cntrl%imu_dftsigma%kname, imu_dftsigma)
    call get_cntrl_key(cntrl%jmu_dftsigma%kname, jmu_dftsigma)
  end if

  l_reaktor = .false.
  call get_cntrl_key(cntrl%reaktor%kname, l_reaktor, l_iscritical = .false. )

  call get_cntrl_key(cntrl%scf_conv%kname, scf_conv)  
  
  k_mat_mode = 1
  l_diag_skip_check = .false.
  nkz_diag = 1 
  diag_dim = 3
  l_diag_fixH = .true.
  l_diag_fixK = .false.
  hartree_shift = 0d0
  call get_cntrl_key(cntrl%diag_fixH%kname, l_diag_fixH, l_iscritical = .false.)
  call get_cntrl_key(cntrl%diag_fixK%kname, l_diag_fixK, l_iscritical = .false.)
  call get_cntrl_key(cntrl%k_mat_mode%kname, k_mat_mode, l_iscritical = .false.)
  if (k_mat_mode .eq. 2) then
    call get_cntrl_key(cntrl%diag_skip_check%kname, l_diag_skip_check)
    call get_cntrl_key(cntrl%diag_fix_ef%kname, l_diag_fix_ef)
    call get_cntrl_key(cntrl%diag_fixH%kname, l_diag_fixH)
    call get_cntrl_key(cntrl%diag_fixK%kname, l_diag_fixK)
    call get_cntrl_key(cntrl%diag_dim%kname, diag_dim)
    call get_cntrl_key(cntrl%diag_nkz%kname, nkz_diag)
    call get_cntrl_key(cntrl%hartree_shift%kname, hartree_shift, l_iscritical = .false.)
  end if

  call get_cntrl_key(cntrl%calculate_ti%kname, lget_ti)

  
  call get_cntrl_key(cntrl%solver_mode%kname, solver_mode)

  nsim_rhs = 1
  ngroups = 1
  call get_cntrl_key(cntrl%nsim_rhs%kname, nsim_rhs, l_iscritical=.false.)
  call get_cntrl_key(cntrl%ngroups%kname, ngroups, l_iscritical=.false.)
  nsim_rhs = max(1, nsim_rhs)
  ngroups = max(1, ngroups)
  call get_cntrl_key(cntrl%ep_loe%kname, l_ep_loe)
  if (l_ep_loe) then
    l_ep = .true.
    call get_cntrl_key(cntrl%ep_active_iat1%kname, ep_active_atoms(1))
    call get_cntrl_key(cntrl%ep_active_iat2%kname, ep_active_atoms(2))
    call get_cntrl_key(cntrl%ep_bias_min%kname, ep_bias(1))
    call get_cntrl_key(cntrl%ep_bias_max%kname, ep_bias(2))
    call get_cntrl_key(cntrl%ep_bias_n%kname, ep_bias_n)
    call get_cntrl_key(cntrl%ep_reuse_lambda%kname, l_ep_reuse_lambda)
    call get_cntrl_key(cntrl%eta_ph_cc%kname, eta_ph_cc)
    call get_cntrl_key(cntrl%temperature_ph%kname, temperature_ph)
!~     call get_cntrl_key(cntrl%eta_ph_elec%kname,eta_ph_elec)
  end if

  call get_cntrl_key(cntrl%loadsave_gf%kname, l_loadsavegf)

  call get_cntrl_key(cntrl%dump_nzs%kname, l_dump_nzs, l_iscritical=.false.)

  l_k_on_demand = .true.
  call get_cntrl_key(cntrl%k_on_demand%kname, l_k_on_demand, l_iscritical=.false.)
  
  l_dftu = .false.
  call get_cntrl_key(cntrl%dftu%kname, l_dftu, l_iscritical=.false.)
  if (l_dftu) then
    dftu_file_c=""
    dftu_file_l=""
    dftu_file_r=""
    call get_cntrl_key(cntrl%dftu_file_c%kname, dftu_file_c, l_iscritical=.false.)
    call get_cntrl_key(cntrl%dftu_file_l%kname, dftu_file_l, l_iscritical=.false.)
    call get_cntrl_key(cntrl%dftu_file_r%kname, dftu_file_r, l_iscritical=.false.)
    call get_cntrl_key(cntrl%dftu_projector%kname, dftu_projector, l_iscritical=.true.)
    call get_cntrl_key(cntrl%dftu_sparse_projector%kname, l_dftu_sparse_projector, &
   &  l_iscritical=.true.)
  end if

  spin_str="u"
  call get_cntrl_key(cntrl%spin%kname, spin_str, l_iscritical=.false.)
  if (trim(spin_str).eq."u") then
    spin_channel=1
  else if (trim(spin_str).eq."d") then
    spin_channel=2
  else
    write (errormsg, fmt='(A,A)') "spin should be u or d not",trim(spin_str)
    call error()
  end if
  
  close (iunit_control)

  if (inode .eq. 0) write (6, fmt='(A)') repeat("-", 80)

end subroutine init_control
