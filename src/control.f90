module control

  use kinds

  implicit none

  type ikey
    character(strln) :: kname
    integer :: def
  end type ikey

  type rkey
    character(strln) :: kname
    real(dp) :: def
  end type rkey

  type lkey
    character(strln) :: kname
    logical :: def
  end type lkey

  type skey
    character(strln) :: kname
    character(strln) :: def
  end type skey

  type cntrl_

    type(skey) :: ecc_dir
    type(skey) :: elec_l_dir
    type(skey) :: elec_r_dir

    type(lkey) :: use_sigma_l
    type(lkey) :: use_sigma_r
    type(rkey) :: ef_c_fix

    type(skey) :: trans_file

    type(skey) :: cell_l_file
    type(skey) :: cell_r_file

    type(rkey) :: eta_cc
    type(rkey) :: eta_elec
    type(rkey) :: eps_geo
    type(rkey) :: eps_real
    type(rkey) :: conv_dez
    type(rkey) :: estart
    type(rkey) :: eend
    
    type(ikey) :: k_mat_mode

    type(ikey) :: kx_dez
    type(ikey) :: ky_dez
    type(ikey) :: n_energy_steps
    type(ikey) :: i_energy_start
    type(ikey) :: i_energy_end

    type(rkey) :: ef_c
    type(rkey) :: ef_l
    type(rkey) :: ef_r
!~     type(rkey) :: bias
    type(rkey) :: bias_l
    type(rkey) :: bias_r

    type(ikey) :: integrator
    type(ikey) :: nint_order
    type(ikey) :: maxsub
    type(ikey) :: lr_select
    type(lkey) :: double_contour
    type(rkey) :: epsfermi
    type(rkey) :: delta_imag
    type(rkey) :: eps_int
    type(rkey) :: eps_int_contour
    type(rkey) :: elow
    type(rkey) :: temperature

    type(lkey) :: oneshot

    type(lkey) :: currentdensity

    type(lkey) :: dftsigma
    type(rkey) :: d_occ
    type(rkey) :: d_virt
    type(ikey) :: imu_dftsigma
    type(ikey) :: jmu_dftsigma

    type(skey) :: f_file_ecc_petsc
    type(skey) :: s_file_ecc_petsc
    type(skey) :: d_file_ecc_petsc

    type(skey) :: f_file_elec_l_petsc
    type(skey) :: s_file_elec_l_petsc
    type(skey) :: d_file_elec_l_petsc

    type(skey) :: f_file_elec_r_petsc
    type(skey) :: s_file_elec_r_petsc
    type(skey) :: d_file_elec_r_petsc

    type(rkey) :: scf_conv

    type(lkey) :: calculate_ti

    type(lkey) :: ep_loe

    type(ikey) :: ep_active_iat1
    type(ikey) :: ep_active_iat2

    type(rkey) :: ep_bias_min
    type(rkey) :: ep_bias_max
    type(ikey) :: ep_bias_n

    type(lkey) :: ep_reuse_lambda

    type(rkey) :: eta_ph_cc
    type(rkey) :: eta_ph_elec
    type(rkey) :: temperature_ph

    type(ikey) :: solver_mode
    type(ikey) :: nsim_rhs
    type(ikey) :: ngroups

    type(lkey) :: loadsave_gf
    type(lkey) :: no_noneq
    
    type(lkey) :: reaktor

    type(lkey) :: k_on_demand

    type(lkey) :: dump_nzs
    
    type(lkey) :: diag_fix_ef
    type(lkey) :: diag_fixH
    type(lkey) :: diag_fixK
    type(lkey) :: diag_skip_check
    type(ikey) :: diag_dim
    type(ikey) :: diag_nkz

    type(lkey) :: dftu
    type(skey) :: dftu_file_c
    type(skey) :: dftu_file_l
    type(skey) :: dftu_file_r
    type(ikey) :: dftu_projector
    type(lkey) :: dftu_sparse_projector
    
    type(skey) :: spin
    
    type(rkey) :: hartree_shift

  end type cntrl_

  interface get_cntrl_key

    subroutine get_skey(keyname, keystr, l_iscritical)
      use kinds
      character(strln) :: keyname
      character(strln) :: keystr
      logical, optional :: l_iscritical
    end subroutine get_skey

    subroutine get_rkey(keyname, kvalue, l_iscritical)
      use kinds
      character(strln) :: keyname
      real(dp) :: kvalue
      logical, optional :: l_iscritical
    end subroutine get_rkey

    subroutine get_ikey(keyname, kvalue, l_iscritical)
      use kinds
      character(strln) :: keyname
      integer :: kvalue
      logical, optional :: l_iscritical
    end subroutine get_ikey

    subroutine get_lkey(keyname, kvalue, l_iscritical)
      use kinds
      character(strln) :: keyname
      logical :: kvalue
      logical, optional :: l_iscritical
    end subroutine get_lkey

  end interface get_cntrl_key

contains

  subroutine init_cntrl(cntrl)

    implicit none

    type(cntrl_) :: cntrl

    cntrl%ecc_dir%kname = "$ecc_dir="
    cntrl%elec_l_dir%kname = "$elec_left_dir="
    cntrl%elec_r_dir%kname = "$elec_right_dir="
    
    cntrl%use_sigma_l%kname = "$use_sigma_left="
    cntrl%use_sigma_r%kname = "$use_sigma_right="
    cntrl%ef_c_fix%kname = "$ef_c="

    cntrl%trans_file%kname = "$trans_file="

    cntrl%eta_cc%kname = "$eta_cc="

    cntrl%eps_geo%kname = "$eps_geo="
    cntrl%eps_real%kname = "$eps_real="

    cntrl%eta_elec%kname = "$eta_elec="
    cntrl%conv_dez%kname = "$conv_dez="
        
    cntrl%k_mat_mode%kname = "$k_mat_mode="

    cntrl%estart%kname = "$e_start="
    cntrl%eend%kname = "$e_end="
    cntrl%n_energy_steps%kname = "$n_steps="
    cntrl%i_energy_start%kname = "$i_start="
    cntrl%i_energy_end%kname = "$i_end="

    cntrl%kx_dez%kname = "$nkx_dez="
    cntrl%ky_dez%kname = "$nky_dez="

!~     cntrl%bias%kname = "$bias="
    cntrl%bias_l%kname = "$bias_l="
    cntrl%bias_r%kname = "$bias_r="
    cntrl%integrator%kname = "$integrator="
    cntrl%nint_order%kname = "$nint_order="
    cntrl%maxsub%kname = "$maxsub="
    cntrl%lr_select%kname = "$lr_select="
    cntrl%double_contour%kname = "$double_contour="
    cntrl%epsfermi%kname = "$epsfermi="
    cntrl%delta_imag%kname = "$delta_imag="
    cntrl%eps_int%kname = "$eps_int="
    cntrl%eps_int_contour%kname = "$eps_int_contour="
    cntrl%elow%kname = "$elow="
    cntrl%temperature%kname = "$temperature="

    cntrl%oneshot%kname = "$oneshot="

    cntrl%currentdensity%kname = "$current_density="

    cntrl%dftsigma%kname = "$dftsigma="
    cntrl%d_occ%kname = "$dftsigma_occ="
    cntrl%d_virt%kname = "$dftsigma_virt="
    cntrl%imu_dftsigma%kname = "$dftsigma_imu="
    cntrl%jmu_dftsigma%kname = "$dftsigma_jmu="

    cntrl%f_file_ecc_petsc%kname = "$f_file_ecc_petsc="
    cntrl%s_file_ecc_petsc%kname = "$s_file_ecc_petsc="
    cntrl%d_file_ecc_petsc%kname = "$d_file_ecc_petsc="

    cntrl%f_file_elec_l_petsc%kname = "$f_file_elec_l_petsc="
    cntrl%s_file_elec_l_petsc%kname = "$s_file_elec_l_petsc="
    cntrl%d_file_elec_l_petsc%kname = "$d_file_elec_l_petsc="

    cntrl%f_file_elec_r_petsc%kname = "$f_file_elec_r_petsc="
    cntrl%s_file_elec_r_petsc%kname = "$s_file_elec_r_petsc="
    cntrl%d_file_elec_r_petsc%kname = "$d_file_elec_r_petsc="

    cntrl%scf_conv%kname = "$scf_conv="

    cntrl%calculate_ti%kname = "$calculate_ti="

    cntrl%ep_loe%kname = "$ep_loe="
    cntrl%ep_active_iat1%kname = "$ep_active_iat1="
    cntrl%ep_active_iat2%kname = "$ep_active_iat2="

    cntrl%ep_bias_min%kname = "$ep_bias_min="
    cntrl%ep_bias_max%kname = "$ep_bias_max="
    cntrl%ep_bias_n%kname = "$ep_bias_n="

    cntrl%ep_reuse_lambda%kname = "$ep_reuse_lambda="

    cntrl%eta_ph_cc%kname = "$eta_ph_cc="
    cntrl%eta_ph_elec%kname = "$eta_ph_elec="
    cntrl%temperature_ph%kname = "$temperature_ph="

    cntrl%solver_mode%kname = "$solver_mode="
    cntrl%nsim_rhs%kname = "$nsim_rhs="
    cntrl%ngroups%kname = "$ngroups="

    cntrl%loadsave_gf%kname = "$reuse_surface_gf="
    
    cntrl%no_noneq%kname = "$no_noneq="
    
    cntrl%reaktor%kname = "$reaktor="

    cntrl%k_on_demand%kname = "$k_on_demand="

    cntrl%dump_nzs%kname = "$dump_nzs="
    
    cntrl%diag_fix_ef%kname = "$diag_fix_ef="
    cntrl%diag_fixH%kname = "$diag_fixH="
    cntrl%diag_fixK%kname = "$diag_fixK="
    cntrl%diag_dim%kname = "$diag_dim="
    cntrl%diag_skip_check%kname = "$diag_skip_check="
    cntrl%diag_nkz%kname = "$diag_nkz="
    
    
    cntrl%dftu%kname = "$dftu="
    cntrl%dftu_file_c%kname = "$dftu_file_c="
    cntrl%dftu_file_l%kname = "$dftu_file_l="
    cntrl%dftu_file_r%kname = "$dftu_file_r="
    cntrl%dftu_projector%kname = "$dftu_projector="
    cntrl%dftu_sparse_projector%kname = "$dftu_sparse_projector="
    
    cntrl%spin%kname = "$spin="
    
    cntrl%hartree_shift%kname = "$hartree_shift="

  end subroutine init_cntrl

  subroutine get_key(keyname, keystr, keyint, keyreal, keylogic, l_iscritical)
    use petsc, only : PETSC_COMM_WORLD
    use kinds
    use misc
    use globals, only: inode, iunit_control
    implicit none

    character(strln) :: keyname
    character(strln), optional :: keystr
    integer, optional :: keyint
    real(dp), optional :: keyreal
    logical, optional :: keylogic
    logical :: l_iscritical

    integer :: io, ierr
    character(strln) :: instr
    logical ::found
  
    character(256) :: outstr

    rewind (iunit_control)

    io = 0
    found = .false.

    do while (io .eq. 0)

      read (iunit_control, fmt='(A256)', IOSTAT=io) instr
      if (io .ne. 0) exit
      instr = adjustl(instr)
      if (index(instr, trim(keyname)) .ge. 1) then
        found = .true.
        instr = instr(1:index(instr//"#", "#") - 1)
        instr = instr(index(instr, trim(keyname)) + len(trim(keyname)):strln)
        if (present(keystr)) call str2(instr, keystr)
        if (present(keyint)) call str2(instr, keyint)
        if (present(keyreal)) call str2(instr, keyreal)
        if (present(keylogic)) call str2(instr, keylogic)
        write (outstr, fmt='(A)') trim(keyname)//" "//trim(instr)        
        call PetscPrintf(PETSC_COMM_WORLD, trim(outstr)//New_line('A'), ierr)
        return
      end if

    end do

    if ((.not. found) .and. (l_iscritical)) then
      close (iunit_control)
      write (0, *) "non optional keyword ", trim(keyname), " not found in transp.ini"
      stop
    else
      if (present(keystr)) write(instr, fmt='(A)') keystr
      if (present(keyint)) write(instr, fmt='(i16)') keyint
      if (present(keyreal)) write(instr, fmt='(e24.12)') keyreal
      if (present(keylogic)) write(instr, fmt='(l)') keylogic
!~       write (outstr, fmt='(A)') trim(keyname)//" "//trim(instr)        
      write (outstr, fmt=*) trim(keyname)
      call PetscPrintf(PETSC_COMM_WORLD, trim(outstr)//New_line('A'), ierr)
      write (outstr, fmt=*) trim(instr)        
      call PetscPrintf(PETSC_COMM_WORLD, trim(outstr)//New_line('A'), ierr)
    end if

  end subroutine get_key

end module control

subroutine get_skey(keyname, kvalue, l_iscritical)
  use kinds
  use control, only: get_key
  implicit none
  character(strln) :: keyname
  character(strln) :: kvalue
  logical, optional :: l_iscritical

  logical :: l_critical

  l_critical = .true.
  if (present(l_iscritical)) l_critical = l_iscritical

  call get_key(keyname, keystr=kvalue, l_iscritical=l_critical)

end subroutine get_skey

subroutine get_rkey(keyname, kvalue, l_iscritical)
  use kinds
  use control, only: get_key
  implicit none
  character(strln) :: keyname
  real(dp) :: kvalue
  logical, optional :: l_iscritical

  logical :: l_critical

  l_critical = .true.
  if (present(l_iscritical)) l_critical = l_iscritical

  call get_key(keyname, keyreal=kvalue, l_iscritical=l_critical)

end subroutine get_rkey

subroutine get_ikey(keyname, kvalue, l_iscritical)
  use kinds
  use control, only: get_key
  implicit none
  character(strln) :: keyname
  integer :: kvalue
  logical, optional :: l_iscritical

  logical :: l_critical

  l_critical = .true.
  if (present(l_iscritical)) l_critical = l_iscritical

  call get_key(keyname, keyint=kvalue, l_iscritical=l_critical)

end subroutine get_ikey

subroutine get_lkey(keyname, kvalue, l_iscritical)
  use kinds
  use control, only: get_key
  implicit none
  character(strln) :: keyname
  logical :: kvalue
  logical, optional :: l_iscritical

  logical :: l_critical

  l_critical = .true.
  if (present(l_iscritical)) l_critical = l_iscritical

  call get_key(keyname, keylogic=kvalue, l_iscritical=l_critical)

end subroutine get_lkey

