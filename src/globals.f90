module globals
#include <petsc/finclude/petscmat.h>
  use petscmat
  use kinds
  implicit none

  complex(dp), parameter :: zione = cmplx(0d0, 1d0, 8)

  integer :: iigc

  real(dp) :: pi = atan(1d0)*4d0, eh = 27.211396132d0
  real(dp) :: kb = (8.6173303d-5)/27.211396132d0

  logical :: ksym = .true.

! mpi
  logical :: l_ionode

! file handeling
  logical :: l_is_big_endian

! petsc
  integer :: solver_mode, nsim_rhs
  integer :: inode, nprocs, psubcomm, inode_sub, inode_group, ngroups, group_range(2)
  integer, allocatable :: nrow_pproc(:, :), nrow_pproc_elec_l(:, :), nrow_pproc_elec_r(:, :), &
    nodes_group(:)
  logical :: l_loadsavegf

  MatType :: mattype_surf, mattype_cc, mattype_surf2, mattype_dense, mattype_sparse, &
    mattype_electrode
  MatSolverType :: matsolvertype_surf, matsolvertype_cc
  Integer, allocatable :: cols_loc(:), nzrow_loc(:)

  integer :: nsim_inv

! current k point
  integer :: iik
  logical :: l_k_on_demand

! integration control and parameters
  integer :: nint_order, maxsub, integrator, contour_select, ineq, int_counter, nint_order_el, select_lr, lr_select
  real(dp) :: epsfermi, delta_imag, eps_int, elow, eps_int_contour, mu_fermi_pol
  logical :: ldouble_contour, lnoneq_int, l_output_progress, l_reaktor, l_no_noneq

! transport parameters bias,temperature,chemical potentials
  real(dp) :: vb, mul, mur, temperature_el, ef_l, ef_c, ef_r, vl, vr, ef_c_fix

! transport parameters phonons
  real(dp) :: eta_ph_cc, eta_ph_elec, temperature_ph

! transport parameters for ep coupling
  integer :: ep_active_atoms(2), n_ep_active_atoms, nat3, n_ep_modes_max, ep_bias_n
  real(dp) :: ep_bias(2)
  integer, allocatable :: n_ep_modes_k(:)
  logical :: l_ep_loe, l_ep, l_ep_reuse_lambda

! transport control
  real(dp) :: eta_cc, eta_elec, conv_dez, eps_geo, eps_real, estart, eend, kk(3)
  character(strln) :: ecc_dir, elec_l_dir, elec_r_dir, trans_file

  integer, allocatable :: species_ecc(:), species_elec_l(:), species_elec_r(:)
  integer :: nkx_dez, nky_dez, nkz_diag, n_energy_steps, i_energy_start, i_energy_end
  integer :: nfiles_c, nfiles_l, nfiles_r
  integer :: iunit_control

  logical :: oldnegf, oneshot, calc_current_density, converged, init_scf, lget_ti

! lattice parameters
  real(dp) :: dlat_l(3, 3), rlat_l(3, 3), dlat_r(3, 3), rlat_r(3, 3), dlat_c(3, 3), &
 &  rlat_c(3, 3), dlat_l3(3,3), dlat_r3(3, 3)
  real(dp), allocatable :: xyz_ecc(:, :), xyz_elec_l(:, :), xyz_elec_r(:, :)
  real(dp), allocatable :: xyz_elec_l3(:, :), xyz_elec_r3(:, :)
  real(dp), allocatable :: kp_l(:, :), kp_r(:, :), kp_c(:,:)
  integer, allocatable :: wkp_r(:), wkp_l(:), wkp_c(:)
  logical :: kpoints_from_file
  integer :: nktot
  integer :: nred_l(2), nred_r(2)

! system info
  character(1), allocatable :: lcr_info(:)
  integer :: nat_ecc, nat_elec_l, nat_elec_r, nmu_ecc, nmu_elec_l, nmu_elec_r, nmu_c, nmu_l, nmu_r
  integer :: nmat_l3, nmat_r3, nat_elec_l3, nat_elec_r3
  integer :: nx_l_max, ny_l_max, nx_r_max, ny_r_max, nk_l, nk_r, nk_c, nx_max_c, ny_max_c
  integer :: ncell_c(3), ncell_l(3), ncell_r(3), nmat_c, nmat_l, nmat_r
  real(dp) :: n_electrons_ecc(2), n_electrons_l(2), n_electrons_r(2)
! mapping to conquest format for read into conquest
  integer, allocatable :: imu_ecc(:), imu_elec_l(:), imu_elec_r(:), atoms_c(:), imu_to_at(:)
  integer, allocatable :: tomat(:, :), inao_c(:), neigh_c(:), ndim_c(:), nlist_c(:)
  real(dp), allocatable :: tomatr(:, :), tomatr_l(:, :), tomatr_r(:, :)
  integer, allocatable :: tomat_l(:, :), inao_l(:), neigh_l(:), ndim_l(:), nlist_l(:), atoms_l(:)
  integer, allocatable :: tomat_r(:, :), inao_r(:), neigh_r(:), ndim_r(:), nlist_r(:), atoms_r(:)
! ----
  integer :: ndim_total, nzmax_row, spin_channel
  character(strln) :: spin_str
  

  complex(dp) :: zdosl(2), zdosr(2)

! DFT+Sigma
  real(dp) ::  d_occ, d_virt
  integer :: imu_dftsigma, jmu_dftsigma
  logical :: dftsigma

! scf
  real(dp) :: scf_conv

! debug output
  logical :: l_dump_nzs

! k mat mode
  integer :: k_mat_mode
  real(dp) :: hartree_shift

! diag
  logical :: l_diag_fix_ef, l_diag_fixH, l_diag_fixK, l_diag_skip_check
  integer :: diag_dim
  
! system
  logical l_use_sigma_l, l_use_sigma_r
  
  
! species

  type t_species
    integer :: nbf, nUshells_at
    integer, allocatable :: l(:), ibf(:), U_shell(:)
    logical, allocatable :: l_U(:)
    real(dp), allocatable :: U(:)
  end type t_species
  
  type(t_species), allocatable :: species_c(:), species_r(:), species_l(:)
  integer :: nspecies_c, nspecies_l, nspecies_r
  
! dft+U
  logical :: l_dftu, l_dftu_sparse_projector
  character(strln) :: dftu_file_c, dftu_file_l, dftu_file_r
  integer :: dftu_projector
  

end module globals
