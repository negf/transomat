
! check if the left/right electrode can be translated onto the l and r ends of the ecc.
! at the time being no reordering, rotations and different dimension between l/r_ecc and electrodes

subroutine check_sys(xyz_ecc, xyz_elec_l, xyz_elec_r, nat_ecc, nat_elec_l, nat_elec_r, lcr_info)
  use kinds
  use petsc_mod, only: pstr_out, petsc_print_master
  use globals, only: eps_geo, inode, l_use_sigma_l, l_use_sigma_r
  use error_handler
  implicit none

  integer :: nat_ecc, nat_elec_l, nat_elec_r
  real(dp) :: xyz_ecc(3, nat_ecc), xyz_elec_l(3, nat_elec_l), xyz_elec_r(3, nat_elec_r)
  character(1) :: lcr_info(nat_ecc)

  real(dp) :: transvec(3)
  integer :: i
  lcr_info = "c"
  
  if (l_use_sigma_l) then
    transvec = xyz_ecc(1:3, 1) - xyz_elec_l(1:3, 1)
!~     write(6,fmt='(A,3e24.12)') "first atom ecc",xyz_ecc(1:3,1)
!~     write(6,fmt='(A,3e24.12)') "first atom el ",xyz_elec_l(1:3,1)
!~     write(6,fmt='(A,3e24.12)') "transvec_r   ",transvec
    do i = 1, nat_elec_l
      if (.not. all(abs(xyz_ecc(1:3, i) - xyz_elec_l(1:3, i) - transvec) .le. eps_geo)) then
        write (pstr_out, fmt='(A)') "missmatch between ecc_l and left electrode"; call petsc_print_master()
        write (pstr_out, fmt='(A,3e24.12)') "transvec ", transvec; call petsc_print_master()
        write (pstr_out, fmt='(A,i8)') "ecc_l,elec_l,abs(diff)", i; call petsc_print_master()
        write (pstr_out, fmt='(3e24.12)') xyz_ecc(1:3,  i); call petsc_print_master()
        write (pstr_out, fmt='(3e24.12)') xyz_elec_l(1:3, i); call petsc_print_master()
        write (pstr_out, fmt='(3e24.12)') transvec; call petsc_print_master()
        write (pstr_out, fmt='(3e24.12)') xyz_ecc(1:3,  i) - xyz_elec_l(1:3, i) - transvec; call petsc_print_master()
        errormsg = "lattice error"
        call error()
      end if
      lcr_info(i) = "l"
    end do
  end if
!~   write(6,fmt='(A)') "left electrode ok"
  if (l_use_sigma_r) then
    transvec = xyz_ecc(1:3, nat_ecc) - xyz_elec_r(1:3, nat_elec_r)
!~   write(6,fmt='(A,3e24.12)') "last atom ecc",xyz_ecc(1:3,nat_ecc)
!~   write(6,fmt='(A,3e24.12)') "last atom er ",xyz_elec_r(1:3,nat_elec_r)
!~   write(6,fmt='(A,3e24.12)') "transvec_r   ",transvec
    do i = 1, nat_elec_r
      if (.not. all(abs(xyz_ecc(1:3, nat_ecc - nat_elec_r + i) - xyz_elec_r(1:3, i) - transvec) .le. eps_geo)) then
        write (pstr_out, fmt='(A)') "missmatch between ecc_r and right electrode"; call petsc_print_master()
        write (pstr_out, fmt='(A,3e24.12)') "transvec ", transvec; call petsc_print_master()
        write (pstr_out, fmt='(A,i8)') "ecc_r,elec_r,abs(diff)", i ; call petsc_print_master()
        write (pstr_out, fmt='(3e24.12)') xyz_ecc(1:3, nat_ecc - nat_elec_r + i); call petsc_print_master()
        write (pstr_out, fmt='(3e24.12)') xyz_elec_r(1:3, i); call petsc_print_master()
        write (pstr_out, fmt='(3e24.12)') transvec; call petsc_print_master()
        write (pstr_out, fmt='(3e24.12)') xyz_ecc(1:3, nat_ecc - nat_elec_r + i) - xyz_elec_r(1:3, i) - transvec; call petsc_print_master()
        errormsg = "lattice error"
        call error()
      end if
      lcr_info(nat_ecc - nat_elec_r + i) = "r"
    end do
!~   write(6,fmt='(A)') "right electrode ok"
  end if

end subroutine check_sys
