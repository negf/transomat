module timing
  use kinds
  implicit none
! timing
  real(dp) :: timing_surf, timing_matsolve, timing_petsc_load, timing_petsc_save
  real(dp) :: timing_matmat_solve(2)

contains

  subroutine init_timings()
    implicit none

    timing_surf = 0d0
    timing_matsolve = 0d0
    timing_petsc_load = 0d0
    timing_petsc_save = 0d0
    timing_matmat_solve = 0d0

  end subroutine init_timings

end module timing
