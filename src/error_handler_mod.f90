module error_handler
  use kinds
  implicit none
  character(strln) :: errormsg

contains

  subroutine error()
#include <slepc/finclude/slepceps.h>
    use slepceps
    implicit none

    integer :: ierr

    write (0, fmt='(A)') trim(errormsg)

    call SlepcFinalize(ierr)
    stop 666

  end subroutine error

end module error_handler
