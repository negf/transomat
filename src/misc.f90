module misc
  implicit none

  interface str2

    subroutine str2int4(instr, outint)
      implicit none
      character(*) :: instr
      integer(4) :: outint
    end subroutine str2int4

    subroutine str2int8(instr, outint)
      implicit none
      character(*) :: instr
      integer(8) :: outint
    end subroutine str2int8

    subroutine str2real(instr, outreal)
      use kinds
      implicit none
      character(*) :: instr
      real(dp) :: outreal
    end subroutine str2real

    subroutine str2str(instr, outstr)
      use kinds
      implicit none
      character(*) :: instr
      character(*) :: outstr
    end subroutine str2str

    subroutine str2logic(instr, outlogic)
      use kinds
      implicit none
      character(*) :: instr
      logical :: outlogic
    end subroutine str2logic

  end interface str2

  interface int2str
    function i2s(i)
      character(256) :: i2s
      integer :: i
    end function i2s
  end interface int2str

end module misc

subroutine str2real(str, outreal)
  use kinds
  implicit none

  character(*) :: str
  real(dp) :: outreal

  read (str, *) outreal

end subroutine str2real

subroutine str2int4(str, outint)

  implicit none

  character(*) :: str
  integer(4) :: outint

  read (str, *) outint

end subroutine str2int4

subroutine str2int8(str, outint)

  implicit none

  character(*) :: str
  integer(8) :: outint

  read (str, *) outint

end subroutine str2int8

subroutine str2str(str, outstr)

  implicit none

  character(*) :: str
  character(*) :: outstr

  outstr = str

end subroutine str2str

subroutine str2logic(str, outlogic)

  implicit none

  character(*) :: str
  logical :: outlogic
  integer :: istr

  read (str, *) istr

  outlogic = 1 .eq. istr

end subroutine str2logic

!~   subroutine int2str(i,str)
!~     character(*) :: str
!~     integer :: i

!~     write(str,fmt='(i8)') i
!~     str=adjustl(str)

!~   end subroutine int2str

function i2s(i)
  character(256) :: i2s
  integer :: i

  write (i2s, fmt='(i8)') i
  i2s = adjustl(i2s)

end function i2s
