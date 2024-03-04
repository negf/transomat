module mathstuff

  implicit none

  interface invert
    subroutine cinvert(A)
      implicit none
      double complex :: A(:, :)
    end subroutine cinvert

    subroutine rinvert(A)
      implicit none
      double precision :: A(:, :)
    end subroutine rinvert

  end interface invert

  interface matnorm_my
    function rmatnorm_my(a)
      implicit none
      double precision :: a(:, :)
      double precision :: rmatnorm_my
    end function rmatnorm_my

    function cmatnorm_my(a)
      implicit none
      double complex :: a(:, :)
      double precision :: cmatnorm_my
    end function cmatnorm_my
  end interface matnorm_my

contains

  function ininterval(i1, i2, j1, j2)
    integer :: i1, i2, j1, j2
    logical ::ininterval

    ininterval = .false.

    if ((j1 .ge. i1) .and. (j1 .le. i2)) ininterval = .true.

    if ((j2 .ge. i1) .and. (j1 .le. i1)) ininterval = .true.

  end function

  function matnorm_myab(a, b)
    implicit none

    double complex :: a(:, :), b(:, :)
    double precision :: matnorm_myab
    double precision :: trace
    double complex :: d
    integer :: n, i, j

    n = size(a, 1)

    trace = 0d0
    do j = 1, n
      do i = 1, n
        d = a(i, j) - b(i, j)
        trace = trace + d*conjg(d)
      end do
    end do
    matnorm_myab = sqrt(trace)

  end function

  function rmatdotab(a, b)
    implicit none

    double precision :: a(:, :), b(:, :)
    double precision :: rmatdotab
    double precision :: trace
    integer :: n, i, j

    n = size(a, 1)

    trace = 0d0
    do j = 1, n
      do i = 1, n
        trace = trace + a(i, j)*b(i, j)
      end do
    end do
    rmatdotab = trace

  end function

!~     subroutine squarerootinv(A)
!~       use lapack95
!~       use blas95
!~       use kinds
!~       implicit none

!~       real(dp), allocatable :: A(:,:)
!~       real(dp), allocatable :: w(:),b(:,:),tmp(:,:)
!~       integer :: n,ierr,i1
!~       real(dp) :: eps

!~       eps=1d-5
!~       n=size(A,1)
!~       allocate(w(n),b(n,n),tmp(n,n),stat=ierr)
!~       if (ierr.ne.0) then
!~         write(0,*) "allocation error ",ierr
!~       end if

!~       b=0d0
!~       call syev(a,w,'V','L',ierr)
!~       do i1=1,n
!~         if (w(i1).lt.0) then
!~           write(0,*) "get S^1/2 negative EV ",i1,w(i1)
!~           stop
!~         end if

!~         if (w(i1).gt.0d0) then
!~           b(i1,i1)=1d0/sqrt(w(i1))
!~         else
!~           b(i1,i1)=1d-16
!~         end if
!~       end do

!~       call gemm(a,b,tmp,'n','n')
!~       call gemm(tmp,a,b,'n','t')
!~       a=b

!~     end subroutine squarerootinv

!~         subroutine squareroot(A)
!~       use lapack95
!~       use blas95
!~       use kinds
!~       implicit none

!~       real(dp), allocatable :: A(:,:)
!~       real(dp), allocatable :: w(:),b(:,:),tmp(:,:)
!~       integer :: n,ierr,i1
!~       real(dp) :: eps

!~       eps=1d-5
!~       n=size(A,1)
!~       allocate(w(n),b(n,n),tmp(n,n),stat=ierr)
!~       if (ierr.ne.0) then
!~         write(0,*) "allocation error ",ierr
!~       end if

!~       b=0d0
!~       call syev(a,w,'V','L',ierr)
!~       do i1=1,n
!~         if (w(i1).lt.0) then
!~           write(0,*) "get S^1/2 negative EV ",i1,w(i1)
!~           w(i1)=0d0

!~         end if

!~         b(i1,i1)=sqrt(w(i1))
!~       end do

!~       call gemm(a,b,tmp,'n','n')
!~       call gemm(tmp,a,b,'n','t')
!~       a=b

!~     end subroutine squareroot

  function vector_product(a, b)
    use kinds
    implicit none
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: vector_product(3)

    vector_product(1) = a(2)*b(3) - a(3)*b(2)
    vector_product(2) = a(3)*b(1) - a(1)*b(3)
    vector_product(3) = a(1)*b(2) - a(2)*b(1)

  end function vector_product

end module

!~     subroutine cinvert(A)
!~       use lapack95
!~       implicit none
!~       double complex :: A(:,:)

!~       integer,allocatable :: ipiv(:)
!~       integer :: n,ierr

!~       n=size(A,1)

!~       allocate(ipiv(n),stat=ierr)
!~       if (ierr.ne.0) then
!~         write(0,*) "allocation error ipiv",ierr
!~         stop
!~       end if

!~       call getrf(A,ipiv,ierr)

!~       if (ierr.ne.0) then
!~         write(0,*) "error in getrf ",ierr
!~         stop
!~       end if

!~       call getri(A,ipiv,ierr)

!~       if (ierr.ne.0) then
!~         write(0,*) "error in getri ",ierr
!~         stop
!~       end if

!~     end subroutine cinvert

!~     subroutine rinvert(A)
!~       use lapack95
!~       implicit none
!~       double precision :: A(:,:)

!~       integer,allocatable :: ipiv(:)
!~       integer :: n,ierr

!~       n=size(A,1)

!~       allocate(ipiv(n),stat=ierr)
!~       if (ierr.ne.0) then
!~         write(0,*) "allocation error ipiv",ierr
!~         stop
!~       end if

!~       call getrf(A,ipiv,ierr)

!~       if (ierr.ne.0) then
!~         write(0,*) "error in getrf ",ierr
!~         stop
!~       end if
!~       call getri(A,ipiv,ierr)

!~       if (ierr.ne.0) then
!~         write(0,*) "error in getri ",ierr
!~         stop
!~       end if

!~     end subroutine rinvert

!~     function cmatnorm_my(a)
!~       implicit none

!~       double complex :: a(:,:)
!~       double precision :: trace,cmatnorm_my
!~       integer :: n,i,j

!~       n=size(a,1)

!~       trace=0d0
!~       do j=1,n
!~         do i=1,n
!~         trace=trace+a(i,j)*conjg(a(i,j))
!~         end do
!~       end do
!~       cmatnorm_my=sqrt(trace)

!~     end function cmatnorm_my

!~     function rmatnorm_my(a)
!~       implicit none

!~       double precision :: a(:,:)
!~       double precision :: trace,rmatnorm_my
!~       integer :: n,i,j

!~       n=size(a,1)

!~       trace=0d0
!~       do j=1,n
!~         do i=1,n
!~         trace=trace+a(i,j)*a(i,j)
!~         end do
!~       end do
!~       rmatnorm_my=sqrt(trace)

!~     end function rmatnorm_my
!~     subroutine cinvert(A)
!~       use lapack95
!~       implicit none
!~       double complex :: A(:,:)

!~       integer,allocatable :: ipiv(:)
!~       integer :: n,ierr

!~       n=size(A,1)

!~       allocate(ipiv(n),stat=ierr)
!~       if (ierr.ne.0) then
!~         write(0,*) "allocation error ipiv",ierr
!~         stop
!~       end if

!~       call getrf(A,ipiv,ierr)

!~       if (ierr.ne.0) then
!~         write(0,*) "error in getrf ",ierr
!~         stop
!~       end if

!~       call getri(A,ipiv,ierr)

!~       if (ierr.ne.0) then
!~         write(0,*) "error in getri ",ierr
!~         stop
!~       end if

!~     end subroutine cinvert

!~     subroutine rinvert(A)
!~       use lapack95
!~       implicit none
!~       double precision :: A(:,:)

!~       integer,allocatable :: ipiv(:)
!~       integer :: n,ierr

!~       n=size(A,1)

!~       allocate(ipiv(n),stat=ierr)
!~       if (ierr.ne.0) then
!~         write(0,*) "allocation error ipiv",ierr
!~         stop
!~       end if

!~       call getrf(A,ipiv,ierr)

!~       if (ierr.ne.0) then
!~         write(0,*) "error in getrf ",ierr
!~         stop
!~       end if
!~       call getri(A,ipiv,ierr)

!~       if (ierr.ne.0) then
!~         write(0,*) "error in getri ",ierr
!~         stop
!~       end if

!~     end subroutine rinvert

!~     function cmatnorm_my(a)
!~       implicit none

!~       double complex :: a(:,:)
!~       double precision :: trace,cmatnorm_my
!~       integer :: n,i,j

!~       n=size(a,1)

!~       trace=0d0
!~       do j=1,n
!~         do i=1,n
!~         trace=trace+a(i,j)*conjg(a(i,j))
!~         end do
!~       end do
!~       cmatnorm_my=sqrt(trace)

!~     end function cmatnorm_my

!~     function rmatnorm_my(a)
!~       implicit none

!~       double precision :: a(:,:)
!~       double precision :: trace,rmatnorm_my
!~       integer :: n,i,j

!~       n=size(a,1)

!~       trace=0d0
!~       do j=1,n
!~         do i=1,n
!~         trace=trace+a(i,j)*a(i,j)
!~         end do
!~       end do
!~       rmatnorm_my=sqrt(trace)

!~     end function rmatnorm_my
