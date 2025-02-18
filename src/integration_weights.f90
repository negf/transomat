module integrations_weights
!~ from FORTRAN90 Source Codes
!~ https://people.sc.fsu.edu/~jburkardt/f_src/f_src.html
!~ released under GNU LGPL
  implicit none
  contains

    subroutine line_ncc_rule(n, a, b, x, w)

    !*****************************************************************************80
    !
    !! LINE_NCC_RULE computes a Newton-Cotes Closed (NCC) quadrature rule.
    !
    !  Discussion:
    !
    !    The integral:
    !
    !      Integral ( A <= X <= B ) F(X) dx
    !
    !    The quadrature rule:
    !
    !      Sum ( 1 <= I <= N ) W(I) * F ( X(I) ).
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    09 April 2014
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) N, the order.
    !
    !    Input, real ( kind = 8 ) A, B, the endpoints of the interval.
    !
    !    Input, real ( kind = 8 ) X(N), the abscissas.
    !
    !    Output, real ( kind = 8 ) W(N), the weights.
    !
      implicit none

      integer(kind=4) n

      real(kind=8) a
      real(kind=8) b
      real(kind=8) d(n)
      integer(kind=4) i
      integer(kind=4) j
      integer(kind=4) k
      real(kind=8) w(n)
      real(kind=8) x(n)
      real(kind=8) y_a
      real(kind=8) y_b
    !
    !  Define the points X.
    !
      call r8vec_linspace(n, a, b, x)
    !
    !  Compute the Lagrange basis polynomial which is 1 at X(I),
    !  and zero at the other nodes.
    !
      do i = 1, n

        d(1:n) = 0.0D+00
        d(i) = 1.0D+00

        do j = 2, n
          do k = j, n
            d(n + j - k) = (d(n + j - k - 1) - d(n + j - k))/(x(n + 1 - k) - x(n + j - k))
          end do
        end do

        do j = 1, n - 1
          do k = 1, n - j
            d(n - k) = d(n - k) - x(n - k - j + 1)*d(n - k + 1)
          end do
        end do
    !
    !  Evaluate the antiderivative of the polynomial at the endpoints.
    !
        y_a = d(n)/real(n, kind=8)
        do j = n - 1, 1, -1
          y_a = y_a*a + d(j)/real(j, kind=8)
        end do
        y_a = y_a*a

        y_b = d(n)/real(n, kind=8)
        do j = n - 1, 1, -1
          y_b = y_b*b + d(j)/real(j, kind=8)
        end do
        y_b = y_b*b

        w(i) = y_b - y_a

      end do

      return
    end subroutine line_ncc_rule

    subroutine r8vec_linspace(n, a, b, x)

    !*****************************************************************************80
    !
    !! R8VEC_LINSPACE creates a vector of linearly spaced values.
    !
    !  Discussion:
    !
    !    An R8VEC is a vector of R8's.
    !
    !    4 points evenly spaced between 0 and 12 will yield 0, 4, 8, 12.
    !
    !    In other words, the interval is divided into N-1 even subintervals,
    !    and the endpoints of intervals are used as the points.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    14 March 2011
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) N, the number of entries in the vector.
    !
    !    Input, real ( kind = 8 ) A, B, the first and last entries.
    !
    !    Output, real ( kind = 8 ) X(N), a vector of linearly spaced data.
    !
      implicit none

      integer(kind=4) n

      real(kind=8) a
      real(kind=8) b
      integer(kind=4) i
      real(kind=8) x(n)

      if (n == 1) then

        x(1) = (a + b)/2.0D+00

      else

        do i = 1, n
          x(i) = (real(n - i, kind=8)*a &
                  + real(i - 1, kind=8)*b) &
                 /real(n - 1, kind=8)
        end do

      end if

      return
    end subroutine r8vec_linspace
  
    subroutine clenshaw_curtis_compute(order, x, w)

    !*****************************************************************************80
    !
    !! CLENSHAW_CURTIS_COMPUTE computes a Clenshaw Curtis quadrature rule.
    !
    !  Discussion:
    !
    !    Our convention is that the abscissas are numbered from left to right.
    !
    !    The rule is defined on [-1,1].
    !
    !    The integral to approximate:
    !
    !      Integral ( -1 <= X <= 1 ) F(X) dX
    !
    !    The quadrature rule:
    !
    !      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    15 February 2009
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) ORDER, the order of the rule.
    !    1 <= ORDER.
    !
    !    Output, real ( kind = 8 ) X(ORDER), the abscissas.
    !
    !    Output, real ( kind = 8 ) W(ORDER), the weights.
    !
      implicit none

      integer(kind=4) order

      real(kind=8) b
      integer(kind=4) i
      integer(kind=4) j
      real(kind=8), parameter :: r8_pi = 3.141592653589793D+00
      real(kind=8) theta
      real(kind=8) w(order)
      real(kind=8) x(order)

      if (order < 1) then
        write (*, '(a)') ' '
        write (*, '(a)') 'CLENSHAW_CURTIS_COMPUTE - Fatal error!'
        write (*, '(a,i8)') '  Illegal value of ORDER = ', order
        stop
      end if

      if (order == 1) then
        x(1) = 0.0D+00
        w(1) = 2.0D+00
        return
      end if

      do i = 1, order
        x(i) = cos(real(order - i, kind=8)*r8_pi &
                   /real(order - 1, kind=8))
      end do

      x(1) = -1.0D+00
      if (mod(order, 2) == 1) then
        x((order + 1)/2) = 0.0D+00
      end if
      x(order) = +1.0D+00

      do i = 1, order

        theta = real(i - 1, kind=8)*r8_pi &
                /real(order - 1, kind=8)

        w(i) = 1.0D+00

        do j = 1, (order - 1)/2

          if (2*j == (order - 1)) then
            b = 1.0D+00
          else
            b = 2.0D+00
          end if

          w(i) = w(i) - b*cos(2.0D+00*real(j, kind=8)*theta) &
                 /real(4*j*j - 1, kind=8)

        end do

      end do

      w(1) = w(1)/real(order - 1, kind=8)
      w(2:order - 1) = 2.0D+00*w(2:order - 1)/real(order - 1, kind=8)
      w(order) = w(order)/real(order - 1, kind=8)

      return
    end subroutine clenshaw_curtis_compute  
  
    subroutine kronrod(n, eps, x, w1, w2)

    !*****************************************************************************80
    !
    !! KRONROD adds N+1 points to an N-point Gaussian rule.
    !
    !  Discussion:
    !
    !    This subroutine calculates the abscissas and weights of the 2N+1
    !    point Gauss Kronrod quadrature formula which is obtained from the
    !    N point Gauss quadrature formula by the optimal addition of N+1 points.
    !
    !    The optimally added points are called Kronrod abscissas.  The
    !    abscissas and weights for both the Gauss and Gauss Kronrod rules
    !    are calculated for integration over the interval [-1,+1].
    !
    !    Since the quadrature formula is symmetric with respect to the origin,
    !    only the nonnegative abscissas are calculated.
    !
    !    Note that the code published in Mathematics of Computation
    !    omitted the definition of the variable which is here called COEF2.
    !
    !  Storage:
    !
    !    Given N, let M = ( N + 1 ) / 2.
    !
    !    The Gauss-Kronrod rule will include 2*N+1 points.  However, by symmetry,
    !    only N + 1 of them need to be listed.
    !
    !    The arrays X, W1 and W2 contain the nonnegative abscissas in decreasing
    !    order, and the weights of each abscissa in the Gauss-Kronrod and
    !    Gauss rules respectively.  This means that about half the entries
    !    in W2 are zero.
    !
    !    For instance, if N = 3, the output is:
    !
    !    I      X               W1              W2
    !
    !    1    0.960491        0.104656         0.000000
    !    2    0.774597        0.268488         0.555556
    !    3    0.434244        0.401397         0.000000
    !    4    0.000000        0.450917         0.888889
    !
    !    and if N = 4, (notice that 0 is now a Kronrod abscissa)
    !    the output is
    !
    !    I      X               W1              W2
    !
    !    1    0.976560        0.062977        0.000000
    !    2    0.861136        0.170054        0.347855
    !    3    0.640286        0.266798        0.000000
    !    4    0.339981        0.326949        0.652145
    !    5    0.000000        0.346443        0.000000
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    03 August 2010
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Robert Piessens, Maria Branders.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Robert Piessens, Maria Branders,
    !    A Note on the Optimal Addition of Abscissas to Quadrature Formulas
    !    of Gauss and Lobatto,
    !    Mathematics of Computation,
    !    Volume 28, Number 125, January 1974, pages 135-139.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) N, the order of the Gauss rule.
    !
    !    Input, real ( kind = 8 ) EPS, the requested absolute accuracy of the
    !    abscissas.
    !
    !    Output, real ( kind = 8 ) X(N+1), the abscissas.
    !
    !    Output, real ( kind = 8 ) W1(N+1), the weights for
    !    the Gauss-Kronrod rule.
    !
    !    Output, real ( kind = 8 ) W2(N+1), the weights for
    !    the Gauss rule.
    !
    implicit none

    integer(kind=4) n

    real(kind=8) ak
    real(kind=8) an
    real(kind=8) b(((n + 1)/2) + 1)
    real(kind=8) bb
    real(kind=8) c
    real(kind=8) coef
    real(kind=8) coef2
    real(kind=8) d
    real(kind=8) eps
    logical even
    integer(kind=4) i
    integer(kind=4) k
    integer(kind=4) l
    integer(kind=4) ll
    integer(kind=4) m
    real(kind=8) s
    real(kind=8) tau((n + 1)/2)
    real(kind=8) w1(n + 1)
    real(kind=8) w2(n + 1)
    real(kind=8) x(n + 1)
    real(kind=8) x1
    real(kind=8) xx
    real(kind=8) y

    m = (n + 1)/2
    even = (2*m == n)

    d = 2.0D+00
    an = 0.0D+00
    do k = 1, n
    an = an + 1.0D+00
    d = d*an/(an + 0.5D+00)
    end do
    !
    !  Calculation of the Chebyshev coefficients of the orthogonal polynomial.
    !
    tau(1) = (an + 2.0D+00)/(an + an + 3.0D+00)
    b(m) = tau(1) - 1.0D+00
    ak = an

    do l = 1, m - 1

    ak = ak + 2.0D+00
    tau(l + 1) = ((ak - 1.0D+00)*ak &
                  - an*(an + 1.0D+00))*(ak + 2.0D+00)*tau(l) &
                 /(ak*((ak + 3.0D+00)*(ak + 2.0D+00) &
                       - an*(an + 1.0D+00)))
    b(m - l) = tau(l + 1)

    do ll = 1, l
      b(m - l) = b(m - l) + tau(ll)*b(m - l + ll)
    end do

    end do

    b(m + 1) = 1.0D+00
    !
    !  Calculation of approximate values for the abscissas.
    !
    bb = sin(1.570796D+00/(an + an + 1.0D+00))
    x1 = sqrt(1.0D+00 - bb*bb)
    s = 2.0D+00*bb*x1
    c = sqrt(1.0D+00 - s*s)
    coef = 1.0D+00 - (1.0D+00 - 1.0D+00/an)/(8.0D+00*an*an)
    xx = coef*x1
    !
    !  Coefficient needed for weights.
    !
    !  COEF2 = 2^(2*n+1) * n! * n! / (2n+1)!
    !
    coef2 = 2.0D+00/real(2*n + 1, kind=8)
    do i = 1, n
    coef2 = coef2*4.0D+00*real(i, kind=8)/real(n + i, kind=8)
    end do
    !
    !  Calculation of the K-th abscissa (a Kronrod abscissa) and the
    !  corresponding weight.
    !
    do k = 1, n, 2

    call abwe1(n, m, eps, coef2, even, b, xx, w1(k))
    w2(k) = 0.0D+00

    x(k) = xx
    y = x1
    x1 = y*c - bb*s
    bb = y*s + bb*c

    if (k == n) then
      xx = 0.0D+00
    else
      xx = coef*x1
    end if
    !
    !  Calculation of the K+1 abscissa (a Gaussian abscissa) and the
    !  corresponding weights.
    !
    call abwe2(n, m, eps, coef2, even, b, xx, w1(k + 1), w2(k + 1))

    x(k + 1) = xx
    y = x1
    x1 = y*c - bb*s
    bb = y*s + bb*c
    xx = coef*x1

    end do
    !
    !  If N is even, we have one more Kronrod abscissa to compute,
    !  namely the origin.
    !
    if (even) then
    xx = 0.0D+00
    call abwe1(n, m, eps, coef2, even, b, xx, w1(n + 1))
    w2(n + 1) = 0.0D+00
    x(n + 1) = xx
    end if

    return
    end subroutine kronrod
    subroutine kronrod_adjust(a, b, n, x, w1, w2)

    !*****************************************************************************80
    !
    !! KRONROD_ADJUST adjusts a Gauss-Kronrod rule from [-1,+1] to [A,B].
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    23 August 2015
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) A, B, the endpoints of the new interval.
    !
    !    Input, integer ( kind = 4 ) N, the order of the rule.
    !
    !    Input/output, real ( kind = 8 ) X(N+1), W1(N+1), W2(N+1), the abscissas
    !    and weights.
    !
    implicit none

    integer(kind=4) n

    real(kind=8) a
    real(kind=8) b
    real(kind=8) w1(n + 1)
    real(kind=8) w2(n + 1)
    real(kind=8) x(n + 1)

    x(1:n + 1) = ((1.0D+00 - x(1:n + 1))*a &
                + (1.0D+00 + x(1:n + 1))*b) &
               /2.0D+00

    w1(1:n + 1) = ((b - a)/2.0D+00)*w1(1:n + 1)
    w2(1:n + 1) = ((b - a)/2.0D+00)*w2(1:n + 1)

    return
    end subroutine kronrod_adjust
    subroutine abwe1(n, m, eps, coef2, even, b, x, w)

    !*****************************************************************************80
    !
    !! ABWE1 calculates a Kronrod abscissa and weight.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    03 August 2010
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Robert Piessens, Maria Branders.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Robert Piessens, Maria Branders,
    !    A Note on the Optimal Addition of Abscissas to Quadrature Formulas
    !    of Gauss and Lobatto,
    !    Mathematics of Computation,
    !    Volume 28, Number 125, January 1974, pages 135-139.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) N, the order of the Gauss rule.
    !
    !    Input, integer ( kind = 4 ) M, the value of ( N + 1 ) / 2.
    !
    !    Input, real ( kind = 8 ) EPS, the requested absolute accuracy of the
    !    abscissas.
    !
    !    Input, real ( kind = 8 ) COEF2, a value needed to compute weights.
    !
    !    Input, logical EVEN, is TRUE if N is even.
    !
    !    Input, real ( kind = 8 ) B(M+1), the Chebyshev coefficients.
    !
    !    Input/output, real ( kind = 8 ) X; on input, an estimate for
    !    the abscissa, and on output, the computed abscissa.
    !
    !    Output, real ( kind = 8 ) W, the weight.
    !
    implicit none

    integer(kind=4) m

    real(kind=8) ai
    real(kind=8) b(m + 1)
    real(kind=8) b0
    real(kind=8) b1
    real(kind=8) b2
    real(kind=8) coef2
    real(kind=8) d0
    real(kind=8) d1
    real(kind=8) d2
    real(kind=8) delta
    real(kind=8) dif
    real(kind=8) eps
    logical even
    real(kind=8) f
    real(kind=8) fd
    integer(kind=4) i
    integer(kind=4) iter
    integer(kind=4) k
    integer(kind=4) ka
    integer(kind=4) n
    real(kind=8) w
    real(kind=8) x
    real(kind=8) yy

    if (x == 0.0D+00) then
    ka = 1
    else
    ka = 0
    end if
    !
    !  Iterative process for the computation of a Kronrod abscissa.
    !
    do iter = 1, 50

    b1 = 0.0D+00
    b2 = b(m + 1)
    yy = 4.0D+00*x*x - 2.0D+00
    d1 = 0.0D+00

    if (even) then
      ai = m + m + 1
      d2 = ai*b(m + 1)
      dif = 2.0D+00
    else
      ai = m + 1
      d2 = 0.0D+00
      dif = 1.0D+00
    end if

    do k = 1, m
      ai = ai - dif
      i = m - k + 1
      b0 = b1
      b1 = b2
      d0 = d1
      d1 = d2
      b2 = yy*b1 - b0 + b(i)
      if (.not. even) then
        i = i + 1
      end if
      d2 = yy*d1 - d0 + ai*b(i)
    end do

    if (even) then
      f = x*(b2 - b1)
      fd = d2 + d1
    else
      f = 0.5D+00*(b2 - b0)
      fd = 4.0D+00*x*d2
    end if
    !
    !  Newton correction.
    !
    delta = f/fd
    x = x - delta

    if (ka == 1) then
      exit
    end if

    if (abs(delta) <= eps) then
      ka = 1
    end if

    end do
    !
    !  Catch non-convergence.
    !
    if (ka /= 1) then
    write (*, '(a)') ' '
    write (*, '(a)') 'ABWE1 - Fatal error!'
    write (*, '(a)') '  Iteration limit reached.'
    write (*, '(a,g14.6)') '  EPS is ', eps
    write (*, '(a,g14.6)') '  Last DELTA was ', delta
    stop 1
    end if
    !
    !  Computation of the weight.
    !
    d0 = 1.0D+00
    d1 = x
    ai = 0.0D+00
    do k = 2, n
    ai = ai + 1.0D+00
    d2 = ((ai + ai + 1.0D+00)*x*d1 - ai*d0)/(ai + 1.0D+00)
    d0 = d1
    d1 = d2
    end do

    w = coef2/(fd*d2)

    return
    end subroutine abwe1
    subroutine abwe2(n, m, eps, coef2, even, b, x, w1, w2)

    !*****************************************************************************80
    !
    !! ABWE2 calculates a Gaussian abscissa and two weights.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    30 April 2013
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Robert Piessens, Maria Branders.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Robert Piessens, Maria Branders,
    !    A Note on the Optimal Addition of Abscissas to Quadrature Formulas
    !    of Gauss and Lobatto,
    !    Mathematics of Computation,
    !    Volume 28, Number 125, January 1974, pages 135-139.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) N, the order of the Gauss rule.
    !
    !    Input, integer ( kind = 4 ) M, the value of ( N + 1 ) / 2.
    !
    !    Input, real ( kind = 8 ) EPS, the requested absolute accuracy of the
    !    abscissas.
    !
    !    Input, real ( kind = 8 ) COEF2, a value needed to compute weights.
    !
    !    Input, logical EVEN, is TRUE if N is even.
    !
    !    Input, real ( kind = 8 ) B(M+1), the Chebyshev coefficients.
    !
    !    Input/output, real ( kind = 8 ) X; on input, an estimate for
    !    the abscissa, and on output, the computed abscissa.
    !
    !    Output, real ( kind = 8 ) W1, the Gauss-Kronrod weight.
    !
    !    Output, real ( kind = 8 ) W2, the Gauss weight.
    !
    implicit none

    integer(kind=4) m

    real(kind=8) ai
    real(kind=8) an
    real(kind=8) b(m + 1)
    real(kind=8) coef2
    real(kind=8) delta
    real(kind=8) eps
    logical even
    integer(kind=4) i
    integer(kind=4) iter
    integer(kind=4) k
    integer(kind=4) ka
    integer(kind=4) n
    real(kind=8) p0
    real(kind=8) p1
    real(kind=8) p2
    real(kind=8) pd0
    real(kind=8) pd1
    real(kind=8) pd2
    real(kind=8) w1
    real(kind=8) w2
    real(kind=8) x
    real(kind=8) yy

    if (x == 0.0D+00) then
    ka = 1
    else
    ka = 0
    end if
    !
    !  Iterative process for the computation of a Gaussian abscissa.
    !
    do iter = 1, 50

    p0 = 1.0D+00
    p1 = x
    pd0 = 0.0D+00
    pd1 = 1.0D+00
    !
    !  When N is 1, we need to initialize P2 and PD2 to avoid problems with DELTA.
    !
    if (n <= 1) then
      if (epsilon(x) < abs(x)) then
        p2 = (3.0D+00*x*x - 1.0D+00)/2.0D+00
        pd2 = 3.0D+00*x
      else
        p2 = 3.0D+00*x
        pd2 = 3.0D+00
      end if
    end if

    ai = 0.0D+00
    do k = 2, n
      ai = ai + 1.0D+00
      p2 = ((ai + ai + 1.0D+00)*x*p1 - ai*p0)/(ai + 1.0D+00)
      pd2 = ((ai + ai + 1.0D+00)*(p1 + x*pd1) - ai*pd0) &
            /(ai + 1.0D+00)
      p0 = p1
      p1 = p2
      pd0 = pd1
      pd1 = pd2
    end do
    !
    !  Newton correction.
    !
    delta = p2/pd2
    x = x - delta

    if (ka == 1) then
      exit
    end if

    if (abs(delta) <= eps) then
      ka = 1
    end if

    end do
    !
    !  Catch non-convergence.
    !
    if (ka /= 1) then
    write (*, '(a)') ' '
    write (*, '(a)') 'ABWE2 - Fatal error!'
    write (*, '(a)') '  Iteration limit reached.'
    write (*, '(a,g14.6)') '  EPS is ', eps
    write (*, '(a,g14.6)') '  Last DELTA was ', delta
    stop 1
    end if
    !
    !  Computation of the weight.
    !
    an = n

    w2 = 2.0D+00/(an*pd2*p0)

    p1 = 0.0D+00
    p2 = b(m + 1)
    yy = 4.0D+00*x*x - 2.0D+00
    do k = 1, m
    i = m - k + 1
    p0 = p1
    p1 = p2
    p2 = yy*p1 - p0 + b(i)
    end do

    if (even) then
    w1 = w2 + coef2/(pd2*x*(p2 - p1))
    else
    w1 = w2 + 2.0D+00*coef2/(pd2*(p2 - p0))
    end if

    return
    end subroutine abwe2
    
end module integrations_weights
