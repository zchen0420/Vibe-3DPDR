!=======================================================================
!
!     Copied from Numerical Recipes
!
! Given a tabulated function YA (of size MxN) and tabulated independent
! variables X1A (M values) and X2A (N values), this routine constructs
! one-dimensional natural cubic splines of the rows of YA and returns
! the second derivatives in the array Y2A.
!
!-----------------------------------------------------------------------
subroutine splie2(x1a,x2a,ya,m,n,y2a)

  use definitions
  use healpix_types

  implicit none

  integer(kind=i4b), intent(in) :: m,n
  real(kind=dp), intent(in)     :: x1a(1:m),x2a(1:n),ya(1:m,1:n)
  real(kind=dp), intent(out)    :: y2a(1:m,1:n)

  integer(kind=i4b) :: i,j
  real(kind=dp) :: ytmp(1:n),y2tmp(1:n)

  do i=1,m
    do j=1,n
      ytmp(j)=ya(i,j)
    end do
    !        Values of 1.0D30 indicate a natural spline
    call spline(x2a,ytmp,n,1.0d30,1.0d30,y2tmp)
    do j=1,n
      y2a(i,j)=y2tmp(j)
    end do
  end do

  return
end
!=======================================================================

!=======================================================================
!
!     Calculate the cubic spline for a set of points (X,Y)
!     (c.f. Numerical Recipes, Chapter 3.3: Spline Routine)
!
!     Given the arrays X and Y (size N) containing a tabulated
!     function, i.e., Y(I)=f(X(I)), with X(1) < X(2) < ... < X(N),
!     and given values YP1 and YPN for the first derivative of the
!     interpolating function at points 1 and N, respectively, this
!     routine returns an array Y2 of length N, which contains the
!     second derivatives of the interpolating function at the
!     tabulated points X(I). If YP1 and/or YPN are equal to 1.0E+30
!     or larger, the routine is signalled to set the corresponding
!     boundary condition for a natural spline, with zero second
!     derivative at that boundary.
!
!     I/O parameters:
!     Input   X   = vector for independent variable; dimension X(1:N)
!     Input   Y   = vector for x-dependent variable; dimension Y(1:N)
!     Input   N   = dimension of vectors containing tabulated function
!     Input   YP1 = 1st derivative of the function at point 1
!     Input   YPN = 1st derivative of the function at point N
!     Output  Y2  = 2nd derivative of the function; dimension Y2(1:N)
!
!-----------------------------------------------------------------------
subroutine spline(x,y,n,yp1,ypn,y2)

  use definitions
  use healpix_types

  implicit none

  integer(kind=i4b), intent(in) :: n
  real(kind=dp), intent(in)     :: x(1:n),y(1:n)
  real(kind=dp), intent(in)     :: yp1,ypn
  real(kind=dp), intent(out)    :: y2(1:n)

  integer(kind=i4b) :: i
  real(kind=dp) :: p,qn,sig,u(1:n),un

  if(yp1.ge.1.0d30) then
    !        The lower boundary condition is either set to be "natural"
    y2(1)=0.0d0
    u(1)=0.0d0
  else
    !        or to have a specified first derivative
    y2(1)=-0.5d0
    u(1)=(3.0d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
  end if

  !     This is the decomposition loop of the tridiagonal algorithm
  !     Y2 and U are used for temporary storage of the decomposed factors
  do i=2,n-1
    sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
    p=sig*y2(i-1)+2.0d0
    y2(i)=(sig-1.0d0)/p
    u(i)=(6.0d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))&
        &              /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
  end do

  if(ypn.ge.1.0d30) then
    !        The upper boundary condition is either set to be "natural"
    qn=0.0d0
    un=0.0d0
  else
    !        or to have a specified first derivative
    qn=0.5d0
    un=(3.0d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
  end if

  y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.0d0)
  !     This is the back-substitution loop of the tridiagonal algorithm
  do i=n-1,1,-1
    y2(i)=y2(i)*y2(i+1)+u(i)
  end do

  return
end
!=======================================================================

!=======================================================================
!
!     Given X1A, X2A, YA, M, N (as described in SPLIE2) and Y2A (as
!     produced by that routine), and given a desired interpolating
!     point (X1,X2), this routine returns an interpolated function
!     value Y by performing a bicubic spline interpolation.
!
!-----------------------------------------------------------------------
subroutine splin2(x1a,x2a,ya,y2a,m,n,x1,x2,y)

  use definitions
  use healpix_types

  implicit none

  integer(kind=i4b), intent(in) :: m,n
  real(kind=dp), intent(in)     :: x1a(1:m),x2a(1:n),ya(1:m,1:n),y2a(1:m,1:n)
  real(kind=dp), intent(in)     :: x1,x2
  real(kind=dp), intent(out)    :: y

  integer(kind=i4b) :: i,j
  real(kind=dp) :: ytmp(1:n),y2tmp(1:n),yytmp(1:m),yy2tmp(1:m)

  !     Perform M evaluations of the row splines constructed by
  !     SPLIE2 using the one-dimensional spline evaluator SPLINT
  do i=1,m
    do j=1,n
      ytmp(j)=ya(i,j)
      y2tmp(j)=y2a(i,j)
    end do
    call splint(x2a,ytmp,y2tmp,n,x2,yytmp(i))
  end do

  !     Construct the one-dimensional column spline and evaluate it
  !     Values of 1.0D30 indicate a natural spline
  call spline(x1a,yytmp,m,1.0d30,1.0d30,yy2tmp)
  call splint(x1a,yytmp,yy2tmp,m,x1,y)

  return
end
!=======================================================================

!=======================================================================
!
!     Perform a cubic spline interpolation evaluated at the point X
!     (c.f. Numerical Recipes, Chapter 3.3: Splint Routine,
!           Numerical Recipes, Chapter 3.4: Hunt Routine)
!
!     Given the arrays XA and YA (size N) containing a tabulated
!     function, i.e., YA(I) = f(XA(I)), with the XA(I)'s in order,
!     and given the array Y2A produced by the SPLINE routine, this
!     routine returns a cubic spline interpolated value Y.
!
!     I/O parameters:
!     Input   XA  = vector for independent variable; dimension XA(1:N)
!     Input   YA  = vector for x-dependent variable; dimension YA(1:N)
!     Input   Y2A = 2nd derivative of the function; dimension Y2A(1:N)
!     Input   N   = dimension of input vectors
!     Input   X   = x-value at which Y is to be interpolated
!     Output  Y   = result of interpolation
!
!-----------------------------------------------------------------------
subroutine splint(xa,ya,y2a,n,x,y)

  use definitions
  use healpix_types

  implicit none

  integer(kind=i4b), intent(in) :: n
  real(kind=dp), intent(in)     :: xa(1:n),ya(1:n),y2a(1:n)
  real(kind=dp), intent(in)     :: x
  real(kind=dp), intent(out)    :: y

  logical :: ascnd
  integer(kind=i4b) :: jlo,jhi,jmid,inc
  real(kind=dp) :: a,b

  jlo=0
  jhi=0

  !     ASCND is TRUE if the table values are in ascending order, FALSE otherwise
  ascnd=xa(n).gt.xa(1)

  !     Find the interval XA(JLO) <= X <= XA(JLO+1) = XA(JHI)
  if(jlo.le.0 .or. jlo.gt.n) then
    !        Input guess not useful, go immediately to bisection
    jlo=0
    jhi=n+1
    goto 300
  end if

  !     Set the hunting increment
  inc=1

  if(x.ge.xa(jlo) .eqv. ascnd) then
    !        Hunt up:
    100     jhi=jlo+inc
    if(jhi.gt.n) then
      !           Done hunting, since off the end of the table
      jhi=n+1
    else if(x.ge.xa(jhi) .eqv. ascnd) then
      !           Not done hunting...
      jlo=jhi
      !           ...so double the increment...
      inc=inc+inc
      !           ...and try again
      goto 100
    end if
    !     Done hunting, value bracketed
  else
    jhi=jlo
    !        Hunt down:
    200     jlo=jhi-inc
    if(jlo.lt.1) then
      !           Done hunting, since off the end of the table
      jlo=0
    else if(x.lt.xa(jlo) .eqv. ascnd) then
      !           Not done hunting...
      jhi=jlo
      !           ...so double the increment...
      inc=inc+inc
      !           ...and try again
      goto 200
    end if
    !     Done hunting, value bracketed
  end if

  300  if((jhi-jlo).ne.1) then
    !        Hunt is done, so begin the final bisection phase
    jmid=(jhi+jlo)/2
    if(x.gt.xa(jmid) .eqv. ascnd) then
      jlo=jmid
    else
      jhi=jmid
    end if
    goto 300
  end if

  if(jlo.eq.0) then
    jlo=1
    jhi=2
  end if
  if(jlo.eq.n) then
    jlo=n-1
    jhi=n
  end if

  !     JLO and JHI now bracket the input value X
  !     The cubic spline polynomial is now evaluated
  a=(xa(jhi)-x)/(xa(jhi)-xa(jlo))
  b=(x-xa(jlo))/(xa(jhi)-xa(jlo))
  y=a*ya(jlo)+b*ya(jhi)+((a**3-a)*y2a(jlo)+(b**3-b)*y2a(jhi))&
      &  *((xa(jhi)-xa(jlo))**2)/6.0d0

  return
end
!=======================================================================
