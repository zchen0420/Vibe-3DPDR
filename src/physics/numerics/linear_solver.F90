subroutine gauss_jordan(a,n,np,b,call_writes)

  use definitions
  use healpix_types
  implicit none
  !      integer(kind=i4b), intent(in) :: coolant
  integer(kind=i4b):: i,j,k,l,ll,irow,icol
  integer(kind=i4b), intent(in):: n,np !,M,MPP
  integer(kind=i4b), parameter :: nmax=100
  integer(kind=i4b):: ipiv(1:nmax),indxr(1:nmax),indxc(1:nmax)
  real(kind=dp), intent(inout) :: a(1:np,1:np)
  real(kind=dp), intent(inout) :: b(1:np) !,1:MPP)
  real(kind=dp) :: big,dum,pivinv
  logical,intent(out)::call_writes

  icol=0
  irow=0
  ipiv=0
  do i=1,n
    big=0.0d0
    do j=1,n
      if(ipiv(j).ne.1) then
        do k=1,n
          if(ipiv(k).eq.0) then
            if(abs(a(j,k)).ge.big) then
              big=abs(a(j,k))
              irow=j
              icol=k
            end if
          else if(ipiv(k).gt.1) then
            print *,'ERROR! Singular matrix in GAUSS_JORDAN'
            call_writes=.true.
            return
            !                     write(6,*) 'Crashed in first loop'
            !                     write(6,*) 'grid point = ',p, ' coolant = ',coolant
            !                     write(6,*) 'thermal%gas_temperature = ',thermal%gas_temperature(p)
            !                     STOP
          end if
        end do
      end if
    end do
    ipiv(icol)=ipiv(icol)+1
    if(irow.ne.icol) then
      do l=1,n
        dum=a(irow,l)
        a(irow,l)=a(icol,l)
        a(icol,l)=dum
      end do
      !            DO L=1,M
      !               DUM=B(IROW,L)
      !               B(IROW,L)=B(ICOL,L)
      !               B(ICOL,L)=DUM
      !            ENDDO
      !================================
      dum=b(irow)
      b(irow)=b(icol)
      b(icol)=dum
      !================================
    end if
    indxr(i)=irow
    indxc(i)=icol
    if(a(icol,icol).eq.0.0d0) then
      print *,'ERROR! Singular matrix found by GAUSS_JORDAN'
      call_writes=.true.
      return
      !            write(6,*) 'Crashed in second loop'
      !            write(6,*) 'grid point = ',p, ' coolant = ',coolant
      !            write(6,*) 'thermal%gas_temperature = ',thermal%gas_temperature(p)
      !            STOP
    end if
    pivinv=1.0d0/a(icol,icol)
    a(icol,icol)=1.0d0
    do l=1,n
      a(icol,l)=a(icol,l)*pivinv
    end do
    !         DO L=1,M
    !            B(ICOL,L)=B(ICOL,L)*PIVINV
    !         ENDDO
    !=======================================
    b(icol)=b(icol)*pivinv
    !=======================================
    do ll=1,n
      if(ll.ne.icol) then
        dum=a(ll,icol)
        a(ll,icol)=0.0d0
        do l=1,n
          a(ll,l)=a(ll,l)-a(icol,l)*dum
        end do
        !               DO L=1,M
        !                  B(LL,L)=B(LL,L)-B(ICOL,L)*DUM
        !               ENDDO
        !=============================================
        b(ll)=b(ll)-b(icol)*dum
        !=============================================
      end if
    end do
  end do
  do l=n,1,-1
    if(indxr(l).ne.indxc(l)) then
      do k=1,n
        dum=a(k,indxr(l))
        a(k,indxr(l))=a(k,indxc(l))
        a(k,indxc(l))=dum
      end do
    end if
  end do
  return
end subroutine
