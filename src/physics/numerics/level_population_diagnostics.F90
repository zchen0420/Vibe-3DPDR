subroutine gauss_jordan_writes(a,n,np,b,ill)

  use definitions
  use healpix_types
  use maincode_module, only : thermal

  implicit none
  integer(kind=i4b), intent(in) :: ill !,coolant
  integer(kind=i4b):: i,j,k,l,ll,irow,icol
  integer(kind=i4b), intent(in):: n,np !,M,MPP
  integer(kind=i4b), parameter :: nmax=100
  integer(kind=i4b):: ipiv(1:nmax),indxr(1:nmax),indxc(1:nmax)
  real(kind=dp), intent(inout) :: a(1:np,1:np)
  real(kind=dp), intent(inout) :: b(1:np) !,1:MPP)
  real(kind=dp) :: big,dum,pivinv
  integer(kind=i4b) :: diagnostic_point

  diagnostic_point = 1

  write(6,*) 'b'
  do i=1,np
    write(6,*) b(i),i
  end do

  write(6,*) 'a'
  do i=1,np
    do j=1,np
      write(6,*) a(i,j)
    end do
  end do


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
            end if !ABS(A
          else if(ipiv(k).gt.1) then
            print *,'ERROR! Singular matrix in GAUSS_JORDAN'
            write(6,*) 'Crashed in first loop'
            !                     write(6,*) 'grid point = ',p, ' coolant = ',coolant
            write(6,*) 'thermal%gas_temperature = ',thermal%gas_temperature(diagnostic_point)
            stop
          end if !IPIV(K).EQ.0
        end do !K=1,N
      end if !IPIV(J).NE.1
    end do !J=1,N
    ipiv(icol)=ipiv(icol)+1
    if(irow.ne.icol) then
      do l=1,n
        dum=a(irow,l)
        a(irow,l)=a(icol,l)
        a(icol,l)=dum
      end do !L=1,N
      !            DO L=1,M
      !               DUM=B(IROW,L)
      !               B(IROW,L)=B(ICOL,L)
      !               B(ICOL,L)=DUM
      !            ENDDO
      !================================
      dum=b(irow)
      if (i.eq.ill) write(6,*) 'dum=',dum,'A'
      !write(6,*) 'DUM=',DUM
      b(irow)=b(icol)
      if (i.eq.ill) write(6,*) 'b(',irow,')=',b(irow),'B'
      !write(6,*) 'irow=',irow
      !write(6,*) 'B(irow)=',b(irow)
      b(icol)=dum
      if (i.eq.ill) write(6,*) 'b(',icol,')=',b(icol),'C'
      !write(6,*) 'icol=',icol
      !write(6,*) 'b(icol)=',b(icol)
      !================================
    end if !IROW.NE.ICOL
    indxr(i)=irow
    indxc(i)=icol
    if(a(icol,icol).eq.0.0d0) then
      print *,'ERROR! Singular matrix found by GAUSS_JORDAN'
      write(6,*) 'Crashed in second loop'
      !            write(6,*) 'grid point = ',p, ' coolant = ',coolant
      write(6,*) 'thermal%gas_temperature = ',thermal%gas_temperature(diagnostic_point)
      stop
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
    if (i.eq.ill) write(6,*) b(icol),pivinv,'D'
    b(icol)=b(icol)*pivinv
    if (i.eq.ill) write(6,*) b(icol),'E'
    !write(6,*) 'pivinv=',pivinv
    !write(6,*) 'b(icol)=',b(icol)
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
        if (i.eq.ill) then
          write(6,*) 'll=',ll,'F'
          write(6,*) 'b(ll)=',b(ll),'G'
          write(6,*) 'b(icol)=',b(icol),'H'
          write(6,*) 'dum=',dum,'I'
        end if
        b(ll)=b(ll)-b(icol)*dum
        if (i.eq.ill) write(6,*) 'b(ll) after=',b(ll),'J'
        !=============================================
      end if !LL.NE.ICOL
    end do !LL=1,N
  end do ! I=1,N
  do l=n,1,-1
    if(indxr(l).ne.indxc(l)) then
      do k=1,n
        dum=a(k,indxr(l))
        a(k,indxr(l))=a(k,indxc(l))
        a(k,indxc(l))=dum
      end do !K=1,N
    end if !INDXR(L).NE.INDXC(L)
  end do !L=N,1,-1
  do i=1,n
    write(6,*) b(i),i
  end do
  return
end subroutine
