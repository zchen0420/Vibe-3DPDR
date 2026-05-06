!C***********************************************************************
!C     Read in the chemical reaction rates and the species, masses and
!C     initial abundances (if specified). The  rates and species files
!C     are assumed to have comma separated values (CSV) format. This is
!C     in line with the Rate05 formatting, removing the need for file-
!C     dependent FORMAT statements.
!C***********************************************************************
!
subroutine read_rates(nreac,reac,prod,alpha,beta,gamma,rate,&
      &duplicate,rtmin,rtmax)
  !T.Bell
  use definitions
  use healpix_types
  use global_module

  implicit none
  integer(kind=i4b), intent(in) :: nreac
  integer(kind=i4b), intent(out) :: duplicate(1:nreac)
  real(kind=dp), intent(out) :: alpha(1:nreac),beta(1:nreac),&
      &gamma(1:nreac),rate(1:nreac),rtmin(1:nreac),rtmax(1:nreac)
  character(len=10), intent(out) :: reac(1:nreac,1:3),prod(1:nreac,1:4)

  integer(kind=i4b) :: i,j,n,ratefile
  character(len=1) :: clem
  ratefile = 2

  !C     Initialize the variables and read in the ratefile data. Check that
  !C     the value of NREAC agrees with the number of reactions in the file
  !C     and produce an error message if not.
  reac="          "
  prod="          "
  alpha=0.0d0
  beta=0.0d0
  gamma=0.0d0
  rtmin=0.0d0
  rtmax=0.0d0
  duplicate=0

  rate=0.0d0


#ifdef REDUCED
  open(ratefile,file="data/rates_reduced.d",status="OLD")
#endif
#ifdef FULL
  open(ratefile,file="data/rates_full.d",status="OLD")
#endif
#ifdef MYNETWORK
  open(ratefile,file="data/rates_mynetwork.d",status="OLD")
#endif
  rewind(ratefile)
  do i=1,nreac
    read(ratefile,*,end=1) n,(reac(i,j),j=1,3),(prod(i,j),j=1,4),&
        &                          alpha(i),beta(i),gamma(i), &
        &                          clem,rtmin(i),rtmax(i)
    if(clem.ne."") clem=""

    !C     Check for duplicate reactions and set the DUPLICATE counter to the
    !C     appropriate value. Adjust their minimum temperatures so that the
    !C     temperature ranges are adjacent.
    if(i.gt.1) then
      if(reac(i,1).eq.reac(i-1,1) .and. &
          &     reac(i,2).eq.reac(i-1,2) .and. reac(i,3).eq.reac(i-1,3) .and. &
          &     prod(i,1).eq.prod(i-1,1) .and. prod(i,2).eq.prod(i-1,2) .and. &
          &     prod(i,3).eq.prod(i-1,3) .and. prod(i,4).eq.prod(i-1,4)) then
      if(duplicate(i-1).eq.0) duplicate(i-1)=1
      duplicate(i)=duplicate(i-1)+1
      rtmin(i)=rtmax(i-1)
    else
      duplicate(i)=0
    end if
  else
    duplicate(i)=0
  end if

  !C     Check for negative gamma values as they could cause problems when
  !C     calculating abundances. Produce a warning message if they occur.
  if(gamma(i).lt.0.0d0) then
    write(6,*) 'Negative gamma factor in rate',n
    write(10,"('Negative gamma factor in rate',I5,' (',F8.1,')')")&
        &         n,gamma(i)
  end if
end do
i=i-1
read(ratefile,*,end=1)
i=i+1
1    if(i.ne.nreac) then
  write(6,*) 'ERROR! Number of reactions (NREAC) does not match ', &
      &           'the number of entries in the ratefile'
  stop
end if

close(ratefile)
return
end subroutine
!C-----------------------------------------------------------------------
