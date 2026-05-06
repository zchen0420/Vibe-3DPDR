!C=======================================================================
!C
!C     !CI photoionization chemistry%rate taking into account grain extinction
!C     and shielding by !CI and H2 lines, adopting the treatment of
!C     Kamp & Bertoldi (2000, A&A, 353, 276, Equation 8)
!C
!C-----------------------------------------------------------------------
!C
!C     Input parameters:
!C     K0   = Unattenuated photoionization chemistry%rate (in cm^3/s)
!C     G0   = Incident FUV field (in Draine units)
!C     AV   = visual extinction (in magnitudes)
!C     KAV  = tau(λ)/tau(V) correction factor
!C     N!CI  = !CI geometry%column_density density (in cm^-2)
!C     NH2  = H2 geometry%column_density density (in cm^-2)
!C     TGAS = gas temperature (in K)
!C
!C     Program variables:
!C     !CIPDRATE = !CI photoionization chemistry%rate taking into
!C                account shielding and grain extinction
!C     TAU!C     = optical depth in the !CI absorption band
!C
!C-----------------------------------------------------------------------
function cipdrate(k0,g0,av,kav,nci,nh2,tgas)

  !      IMPLICIT NONE
  !      real(kind=dp) :: K0,G0,AV,KAV,NCI,NH2,TGAS
  !      real(kind=dp) :: TAUC

  use definitions
  use healpix_types
  !     use global_module, only : nh2
  implicit none
  real(kind=dp) :: cipdrate
  real(kind=dp), intent(in) :: k0,g0,av,kav,nci,tgas
  !     integer(kind=i4b), intent(in) :: nh2
  real(kind=dp), intent(in) :: nh2
  real(kind=dp) :: tauc


  !C     !Calculate the optical depth in the !CI absorption band, accounting
  !C     for grain extinction and shielding by !CI and overlapping H2 lines
  tauc=kav*av+1.1d-17*nci+(0.9d0*tgas**0.27d0*(nh2/1.59d21)**0.45d0)

  !C     Calculate the CI photoionization chemistry%rate (CIPDRATE)
  cipdrate=k0*g0*exp(-tauc)/2.0

  return
end
!C=======================================================================

!C=======================================================================
!C
!C     SI photoionization chemistry%rate -- needs to be implemented!
!C     For now, use the standard expression for photorates
!C
!C-----------------------------------------------------------------------
!C
!C     Input parameters:
!C     K0   = Unattenuated photoionization chemistry%rate (in cm^3/s)
!C     G0   = Incident FUV field (in Draine units)
!C     AV   = visual extinction (in magnitudes)
!C     KAV  = tau(λ)/tau(V) correction factor
!C     NSI  = SI geometry%column_density density (in cm^-2)
!C
!C     Program variables:
!C     SIPDRATE = SI photoionization chemistry%rate taking into
!C                account shielding and grain extinction
!C     TAUS     = optical depth in the SI absorption band
!C
!C-----------------------------------------------------------------------
function sipdrate(k0,g0,av,kav) !,NSI)

  !      IMPLICIT NONE
  !      real(kind=dp) K0,G0,AV,KAV,NSI
  !      real(kind=dp) TAUS

  use definitions
  use healpix_types
  implicit none
  real(kind=dp) :: sipdrate
  real(kind=dp), intent(in) :: k0,g0,av,kav !,NSI
  real(kind=dp) :: taus
  !C     Calculate the optical depth in the SI absorption band, accounting
  !C     for grain extinction and shielding by ???
  taus=kav*av

  !C     Calculate the SI photoionization chemistry%rate (SIPDRATE)
  sipdrate=k0*g0*exp(-taus)/2.0

  return
end
