!C=======================================================================
!C
!C     !CO photodissociation chemistry%rate taking into
!C     account shielding and grain extinction
!C
!C-----------------------------------------------------------------------
!C
!C     Input parameters:
!C     K0  = Unattenuated photodissociation chemistry%rate (in cm^3/s)
!C     G0  = Incident FUV field (in Draine units)
!C     AV  = visual extinction (in magnitudes)
!C     N!CO = !CO geometry%column_density density (in cm^-2)
!C     NH2 = H2 geometry%column_density density (in cm^-2)
!C
!C     Program variables:
!C     !COPDRATE = !CO photodissociation chemistry%rate taking into
!C                account self-shielding and grain extinction
!C     LAMBDA   = wavelength (in Å) of a typical transition
!C
!C     Functions called:
!C     LBAR     = function to determine the wavelength
!C     !COSHIELD = !CO shielding function
!C     S!CATTER  = attenuation due to scattering by dust
!C
!C-----------------------------------------------------------------------
function copdrate(k0,g0,av,nco,nh2)

  !      IMPLICIT NONE
  !      real(kind=dp) :: K0,G0,AV,NCO,NH2
  !      real(kind=dp) :: LAMBDA,LBAR
  !      real(kind=dp) :: COSHIELD,SCATTER

  use definitions
  use healpix_types
  !     use global_module, only : nh2
  implicit none
  real(kind=dp) :: copdrate
  real(kind=dp), intent(in) :: k0, g0, av, nco
  !     integer(kind=i4b), intent(in) :: nh2
  real(kind=dp), intent(in) :: nh2
  real(kind=dp) :: lambda, lbar, coshield, scatter


  lambda=lbar(nco,nh2)

  !C     Calculate the CO photodissociation chemistry%rate (COPDRATE)
  copdrate=k0*g0*coshield(nco,nh2)*scatter(av,lambda)/2.0
  return
end
!C=======================================================================
!C
!C     12!CO line shielding, using the computed values listed in
!C     van Dishoeck & Black (1988, ApJ, 334, 771, Table 5)
!C
!C     Appropriate shielding factors are determined by performing a
!C     2-dimensional spline interpolation over the values listed in
!C     Table 5 of van Dishoeck & Black, which include contributions
!C     from self-shielding and H2 screening
!C
!C-----------------------------------------------------------------------
!C
!C     Input parameters:
!C     N!CO = !CO geometry%column_density density (in cm^-2)
!C     NH2 = H2 geometry%column_density density (in cm^-2)
!C
!C     Program variables:
!C     !COSHIELD  = total 12!CO shielding factor containing
!C                 contributions from both H2 and !CO lines
!C                 from 2D spline interpolation over the grid
!C     S!CO_GRID  = log10 values of the 12!CO shielding factors
!C                 from van Dishoeck & Black (1988) as a function
!C                 of !CO geometry%column_density density (1st index) and H2 geometry%column_density
!C                 density (2nd index)
!C     S!CO_DERIV = 2nd derivative of S!CO_GRID values from SPLIE2
!C     N!CO_GRID  = log10 values of !CO geometry%column_density densities (in cm^-2)
!C     NH2_GRID  = log10 values of H2 geometry%column_density densities (in cm^-2)
!C     DIM!CO     = number of !CO geometry%column_density densities
!C     DIMH2     = number of H2 geometry%column_density densities
!C     START     = .TRUE. when !COSHIELD is first called
!C
!C     Functions called:
!C     SPLIE2 =
!C     SPLIN2 =
!C
!C-----------------------------------------------------------------------
function coshield(nco,nh2)

  !      IMPLICIT NONE
  !      LOGICAL :: START
  !      INTEGER(kind=i4b) :: DIMCO,DIMH2
  !      real(kind=dp) ::  NCO,NH2
  !      real(kind=dp) ::  LOGNCO,LOGNH2
  !      real(kind=dp) ::  NCO_GRID(8),NH2_GRID(6)
  !      real(kind=dp) ::  SCO_GRID(8,6),SCO_DERIV(8,6)
  !      COMMON /STATUS/START
  !      COMMON /COGRID/SCO_GRID,SCO_DERIV,NCO_GRID,NH2_GRID,DIMCO,DIMH2

  use definitions
  use healpix_types
  !     use global_module, only : NCO, NH2
  use uclpdr_module, only : start, nco_grid, nh2_grid, sco_grid, &
      & dimco, sco_deriv, dimh2, sco_deriv
  implicit none
  real(kind=dp) :: coshield
  real(kind=dp) :: lognco, lognh2
  !     integer(kind=i4b), intent(in) :: NCO, NH2
  real(kind=dp), intent(in) :: nh2, nco

  if(start) then
    call splie2(nco_grid,nh2_grid,sco_grid,dimco,dimh2,sco_deriv)
    start=.false.
  end if

  lognco=dlog10(nco+1.0d0)
  lognh2=dlog10(nh2+1.0d0)

  if(lognco.lt.nco_grid(1)) lognco=nco_grid(1)
  if(lognh2.lt.nh2_grid(1)) lognh2=nh2_grid(1)
  if(lognco.gt.nco_grid(dimco)) lognco=nco_grid(dimco)
  if(lognh2.gt.nh2_grid(dimh2)) lognh2=nh2_grid(dimh2)

  call splin2(nco_grid,nh2_grid,sco_grid,sco_deriv,&
      &            dimco,dimh2,lognco,lognh2,coshield)
  coshield=10.0d0**coshield

  return
end
!C=======================================================================
!C
!C     !Calculate the mean wavelength (in Å) of the 33 dissociating bands,
!C     weighted by their fractional contribution to the total shielding
!C     van Dishoeck & Black (1988, ApJ, 334, 771, Equation 4)
!C
!C-----------------------------------------------------------------------
!C
!C     Input parameters:
!C     N!CO = !CO geometry%column_density density (in cm^-2)
!C     NH2 = H2 geometry%column_density density (in cm^-2)
!C
!C     Program variables:
!C     LBAR = mean wavelength (in Å)
!C     U    = log10(N!CO)
!C     W    = log10(NH2)
!C
!C-----------------------------------------------------------------------
function lbar(nco,nh2)

  !      IMPLICIT NONE
  !      real(kind=dp) ::  NCO,NH2
  !      real(kind=dp) ::  U,W


  use definitions
  use healpix_types
  !     use global_module, only : NCO, NH2
  implicit none
  real(kind=dp) :: lbar
  real(kind=dp) :: u,w
  !     integer(kind=i4b), intent(in) :: NCO, NH2
  real(kind=dp), intent(in) :: nh2,nco

  u=dlog10(nco+1.0d0)
  w=dlog10(nh2+1.0d0)

  lbar=(5675.0d0 - 200.6d0*w) &
      &    - (571.6d0 - 24.09d0*w)*u &
      &   + (18.22d0 - 0.7664d0*w)*u**2

  !C     LBAR cannot be larger than the wavelength of band 33 (1076.1Å)
  !C     and cannot be smaller than the wavelength of band 1 (913.6Å)
  if(lbar.lt.913.6d0)  lbar=913.6d0
  if(lbar.gt.1076.1d0) lbar=1076.1d0

  return
end
!C=======================================================================
