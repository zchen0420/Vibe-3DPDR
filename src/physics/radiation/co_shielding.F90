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
  FUNCTION COPDRATE(K0,G0,AV,NCO,NH2)

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


    LAMBDA=LBAR(NCO,NH2)

    !C     Calculate the CO photodissociation chemistry%rate (COPDRATE)
    COPDRATE=K0*G0*COSHIELD(NCO,NH2)*SCATTER(AV,LAMBDA)/2.0
    RETURN
    END
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
            FUNCTION COSHIELD(NCO,NH2)

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
              use uclpdr_module, only : start, NCO_GRID, NH2_GRID, SCO_GRID, &
                  & DIMCO, SCO_DERIV, DIMH2, SCO_DERIV
              implicit none
              real(kind=dp) :: COSHIELD
              real(kind=dp) :: LOGNCO, LOGNH2
              !     integer(kind=i4b), intent(in) :: NCO, NH2
              real(kind=dp), intent(in) :: nh2, nco

              IF(START) THEN
                CALL SPLIE2(NCO_GRID,NH2_GRID,SCO_GRID,DIMCO,DIMH2,SCO_DERIV)
                START=.FALSE.
              ENDIF

              LOGNCO=DLOG10(NCO+1.0D0)
              LOGNH2=DLOG10(NH2+1.0D0)

              IF(LOGNCO.LT.NCO_GRID(1)) LOGNCO=NCO_GRID(1)
              IF(LOGNH2.LT.NH2_GRID(1)) LOGNH2=NH2_GRID(1)
              IF(LOGNCO.GT.NCO_GRID(DIMCO)) LOGNCO=NCO_GRID(DIMCO)
              IF(LOGNH2.GT.NH2_GRID(DIMH2)) LOGNH2=NH2_GRID(DIMH2)

              CALL SPLIN2(NCO_GRID,NH2_GRID,SCO_GRID,SCO_DERIV,&
                  &            DIMCO,DIMH2,LOGNCO,LOGNH2,COSHIELD)
              COSHIELD=10.0D0**COSHIELD

              RETURN
              END
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
                  FUNCTION LBAR(NCO,NH2)

                    !      IMPLICIT NONE
                    !      real(kind=dp) ::  NCO,NH2
                    !      real(kind=dp) ::  U,W


                    use definitions
                    use healpix_types
                    !     use global_module, only : NCO, NH2
                    implicit none
                    real(kind=dp) :: lbar
                    real(kind=dp) :: U,W
                    !     integer(kind=i4b), intent(in) :: NCO, NH2
                    real(kind=dp), intent(in) :: nh2,nco

                    U=DLOG10(NCO+1.0D0)
                    W=DLOG10(NH2+1.0D0)

                    LBAR=(5675.0D0 - 200.6D0*W) &
                        &    - (571.6D0 - 24.09D0*W)*U &
                        &   + (18.22D0 - 0.7664D0*W)*U**2

                    !C     LBAR cannot be larger than the wavelength of band 33 (1076.1Å)
                    !C     and cannot be smaller than the wavelength of band 1 (913.6Å)
                    IF(LBAR.LT.913.6D0)  LBAR=913.6D0
                    IF(LBAR.GT.1076.1D0) LBAR=1076.1D0

                    RETURN
                    END
                    !C=======================================================================
