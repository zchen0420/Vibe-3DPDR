module functions_module

  use definitions
  use healpix_types
  implicit none

  interface
    function H2PDRATE(K0,G0,AV,NH2)
      use definitions
      use healpix_types
      real(kind=dp) :: H2PDRATE
      real(kind=dp), intent(in) :: k0, g0, av
      real(kind=dp), intent(in) :: nh2
      real(kind=dp) :: lambda, scatter, h2shield2
    end function H2PDRATE

    function COPDRATE(K0,G0,AV,NCO,NH2)
      use definitions
      use healpix_types
      real(kind=dp) :: copdrate
      real(kind=dp), intent(in) :: k0, g0, av, nco
      real(kind=dp), intent(in) :: nh2
      real(kind=dp) :: lambda, lbar, coshield, scatter
    end function copdrate

    function CIPDRATE(K0,G0,AV,KAV,NCI,NH2,TGAS)
      use definitions
      use healpix_types
      real(kind=dp) :: cipdrate
      real(kind=dp), intent(in) :: K0,G0,AV,KAV,NCI,TGAS
      real(kind=dp), intent(in) :: nh2
      real(kind=dp) :: tauc
    end function cipdrate

    function SIPDRATE(K0,G0,AV,KAV,NSI)
      use definitions
      use healpix_types
      real(kind=dp) :: sipdrate
      real(kind=dp), intent(in) :: K0,G0,AV,KAV,NSI
      real(kind=dp) :: taus
    end function sipdrate

    function H2SHIELD1(NH2,DOPW,RADW)
      use definitions
      use healpix_types
      real(kind=dp) :: h2shield1
      real(kind=dp), intent(in) :: nh2
      real(kind=dp), intent(in) ::DOPW,RADW
      real(kind=dp) :: FPARA, FOSC, TAUD, R, T, U, JD, JR
    end function h2shield1

    function h2shield2(nh2)
      use definitions
      use healpix_types
      use uclpdr_module, only : start, numh2, COL_GRID, SH2_GRID, SH2_DERIV
      real(kind=dp) :: h2shield2
      real(kind=dp), intent(in) :: nh2
    end function h2shield2

    function COSHIELD(NCO,NH2)
      use definitions
      use healpix_types
      use uclpdr_module, only : start, NCO_GRID, NH2_GRID, SCO_GRID, SCO_DERIV
      real(kind=dp) :: COSHIELD
      real(kind=dp) :: LOGNCO, LOGNH2
      real(kind=dp), intent(in) :: NCO, NH2
    end function COSHIELD

    function SCATTER(AV,LAMBDA)
      use definitions
      use healpix_types
      real(kind=dp) :: scatter
      real(kind=dp), intent(in) :: AV, LAMBDA
      real(kind=dp), dimension(0:5), save :: A = (/&
          &1.000D0,2.006D0,-1.438D0,0.7364D0,-0.5076D0,-0.0592D0/)
      real(kind=dp), dimension(0:5), save :: K = (/&
          &0.7514D0,0.8490D0,1.013D0,1.282D0,2.005D0,5.832D0/)
      real(kind=dp) :: EXPONENT, XLAMBDA
    end function scatter

    function XLAMBDA(LAMBDA)
      use definitions
      use healpix_types
      use uclpdr_module, only : start, N_GRID, L_GRID, X_GRID, X_DERIV
      real(kind=dp) :: xlambda
      real(kind=dp), intent(in) :: lambda
    end function xlambda

    function LBAR(NCO,NH2)
      use definitions
      use healpix_types
      real(kind=dp) :: lbar
      real(kind=dp) :: U,W
      real(kind=dp) :: NCO, NH2
    end function LBAR

    function calculate_heating(density, gas_temperature, dust_temperature, UV_field, &
          & v_turb, nspec, dummyabundance, nreac, rate)
      use definitions
      use healpix_types
      real(kind=dp) :: calculate_heating
      integer(kind=i4b) :: nspec, nreac
      real(kind=dp) :: density, gas_temperature, dust_temperature, UV_field, v_turb
      real(kind=dp) :: dummyabundance(1:nspec), rate(1:nreac)
    end function calculate_heating

#ifdef H2FORM

    FUNCTION H2_FORMATION_RATE(GAS_TEMPERATURE,GRAIN_TEMPERATURE) RESULT(RATE)
      USE DEFINITIONS
      USE HEALPIX_TYPES
      IMPLICIT NONE
      REAL(KIND=DP) :: RATE
      REAL(KIND=DP), INTENT(IN) :: GAS_TEMPERATURE,GRAIN_TEMPERATURE
    END FUNCTION H2_FORMATION_RATE

#endif

  end interface

end module functions_module
