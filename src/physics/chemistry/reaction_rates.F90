module reaction_rates_module
  use chemistry_network_module, only : reaction_network
  use definitions, only : dp
  use global_module
  use healpix_types, only : i4b
  use maincode_module, only : runtime
  use photo_rate_interfaces_module
  use reaction_rate_kernels_module
  implicit none

  private
  public :: reaction_rate_environment
  public :: reaction_rate_indices
  public :: calculate_reaction_rates

  type :: reaction_rate_environment
    real(kind=dp) :: gas_temperature
    real(kind=dp) :: dust_temperature
    real(kind=dp), pointer :: radiation_surface(:) => null()
    real(kind=dp), pointer :: visual_extinction(:) => null()
    real(kind=dp), pointer :: column_density(:,:) => null()
  end type reaction_rate_environment

  type :: reaction_rate_indices
    integer(kind=i4b) :: grain_surface
    integer(kind=i4b) :: h2_photodissociation
    integer(kind=i4b) :: hd_photodissociation
    integer(kind=i4b) :: co_photodissociation
    integer(kind=i4b) :: carbon_photoionization
    integer(kind=i4b) :: silicon_photoionization
  end type reaction_rate_indices

contains

  subroutine calculate_reaction_rates(environment, network, rate, indices)
    type(reaction_rate_environment), intent(in) :: environment
    type(reaction_network), intent(in) :: network
    real(kind=dp), intent(out) :: rate(:)
    type(reaction_rate_indices), intent(out) :: indices

  real(kind=dp) :: temperature
  real(kind=dp) :: dust_temperature
  real(kind=dp) :: PHI_PAH,FLUX,YIELD
  integer(kind=i4b) :: I,J,K
  integer(kind=i4b) :: nreac
  integer(kind=i4b) :: ray_index_lower
  integer(kind=i4b) :: ray_index_upper

  temperature = environment%gas_temperature
  dust_temperature = environment%dust_temperature
  nreac = size(rate)
  ray_index_lower = lbound(environment%radiation_surface,1)
  ray_index_upper = ubound(environment%radiation_surface,1)

  associate(reactant => network%reactant, product => network%product, alpha => network%alpha, &
      &beta => network%beta, gamma => network%gamma, rtmin => network%rtmin, rtmax => network%rtmax, &
      &duplicate => network%duplicate, rad_surface => environment%radiation_surface, &
      &av => environment%visual_extinction, column => environment%column_density)

  !     Initialize the RATE coefficients.
  RATE=0.0D0

  !     Initialize the stored reaction numbers. If they are not assigned
  !     subsequently, any attempt to access that reaction will generate an
  !     error and the code will crash. This is a useful bug catch.
  indices%grain_surface=0
  indices%h2_photodissociation=0
  indices%hd_photodissociation=0
  indices%co_photodissociation=0
  indices%carbon_photoionization=0
  indices%silicon_photoionization=0

  DO I=1,NREAC
    !        Determine the type of reaction
    IF(REACTANT(I,2).EQ."PHOTON") GOTO 1
    IF(REACTANT(I,2).EQ."CRP   ") GOTO 2
    IF(REACTANT(I,2).EQ."CRPHOT") GOTO 3
    IF(REACTANT(I,2).EQ."FREEZE") GOTO 4
    IF(REACTANT(I,2).EQ."ELFRZE") GOTO 5
    IF(REACTANT(I,2).EQ."CRH   ") GOTO 6
    IF(REACTANT(I,2).EQ."PHOTD ") GOTO 7
    IF(REACTANT(I,2).EQ."THERM ") GOTO 8
    IF(REACTANT(I,2)(1:1).EQ."#") GOTO 9
    IF(REACTANT(I,2).EQ."XRAY  ") GOTO 110
    IF(REACTANT(I,2).EQ."XRSEC ") GOTO 120
    IF(REACTANT(I,2).EQ."XRLYA ") GOTO 130
    IF(REACTANT(I,2).EQ."XRPHOT") GOTO 140

    !-----------------------------------------------------------------------

    !     Thermal reactions:

    !     The RATE of H2 formation on grains is calculated separately
    !     by the function H2_FORMATION_RATE (see function for details)
    IF((REACTANT(I,1).EQ."H  " .AND. REACTANT(I,2).EQ."H  "   .AND. &
        & (REACTANT(I,3).EQ."   " .OR.  REACTANT(I,3).EQ."#  ")) .AND. &
        &  (PRODUCT(I,1).EQ."H2 " .AND. &
        &  (PRODUCT(I,2).EQ."   " .OR.  PRODUCT(I,2).EQ."#  "))) THEN
#ifdef H2FORM
    RATE(I)=H2_FORMATION_RATE(TEMPERATURE,DUST_TEMPERATURE)
#else
    !            RATE(I)=3.0D-18*SQRT(TEMPERATURE)*EXP(-(TEMPERATURE/1.0D3))
    RATE(I)=3.0D-18*SQRT(TEMPERATURE)
#endif
    indices%grain_surface=I
    GOTO 10
  ENDIF

  !     Rates for reactions involving PAHs are calculated according to the
  !     treatment of Wolfire et al. (2003, ApJ, 587, 278; 2008, ApJ, 680, 384)
  IF(ANY(REACTANT(I,:).EQ."PAH  ") .OR. ANY(REACTANT(I,:).EQ."PAH0 ") .OR. &
      & ANY(REACTANT(I,:).EQ."PAH+ ") .OR. ANY(REACTANT(I,:).EQ."PAH- ")) THEN
  PHI_PAH=0.4D0
  RATE(I)=ALPHA(I)*(TEMPERATURE/100.0D0)**BETA(I)*PHI_PAH
  GOTO 10
END IF

!C     Check for large negative GAMMA values that might cause discrepant
!C     rates at low temperatures. Set these rates to zero when T < RTMIN.

!CODE RESPONSINBLE FOR O + H+ --> O+ + H
IF(DUPLICATE(I).EQ.0) THEN
  IF(GAMMA(I).LT.-200.0D0 .AND. TEMPERATURE.LT.RTMIN(I)) THEN
    RATE(I)=0.0D0
  ELSE
    !if (i==76) then
    !RATE(i)=ALPHA(i)
    !write(6,*) '1'
    !else
      RATE(I)=arrhenius_reaction_rate(ALPHA(I),BETA(I),GAMMA(I),TEMPERATURE,RTMIN(I))
    !endif
  ENDIF
ELSE IF(DUPLICATE(I).EQ.1) THEN
  J=I
  DO
    IF(TEMPERATURE.LE.RTMAX(J)) THEN
      IF(GAMMA(J).LT.-200.0D0 .AND. TEMPERATURE.LT.RTMIN(J)) THEN
        RATE(J)=0.0D0
      ELSE
        !if (i==76) then
        !RATE(j)=ALPHA(j)
        !write(6,*) '2'
        !else
        RATE(J)=arrhenius_reaction_rate(ALPHA(J),BETA(J),GAMMA(J),TEMPERATURE,RTMIN(J))
        !endif
      ENDIF
      EXIT
    ELSE IF(DUPLICATE(J+1).LT.DUPLICATE(J)) THEN
      IF(GAMMA(J).LT.-200.0D0 .AND. TEMPERATURE.LT.RTMIN(J)) THEN
        RATE(J)=0.0D0
      ELSE
        !if (i==76) then
        !RATE(j)=ALPHA(j)
        !write(6,*) '3'
        !else
        RATE(J)=arrhenius_reaction_rate(ALPHA(J),BETA(J),GAMMA(J),TEMPERATURE,RTMIN(J))
        !endif
      ENDIF
      EXIT
    ELSE
      RATE(J)=0.0D0
      J=J+1
    ENDIF
  ENDDO
ENDIF
GOTO 10

!C-----------------------------------------------------------------------

!C     Photoreactions:

!C     Store the reaction number for H2 photodissociation. The RATE itself
!C     is calculated separately by the function H2PDRATE (within shield.f)
1       IF(REACTANT(I,1).EQ."H2 " .AND. REACTANT(I,3).EQ."   ") THEN
!C           Loop over all rays
do k=ray_index_lower,ray_index_upper
  RATE(I)=RATE(I) + H2PDRATE(ALPHA(I),RAD_SURFACE(K),AV(K),COLUMN(K,species_idx%NH2))
ENDDO
IF(PRODUCT(I,1).EQ."H " .AND. PRODUCT(I,2).EQ."H ") indices%h2_photodissociation=I
GOTO 10
ENDIF

!C     Store the reaction number for HD photodissociation. The RATE itself
!C     is calculated separately by the function H2PDRATE (within shield.f)
IF(REACTANT(I,1).EQ."HD " .AND. REACTANT(I,3).EQ."   ") THEN
  !C           Loop over all rays
  do k=ray_index_lower,ray_index_upper
    RATE(I)=RATE(I) + H2PDRATE(ALPHA(I),RAD_SURFACE(K),AV(K),COLUMN(K,species_idx%NHD))
  ENDDO
  IF(ANY(PRODUCT(I,:).EQ."H ") .AND. ANY(PRODUCT(I,:).EQ."D ")) indices%hd_photodissociation=I
  GOTO 10
ENDIF

!C     Store the reaction number for !CO photodissociation. The RATE itself
!C     is calculated separately by the function !COPDRATE (within shield.f)
IF(REACTANT(I,1).EQ."CO " .AND. REACTANT(I,3).EQ."   " .AND. &
    & ANY(PRODUCT(I,:).EQ."C ") .AND. ANY(PRODUCT(I,:).EQ."O ")) THEN
!C           Loop over all rays
do k=ray_index_lower,ray_index_upper
  RATE(I)=RATE(I) + COPDRATE(ALPHA(I),RAD_SURFACE(K),AV(K),COLUMN(K,species_idx%NCO),COLUMN(K,species_idx%NH2))
ENDDO
indices%co_photodissociation=I
GOTO 10

ENDIF

!C     Store the reaction number for !CI photoionization. The RATE itself
!C     is calculated separately by the function CIPDRATE (within shield.f)
IF(REACTANT(I,1).EQ."C  " .AND. REACTANT(I,3).EQ."   " .AND.&
    &    ((PRODUCT(I,1).EQ."C+ " .AND. PRODUCT(I,2).EQ."e- ") .OR.&
    &     (PRODUCT(I,1).EQ."e- " .AND. PRODUCT(I,2).EQ."C+ "))) THEN
!C           Loop over all rays
do k=ray_index_lower,ray_index_upper
  RATE(I)=RATE(I) + CIPDRATE(ALPHA(I),RAD_SURFACE(K),AV(K),GAMMA(I),COLUMN(K,species_idx%NC),COLUMN(K,species_idx%NH2),TEMPERATURE)
ENDDO
indices%carbon_photoionization=I
GOTO 10

ENDIF

!C     Store the reaction number for SI photoionization. The RATE itself
!C     is calculated separately by the function SIPDRATE (within shield.f)
IF(REACTANT(I,1).EQ."S  " .AND. REACTANT(I,3).EQ."   " .AND.&
    &    ((PRODUCT(I,1).EQ."S+ " .AND. PRODUCT(I,2).EQ."e- ") .OR.&
    &     (PRODUCT(I,1).EQ."e- " .AND. PRODUCT(I,2).EQ."S+ "))) THEN
!C           Loop over all rays
do k=ray_index_lower,ray_index_upper
  RATE(I)=RATE(I) + SIPDRATE(ALPHA(I),RAD_SURFACE(K),AV(K),GAMMA(I),COLUMN(K,species_idx%NS))
ENDDO
indices%silicon_photoionization=I
GOTO 10
ENDIF

IF(DUPLICATE(I).EQ.0) THEN
  !C           Loop over all rays
  do k=ray_index_lower,ray_index_upper
    RATE(I)=RATE(I) + ALPHA(I)*RAD_SURFACE(K)*EXP(-(GAMMA(I)*AV(K)))/2.0
  ENDDO
ELSE IF(DUPLICATE(I).EQ.1) THEN
  J=I
  DO
    IF(TEMPERATURE.LE.RTMAX(J)) THEN
      !C                 Loop over all rays
      do k=ray_index_lower,ray_index_upper
        RATE(J)=RATE(J) + ALPHA(J)*RAD_SURFACE(K)*EXP(-(GAMMA(J)*AV(K)))/2.0
      ENDDO
      EXIT
    ELSE IF(DUPLICATE(J+1).LT.DUPLICATE(J)) THEN
      !C                 Loop over all rays
      do k=ray_index_lower,ray_index_upper
        RATE(J)=RATE(J) + ALPHA(J)*RAD_SURFACE(K)*EXP(-(GAMMA(J)*AV(K)))/2.0
      ENDDO
      EXIT
    ELSE
      RATE(J)=0.0D0
      J=J+1
    ENDIF
  ENDDO
ENDIF
GOTO 10

!C-----------------------------------------------------------------------

!C     Cosmic ray-induced ionization:

2       IF(DUPLICATE(I).EQ.0) THEN
RATE(I)=ALPHA(I)*runtime%cosmic_ray_ionization_rate
ELSE IF(DUPLICATE(I).EQ.1) THEN
  J=I
  DO
    IF(TEMPERATURE.LE.RTMAX(J)) THEN
      RATE(J)=ALPHA(J)*runtime%cosmic_ray_ionization_rate
      EXIT
    ELSE IF(DUPLICATE(J+1).LT.DUPLICATE(J)) THEN
      RATE(J)=ALPHA(J)*runtime%cosmic_ray_ionization_rate
      EXIT
    ELSE
      RATE(J)=0.0D0
      J=J+1
    ENDIF
  ENDDO
ENDIF
GOTO 10

110     RATE(I)=0.0D0; GOTO 10
120     RATE(I)=0.0D0; GOTO 10
130     RATE(I)=0.0D0; GOTO 10
140     RATE(I)=0.0D0; GOTO 10

!C-----------------------------------------------------------------------

!C     Photoreactions due to cosmic ray-induced secondary photons:

3       IF(DUPLICATE(I).EQ.0) THEN
RATE(I)=ALPHA(I)*runtime%cosmic_ray_ionization_rate*(TEMPERATURE/300.0D0)**BETA(I)&
    &           *GAMMA(I)/(1.0D0-OMEGA)
ELSE IF(DUPLICATE(I).EQ.1) THEN
  J=I
  DO
    IF(TEMPERATURE.LE.RTMAX(J)) THEN
      RATE(J)=ALPHA(J)*runtime%cosmic_ray_ionization_rate*(TEMPERATURE/300.0D0)**BETA(J)&
          &                 *GAMMA(J)/(1.0D0-OMEGA)
      EXIT
    ELSE IF(DUPLICATE(J+1).LT.DUPLICATE(J)) THEN
      RATE(J)=ALPHA(J)*runtime%cosmic_ray_ionization_rate*(TEMPERATURE/300.0D0)**BETA(J)&
          &                 *GAMMA(J)/(1.0D0-OMEGA)
      EXIT
    ELSE
      RATE(J)=0.0D0
      J=J+1
    ENDIF
  ENDDO
ENDIF
GOTO 10

!C-----------------------------------------------------------------------

!C     Freeze-out of neutral chemistry%network%species:

4       RATE(I)=freeze_out_reaction_rate(ALPHA(I),BETA(I),GAMMA(I),TEMPERATURE)
GOTO 10

!C-----------------------------------------------------------------------

!C     Freeze-out of singly charged positive ions:

5       RATE(I)=freeze_out_reaction_rate(ALPHA(I),BETA(I),GAMMA(I),TEMPERATURE)
GOTO 10

!C-----------------------------------------------------------------------

!C     Desorption due to cosmic ray heating:

!CC     Treatment of Hasegawa & Herbst (1993, MNRAS, 261, 83, Equation 15)
!C 6       RATE(I)=ALPHA(I)*runtime%cosmic_ray_ionization_rate

!C     Treatment of Roberts et al. (2007, MNRAS, 382, 773, Equation 3)
6       RATE(I)=cosmic_ray_desorption_rate(GAMMA(I),runtime%cosmic_ray_ionization_rate)
GOTO 10

!C-----------------------------------------------------------------------

!C     Photodesorption:

7       YIELD=photodesorption_yield(TEMPERATURE)
!C         FLUX=1.0D8 ! Flux of FUV photons in the unattenuated Habing field (in photons cm^-2 s^-1)
FLUX=1.7D8 ! Flux of FUV photons in the unattenuated Draine field (in photons cm^-2 s^-1)
!C        Loop over all rays
do k=ray_index_lower,ray_index_upper
  RATE(I)=RATE(I) + FLUX*RAD_SURFACE(K)*EXP(-(1.8D0*AV(K)))*2.4D-22*YIELD
ENDDO
GOTO 10

!C-----------------------------------------------------------------------

!C     Thermal desorption:

!C     Treatment of Hasegawa, Herbst & Leung (1992, ApJS, 82, 167, Equations 2 & 3)
8       RATE(I)=thermal_desorption_rate(ALPHA(I),GAMMA(I),DUST_TEMPERATURE)
GOTO 10

!C-----------------------------------------------------------------------

!C     Grain mantle reactions:

9       RATE(I)=ALPHA(I)
GOTO 10

!C-----------------------------------------------------------------------

!C     Check that the RATE is physical (0<RATE(I)<1) and produce an error
!C     message if not. Impose a lower cut-off on all RATE coefficients to
!C     prevent the problem becoming too stiff. Rates less than 1E-99 are
!C     set to zero. Grain-surface reactions and desorption mechanisms are
!C     allowed rates greater than 1.
10      call clamp_reaction_rate(RATE(I),I,REACTANT(I,1))
!C     End of loop over rates
ENDDO

RETURN
end associate
  end subroutine calculate_reaction_rates

end module reaction_rates_module
