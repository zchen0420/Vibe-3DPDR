module heating_rates_module
  use definitions, only : dp
  use global_module
  use heating_rate_kernels_module
  use healpix_types
  use maincode_module, only : runtime
  use reaction_rates_module, only : reaction_rate_indices
  implicit none

  private
  public :: heating_rate_environment
  public :: calculate_heating_rates

  type :: heating_rate_environment
    real(kind=dp) :: density
    real(kind=dp) :: gas_temperature
    real(kind=dp) :: dust_temperature
    real(kind=dp) :: uv_field
    real(kind=dp) :: turbulent_velocity
    real(kind=dp), pointer :: abundance(:) => null()
    real(kind=dp), pointer :: reaction_rate(:) => null()
  end type heating_rate_environment

contains

  subroutine calculate_heating_rates(environment, reaction_indices, heating_rate)
    type(heating_rate_environment), intent(in) :: environment
    type(reaction_rate_indices), intent(in) :: reaction_indices
    real(kind=dp), intent(out) :: heating_rate(1:12)

    real(kind=dp) :: total_heating,photoelectric_heating,pahphotoelec_heating,weingartner_heating, &
        & cionization_heating,h2formation_heating,h2photodiss_heating,fuvpumping_heating, &
        & cosmicray_heating,turbulent_heating,chemical_heating,gasgrain_heating,softxray_heating

    !     Dust PE heating
    integer(kind=i4b) :: iteration
    real(kind=dp) :: habing_field
    real(kind=dp) :: x,xx,xk,xd,gamma,delta
    real(kind=dp) :: deltad,deltauv,y,hnud,hnuh

    !     PAH heating/cooling
    real(kind=dp) :: epsilon,alpha,beta,pah_heating,pah_cooling
    real(kind=dp) :: phi_pah

    !     Weingartner & Draine treatment of photoelectric heating (grains+PAHs)
    real(kind=dp) :: c0,c1,c2,c3,c4,c5,c6

    !     H2* FUV pumping heating
    real(kind=dp) :: ncr_h2

    !!     Soft X-ray heating
    !      REAL(KIND=DP) :: PP1,PP2,F6

    associate(density => environment%density, gas_temperature => environment%gas_temperature, &
          &dust_temperature => environment%dust_temperature, uv_field => environment%uv_field, &
          &v_turb => environment%turbulent_velocity, abundance => environment%abundance, &
          &rate => environment%reaction_rate, nrgr => reaction_indices%grain_surface, &
          &nrh2 => reaction_indices%h2_photodissociation, nrhd => reaction_indices%hd_photodissociation, &
          &nrco => reaction_indices%co_photodissociation, nrci => reaction_indices%carbon_photoionization, &
          &nrsi => reaction_indices%silicon_photoionization)
      !     Convert the FUV field (in Draine units) to the Habing equivalent
      habing_field=1.68d0*uv_field

      !-----------------------------------------------------------------------
      !     Dust photoelectric heating
      !
      !     Use the treatment of Tielens & Hollenbach, 1985, ApJ, 291, 722,
      !     which follows de Jong (1977,1980)
      !
      !     The charge of a dust grain can be found by equating the rate of
      !     photo-ejection of electrons from the dust grain to the rate of
      !     recombination of electrons with the dust grain (Spitzer)
      !
      !     The various parameter values are taken from Table 2 of the paper
      !-----------------------------------------------------------------------

      deltad=1.0d0
      deltauv=1.8d0
      y=0.1d0
      hnud=6.0d0
      hnuh=13.6d0

      xk=kb*gas_temperature/(hnuh*ev)
      xd=hnud/hnuh
      gamma=2.9d-4*y*sqrt(gas_temperature)*habing_field/(abundance(nelect)*density)
      delta=xk-xd+gamma

      !     Iterate to determine X by finding the zero of the function F
      x=0.5d0
      do iteration=1,100
        xx=x-(f(x,delta,gamma)/ff(x,delta))
        if(abs(xx-x).lt.1.0d-2) exit
        x=xx
      end do
      x=xx

      if(iteration.ge.100) then
        write(10,*)'WARNING! Grain parameter X not found in PE heating'
        write(10,*)'Using final value from interation loop: X =',x
      end if

      !     Assume the dust PE heating scales linearly with metallicity
      photoelectric_heating=2.7d-25*deltauv*deltad*density*y*habing_field &
          & *(((1.0d0-x)**2)/x + xk*((x**2)-1.0d0)/(x**2))*metallicity


      !==================PAH PHOTOELECTRIC HEATING===============================
      !-----------------------------------------------------------------------
      !  Grain + PAH photoelectric heating (MRN size distribution; r = 3-100
      !  Å)
      !
      !  Use the treatment of Bakes & Tielens (1994, ApJ, 427, 822) with the
      !  modifications suggested by Wolfire et al. (2003, ApJ, 587, 278) to
      !  account for the revised PAH abundance estimate from Spitzer data.
      !
      !  See also:
      !  Wolfire et al. (1995, ApJ, 443, 152)
      !  Le Page, Snow & Bierbaum (2001, ApJS, 132, 233)
      !-----------------------------------------------------------------------

      !  Adopt the PAH rate scaling factor of Wolfire et al. (2008, ApJ, 680,
      !  384)
      !  Setting this factor to 1.0 gives the standard Bakes & Tielens
      !  expression
      phi_pah=1.0d0 !0.4D0

      alpha=0.944d0
      beta=0.735d0/gas_temperature**0.068
      delta=habing_field*sqrt(gas_temperature)/(abundance(nelect)*density*phi_pah)
      epsilon=4.87d-2/(1.0d0+4.0d-3*delta**0.73) + 3.65d-2*(gas_temperature/1.0d4)**0.7/(1.0d0+2.0d-4*delta)

      pah_heating=1.30d-24*epsilon*habing_field*density
      pah_cooling=4.65d-30*gas_temperature**alpha*(delta**beta)*abundance(nelect)*density*phi_pah*density

      !  Assume the PE heating rate scales linearly with metallicity
      pahphotoelec_heating=(pah_heating-pah_cooling)*metallicity
      !
      !!==============================================================================
      !

      !-----------------------------------------------------------------------
      !     Weingartner & Draine, 2001, ApJS, 134, 263
      !
      !     Includes photoelectric heating due to PAHs, VSGs and larger grains
      !     Assumes a gas-to-dust mass ratio of 100:1
      !-----------------------------------------------------------------------

      c0=5.72d+0
      c1=3.45d-2
      c2=7.08d-3
      c3=1.98d-2
      c4=4.95d-1
      c5=6.92d-1
      c6=5.20d-1

      weingartner_heating=metallicity*1.0d-26*(habing_field*density)*(c0+c1*gas_temperature**c4) &
          & /(1.0d0+c2*(habing_field*sqrt(gas_temperature)/(abundance(nelect)*density))**c5  &
          & *(1.0d0+c3*(habing_field*sqrt(gas_temperature)/(abundance(nelect)*density))**c6))

      !-----------------------------------------------------------------------
      !     Carbon photoionization heating
      !
      !     1 eV on average per carbon ionization
      !     Use the C photoionization rate determined by calculate_reaction_rates.
      !-----------------------------------------------------------------------
      cionization_heating=(1.0*ev)*rate(nrci)*abundance(species_idx%nc)*density
      !-----------------------------------------------------------------------
      !     H2 formation heating
      !
      !     Assume 1.5 eV liberated as heat during H2 formation
      !     See: Hollenbach & Tielens, Review of Modern Physics, 1999, 71, 173
      !     Use the grain-surface rate determined by calculate_reaction_rates.
      !-----------------------------------------------------------------------

      h2formation_heating=(1.5*ev)*rate(nrgr)*density*abundance(species_idx%nh)*density

      !-----------------------------------------------------------------------
      !     H2 photodissociation heating
      !
      !     0.4 eV on average per photodissociated molecule
      !     Use the H2 photodissociation rate determined by calculate_reaction_rates.
      !-----------------------------------------------------------------------

      h2photodiss_heating=(0.4*ev)*rate(nrh2)*abundance(species_idx%nh2)*density

      !-----------------------------------------------------------------------
      !     H2 FUV pumping heating
      !
      !     2.2 eV on average per vibrationally excited H2* molecule
      !     See: Hollenbach & McKee (1979)
      !     Use the H2 photodissociation rate determined by calculate_reaction_rates.
      !     Use the H2 critical density calculation from Hollenbach & McKee (1979)
      !-----------------------------------------------------------------------

      ncr_h2=1.0d6/sqrt(gas_temperature)/(1.6d0*abundance(species_idx%nh)*exp(-((400.0d0/gas_temperature)**2)) &
          & + 1.4d0*abundance(species_idx%nh2)*exp(-(18100.0d0/(gas_temperature+1200.0d0))))

      fuvpumping_heating=(2.2*ev)*9.0d0*rate(nrh2)*abundance(species_idx%nh2)*density/(1.0d0+ncr_h2/density)

      !-----------------------------------------------------------------------
      !     Cosmic-ray ionization heating
      !
      !     8.0 eV of heat deposited per primary ionization (plus some from He ionization)
      !     Use the treatment of Tielens & Hollenbach, 1985, ApJ, 291, 772
      !     See also: Shull & Van Steenberg, 1985, ApJ, 298, 268
      !               Clavel et al. (1978), Kamp & van Zadelhoff (2001)
      !-----------------------------------------------------------------------

      !      COSMICRAY_HEATING=(20.0*EV)*(1.3D-17*runtime%cosmic_ray_ionization_rate)*DENSITY*ABUNDANCE(species_idx%NH2) !20.0 -> 9.4 eV
      cosmicray_heating=cosmic_ray_heating_rate(runtime%cosmic_ray_ionization_rate,density,abundance(species_idx%nh2))
      !      COSMICRAY_HEATING=(9.4*EV)*(1.3D-17*runtime%cosmic_ray_ionization_rate)*DENSITY

      !-----------------------------------------------------------------------
      !     Supersonic turbulent decay heating
      !
      !     Most relevant for the inner parsecs of galaxies (Black)
      !     Black, in Interstellar Processes, 1987, p731
      !     See also: Rodriguez-Fernandez et al., 2001, A&A, 365, 174
      !
      !     V_TURB = turbulent velocity (km/s); Galactic center ~ 15 km/s
      !     L_TURB = turbulent scale length (pc); typically 5 pc
      !-----------------------------------------------------------------------

      turbulent_heating=turbulent_heating_rate(v_turb,density)

      !-----------------------------------------------------------------------
      !     Exothermic chemical reaction heating
      !
      !     See: Clavel et al., 1978, A&A, 65, 435
      !     Recombination reactions: HCO+ (7.51 eV); H3+ (4.76+9.23 eV); H3O+ (1.16+5.63+6.27 eV)
      !     Ion-neutral reactions  : He+ + H2 (6.51 eV); He+ + CO (2.22 eV)
      !     For each reaction, the heating rate should be: n(1) * n(2) * K * E
      !     with n(1) and n(2) the densities, K the rate coefficient [cm^3.s^-1], and E the energy [erg]
      !-----------------------------------------------------------------------

#ifdef REDUCED
      chemical_heating=abundance(species_idx%nh2x)*density*abundance(nelect)*rate(216)*10.9*ev& !H2+ + e-
          & + abundance(species_idx%nh2x)*density*abundance(species_idx%nh)*rate(155)*0.94*ev& !H2+ + H
          & + abundance(species_idx%nhcox)*density*abundance(nelect)*density*(rate(240)*(7.51*ev)) & ! HCO+ + e-
          & + abundance(species_idx%nh3x)*density*abundance(nelect)*density*(rate(217)*(4.76*ev)+rate(218)*(9.23*ev)) & ! H3+  + e-
          & + abundance(species_idx%nh3ox)*density*abundance(nelect)*density*(rate(236)*(1.16*ev)+&
          &rate(237)*(5.63*ev)+rate(238)*(6.27*ev)) & ! H3O+ + e-
          & + abundance(species_idx%nhex)*density*abundance(species_idx%nh2)*density*(rate(50)*(6.51*ev)+rate(170)*(6.51*ev)) & ! He+  + H2
          & + abundance(species_idx%nhex)*density*abundance(species_idx%nco)*density*(rate(89)*(2.22*ev)+rate(90)*(2.22*ev)+&
          &rate(91)*(2.22*ev)) ! He+  + CO
#endif

#ifdef FULL
      chemical_heating=abundance(species_idx%nhcox)*density*abundance(nelect)*density*(rate(730)*(7.51*ev)) & ! HCO+ + e-
          & + abundance(species_idx%nh3x)*density*abundance(nelect)*density*(rate(706)*(4.76*ev)+rate(705)*(9.23*ev)) & ! H3+  + e-
          & + abundance(species_idx%nh3ox)*density*abundance(nelect)*density*(rate(717)*(1.16*ev)+&
          &rate(716)*(5.63*ev)+rate(714)*(6.27*ev)) & ! H3O+ + e-
          & + abundance(species_idx%nhex)*density*abundance(species_idx%nh2)*density*(rate(1227)*(6.51*ev)+rate(265)*(6.51*ev)) & ! He+  + H2
          & + abundance(species_idx%nhex)*density*abundance(species_idx%nco)*density*(rate(1541)*(2.22*ev))
#endif

#ifdef MYNETWORK
      stop "CHEMICAL_HEATING function has to be declared at &
          & [heating_rates.F90] &
          If you are using the pre-set 'mynetwork' network comment &
          & out this STOP [heating_rates.F90]"
      chemical_heating=abundance(species_idx%nhcox)*density*abundance(nelect)*density*(rate(185)*(7.51*ev)) & ! HCO+ + e-
          & + abundance(species_idx%nh3x)*density*abundance(nelect)*density*(rate(173)*(4.76*ev)+rate(172)*(9.23*ev)) & ! H3+  + e-
          & + abundance(species_idx%nh3ox)*density*abundance(nelect)*density*(rate(179)*(1.16*ev)+&
          &rate(178)*(5.63*ev)+rate(176)*(6.27*ev)) & ! H3O+ + e-
          & + abundance(species_idx%nhex)*density*abundance(species_idx%nh2)*density*(rate(297)*(6.51*ev)+rate(67)*(6.51*ev)) & ! He+  + H2
          & + abundance(species_idx%nhex)*density*abundance(species_idx%nco)*density*(rate(364)*(2.22*ev))
#endif

      !-----------------------------------------------------------------------
      !     Gas-grain collisional heating
      !
      !     Use the treatment of Burke & Hollenbach, 1983, ApJ, 265, 223, and
      !     accommodation fitting formula of Groenewegen, 1994, A&A, 290, 531
      !
      !     Other relevant references:
      !     Hollenbach & McKee, 1979, ApJS, 41,555
      !     Tielens & Hollenbach, 1985, ApJ, 291,722
      !     Goldsmith, 2001, ApJ, 557, 736
      !
      !     This process is insignificant for the energy balance of the dust
      !     but can influence the gas temperature. If the dust temperature is
      !     lower than the gas temperature, this becomes a cooling mechanism
      !
      !     In Burke & Hollenbach (1983) the factor:
      !
      !     (8*kb/(pi*mass_proton))**0.5*2*kb = 4.003D-12
      !
      !     This value has been used in the expression below
      !-----------------------------------------------------------------------
      !
      gasgrain_heating=gas_grain_exchange_rate(density,gas_temperature,dust_temperature)

      ! Gas-grain collisional heating (Tielens 2005 ISM book)
      !      GASGRAIN_HEATING=-1d-33*DENSITY**2*sqrt(GAS_TEMPERATURE)*(GAS_TEMPERATURE-DUST_TEMPERATURE)


      !-----------------------------------------------------------------------
      !     Soft X-ray heating
      !
      !     Use the treatment of from Wolfire et al., 1995, ApJ, 443, 152
      !-----------------------------------------------------------------------

      !      IF(DEPTH.EQ.0) THEN
      !         PP1=0.0D0
      !      ELSE
      !         PP1=DLOG10(DENSITY*(DIST(DEPTH)-DIST(DEPTH-1))/1.0D18)
      !      ENDIF
      !
      !      PP2=DLOG10(ABUNDANCE(NELECT))
      !      F6=0.990D0-2.74D-3*PP2+1.13D-3*PP2**2
      !
      !      SOFTXRAY_HEATING=10.0D0**(F6*(-26.5D0-0.920D0*PP1+5.89D-2*PP1**2)
      !     *                +F6*0.96D0*EXP(-(((PP1-0.38D0)/0.87D0)**2)))



      !-----------------------------------------------------------------------
      !     Total heating rate (sum of all contributions)
      !-----------------------------------------------------------------------

      total_heating=&
      !           & + PHOTOELECTRIC_HEATING &
          & + pahphotoelec_heating &
      !           & + WEINGARTNER_HEATING &
          & + cionization_heating &
          & + h2formation_heating &
          & + h2photodiss_heating &
          & + fuvpumping_heating &
          & + cosmicray_heating &
          & + turbulent_heating &
          & + chemical_heating &
      !           & + SOFTXRAY_HEATING &
          & + gasgrain_heating

      call store_heating_rates(heating_rate,photoelectric_heating,pahphotoelec_heating, &
          & weingartner_heating,cionization_heating,h2formation_heating,h2photodiss_heating, &
          & fuvpumping_heating,cosmicray_heating,turbulent_heating,chemical_heating, &
          & gasgrain_heating,total_heating)

      !-----------------------------------------------------------------------

    end associate

  contains ! Dust photoelectric heating functions...

    !=======================================================================
    !     X is the grain charge parameter and is the solution to F(X)=0
    !-----------------------------------------------------------------------
    function f(x,delta,gamma)

      use definitions
      use healpix_types

      implicit none

      real(kind=dp) :: f
      real(kind=dp), intent(in) :: x,delta,gamma

      f=(x**3)+delta*(x**2)-gamma

    end function f
    !-----------------------------------------------------------------------

    !=======================================================================
    !     FF(X) is the derivative of F(X) with respect to X
    !-----------------------------------------------------------------------
    function ff(x,delta)

      use definitions
      use healpix_types

      implicit none

      real(kind=dp) :: ff
      real(kind=dp), intent(in) :: x,delta

      ff=3*(x**2)+delta*(2*x)

    end function ff
    !-----------------------------------------------------------------------

  end subroutine calculate_heating_rates
  !=======================================================================

end module heating_rates_module
