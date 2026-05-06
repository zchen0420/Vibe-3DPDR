!=======================================================================
!
!  Calculate the rate of molecular hydrogen (H2) formation on grains
!  using the treatment of Cazaux & Tielens (2002, ApJ, 575, L29) and
!  Cazaux & Tielens (2004, ApJ, 604, 222).
!
!-----------------------------------------------------------------------
function h2_formation_rate(gas_temperature,grain_temperature) result(rate)

  use definitions
  use healpix_types
  use global_module, only: metallicity, g2d
  implicit none

  real(kind=dp) :: rate
  real(kind=dp), intent(in) :: gas_temperature,grain_temperature

  real(kind=dp) :: thermal_velocity,sticking_coefficient,total_cross_section
  real(kind=dp) :: flux,factor1,factor2,epsilon
  real(kind=dp) :: silicate_formation_efficiency,graphite_formation_efficiency
  real(kind=dp) :: silicate_cross_section,silicate_mu,silicate_e_s,silicate_e_h2
  real(kind=dp) :: silicate_e_hp,silicate_e_hc,silicate_nu_h2,silicate_nu_hc
  real(kind=dp) :: graphite_cross_section,graphite_mu,graphite_e_s,graphite_e_h2
  real(kind=dp) :: graphite_e_hp,graphite_e_hc,graphite_nu_h2,graphite_nu_hc

  !  Mean thermal velocity of hydrogen atoms (cm s^-1)
  thermal_velocity=1.45d5*sqrt(gas_temperature/1.0d2)

  !  Calculate the thermally averaged sticking coefficient of hydrogen atoms on grains,
  !  as given by Hollenbach & McKee (1979, ApJS, 41, 555, eqn 3.7)
  sticking_coefficient=1.0d0/(1.0d0+0.04d0*sqrt(gas_temperature+grain_temperature) &
      & + 0.2d0*(gas_temperature/1.0d2)+0.08d0*(gas_temperature/1.0d2)**2)

  flux=1.0d-10 ! Flux of H atoms in monolayers per second (mLy s^-1)

  total_cross_section=6.273d-22 ! Total mixed grain cross section per H nucleus (cm^-2/nucleus)
  silicate_cross_section=8.473d-22 ! Silicate grain cross section per H nucleus (cm^-2/nucleus)
  graphite_cross_section=7.908d-22 ! Graphite grain cross section per H nucleus (cm^-2/nucleus)

  !  Silicate grain properties
  silicate_mu=0.005d0 ! Fraction of newly formed H2 that stays on the grain surface
  silicate_e_s=110.0d0 ! Energy of the saddle point between a physisorbed and a chemisorbed site (K)
  silicate_e_h2=320.0d0 ! Desorption energy of H2 molecules (K)
  silicate_e_hp=450.0d0 ! Desorption energy of physisorbed H atoms (K)
  silicate_e_hc=3.0d4 ! Desorption energy of chemisorbed H atoms (K)
  silicate_nu_h2=3.0d12 ! Vibrational frequency of H2 molecules in surface sites (s^-1)
  silicate_nu_hc=1.3d13 ! Vibrational frequency of H atoms in their surface sites (s^-1)

  factor1=silicate_mu*flux/(2*silicate_nu_h2*exp(-silicate_e_h2/grain_temperature))

  factor2=1.0d0*(1.0d0+sqrt((silicate_e_hc-silicate_e_s)/(silicate_e_hp-silicate_e_s)))**2 &
      & /4.0d0*exp(-silicate_e_s/grain_temperature)

  epsilon=1.0d0/(1.0d0+silicate_nu_hc/(2*flux)*exp(-1.5*silicate_e_hc/grain_temperature) &
      & *(1.0d0+sqrt((silicate_e_hc-silicate_e_s)/(silicate_e_hp-silicate_e_s)))**2)

  silicate_formation_efficiency=1.0d0/(1.0d0+factor1+factor2)*epsilon

  !  Graphite grain properties
  graphite_mu=0.005d0 ! Fraction of newly formed H2 that stays on the grain surface
  graphite_e_s=260.0d0 ! Energy of the saddle point between a physisorbed and a chemisorbed site (K)
  graphite_e_h2=520.0d0 ! Desorption energy of H2 molecules (K)
  graphite_e_hp=800.0d0 ! Desorption energy of physisorbed H atoms (K)
  graphite_e_hc=3.0d4 ! Desorption energy of chemisorbed H atoms (K)
  graphite_nu_h2=3.0d12 ! Vibrational frequency of H2 molecules in surface sites (s^-1)
  graphite_nu_hc=1.3d13 ! Vibrational frequency of H atoms in their surface sites (s^-1)

  factor1=graphite_mu*flux/(2*graphite_nu_h2*exp(-graphite_e_h2/grain_temperature))

  factor2=1.0d0*(1.0d0+sqrt((graphite_e_hc-graphite_e_s)/(graphite_e_hp-graphite_e_s)))**2 &
      & /4.0d0*exp(-graphite_e_s/grain_temperature)

  epsilon=1.0d0/(1.0d0+graphite_nu_hc/(2*flux)*exp(-1.5*graphite_e_hc/grain_temperature) &
      & *(1.0d0+sqrt((graphite_e_hc-graphite_e_s)/(graphite_e_hp-graphite_e_s)))**2)

  graphite_formation_efficiency=1.0d0/(1.0d0+factor1+factor2)*epsilon

  !!$!  Use the tradional rate, with a simple temperature dependence based on the
  !!$!  thermal velocity of the H atoms in the gas and neglecting any temperature
  !!$!  dependency of the formation and sticking efficiencies
  !!$   RATE=3.0D-18*SQRT(GAS_TEMPERATURE)

  !!$!  Use the treatment of de Jong (1977, A&A, 55, 137, p140 right column).
  !!$!  The second exponential dependence on the gas temperature reduces the
  !!$!  efficiency at high temperatures and so prevents runaway H2 formation
  !!$!  heating at high temperatures:
  !!$!
  !!$!  k_H2 = 3E-18 * T^0.5 * exp(-T/1000)   [cm3/s]
  !!$!
  !!$   RATE=3.0D-18*SQRT(GAS_TEMPERATURE)*EXP(-(GAS_TEMPERATURE/1.0D3))

  !!$!  Use the treatment of Tielens & Hollenbach (1985, ApJ, 291, 722, eqn 4)
  !!$   RATE=0.5D0*THERMAL_VELOCITY*TOTAL_CROSS_SECTION*STICKING_COEFFICIENT

  !  Use the treatment of Cazaux & Tielens (2002, ApJ, 575, L29) and
  !  Cazaux & Tielens (2004, ApJ, 604, 222)
  rate=0.5d0*thermal_velocity*(silicate_cross_section*silicate_formation_efficiency &
      & + graphite_cross_section*graphite_formation_efficiency)*sticking_coefficient*metallicity*100./g2d

  !RATE = RATE / 3.0

  !!$!  Use the expression given by Markus Rollig during the February 2012 Leiden workshop
  !!$   RATE=0.5D0*THERMAL_VELOCITY &
  !!$      & *(SILICATE_CROSS_SECTION/((1.0D0 + 6.998D24/EXP(1.5*SILICATE_E_HC/GRAIN_TEMPERATURE)) &
  !!$      & *(1.0D0 + 1.0D0/(EXP(SILICATE_E_HP/GRAIN_TEMPERATURE) &
  !!$      & *(0.427D0/EXP((SILICATE_E_HP-SILICATE_E_S)/GRAIN_TEMPERATURE) + 2.5336D-14*SQRT(GRAIN_TEMPERATURE))))) &
  !!$      & + GRAPHITE_CROSS_SECTION/((1.0D0 + 4.610D24/EXP(1.5*GRAPHITE_E_HC/GRAIN_TEMPERATURE)) &
  !!$      & *(1.0D0 + 1.0D0/(EXP(GRAPHITE_E_HP/GRAIN_TEMPERATURE) &
  !!$      & *(0.539D0/EXP((GRAPHITE_E_HP-GRAPHITE_E_S)/GRAIN_TEMPERATURE) + 5.6334D-14*SQRT(GRAIN_TEMPERATURE)))))) &
  !!$      & *STICKING_COEFFICIENT

  return
end function h2_formation_rate
!=======================================================================
