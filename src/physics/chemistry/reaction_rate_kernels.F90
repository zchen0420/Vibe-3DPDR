module reaction_rate_kernels_module

  use definitions
  use healpix_types
  use global_module, only : au, grain_radius, kb, pi

  implicit none

contains

  function arrhenius_reaction_rate(alpha, beta, gamma, temperature, minimum_temperature) result(rate)
    real(kind=dp), intent(in) :: alpha, beta, gamma, temperature, minimum_temperature
    real(kind=dp) :: rate

    if (gamma.lt.-200.0D0 .and. temperature.lt.minimum_temperature) then
      rate=0.0D0
    else
      rate=alpha*(temperature/300.0D0)**beta*exp(-(gamma/temperature))
    endif
  end function arrhenius_reaction_rate

  function freeze_out_reaction_rate(alpha, beta, gamma, temperature) result(rate)
    real(kind=dp), intent(in) :: alpha, beta, gamma, temperature
    real(kind=dp) :: rate
    real(kind=dp) :: ion_factor, sticking

    if (beta.eq.0.0D0) then
      ion_factor=1.0D0
    else if (beta.eq.1.0D0) then
      ion_factor=1.0D0+16.71D-4/(grain_radius*temperature)
    else
      ion_factor=0.0D0
    endif

    sticking=0.3D0
    rate=alpha*4.57D4*2.4D-22*sqrt(temperature/gamma)*ion_factor*sticking
  end function freeze_out_reaction_rate

  function cosmic_ray_desorption_rate(gamma, cosmic_ray_ionization_rate) result(rate)
    real(kind=dp), intent(in) :: gamma, cosmic_ray_ionization_rate
    real(kind=dp) :: rate
    real(kind=dp) :: flux, yield

    if (gamma.le.1210.0D0) then
      yield=1.0D5
    else
      yield=0.0D0
    endif

    flux=2.06D-3
    rate=flux*cosmic_ray_ionization_rate*2.4D-22*yield
  end function cosmic_ray_desorption_rate

  function photodesorption_yield(temperature) result(yield)
    real(kind=dp), intent(in) :: temperature
    real(kind=dp) :: yield

    if (temperature.lt.50.0D0) then
      yield=3.5D-3
    else if (temperature.lt.85.0D0) then
      yield=4.0D-3
    else if (temperature.lt.100.0D0) then
      yield=5.5D-3
    else
      yield=7.5D-3
    endif
  end function photodesorption_yield

  function thermal_desorption_rate(alpha, gamma, dust_temperature) result(rate)
    real(kind=dp), intent(in) :: alpha, gamma, dust_temperature
    real(kind=dp) :: rate

    rate=sqrt(2.0D0*1.5D15*kb/(pi**2*au)*alpha/gamma)*exp(-(alpha/dust_temperature))
  end function thermal_desorption_rate

  subroutine clamp_reaction_rate(rate, reaction_id, first_reactant)
    real(kind=dp), intent(inout) :: rate
    integer(kind=i4b), intent(in) :: reaction_id
    character(len=*), intent(in) :: first_reactant

    if (rate.lt.0.0D0) then
      print *,'ERROR! Negative RATE for reaction',reaction_id
      stop
    endif
    if (rate.gt.1.0D0 .and. first_reactant(1:1).ne."G") then
      write(10,*)'WARNING! RATE is too large for reaction',reaction_id
      write(10,*)'RATE =',rate
      rate=1.0D0
    endif
    if (rate.lt.1.0D-99) rate=0.0D0
  end subroutine clamp_reaction_rate

end module reaction_rate_kernels_module
