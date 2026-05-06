module reaction_rate_kernels_module

  use definitions
  use healpix_types
  use global_module, only : au, grain_radius, kb, pi

  implicit none

contains

  function arrhenius_reaction_rate(alpha, beta, gamma, temperature, minimum_temperature) result(rate)
    real(kind=dp), intent(in) :: alpha, beta, gamma, temperature, minimum_temperature
    real(kind=dp) :: rate

    if (gamma.lt.-200.0d0 .and. temperature.lt.minimum_temperature) then
      rate=0.0d0
    else
      rate=alpha*(temperature/300.0d0)**beta*exp(-(gamma/temperature))
    end if
  end function arrhenius_reaction_rate

  function freeze_out_reaction_rate(alpha, beta, gamma, temperature) result(rate)
    real(kind=dp), intent(in) :: alpha, beta, gamma, temperature
    real(kind=dp) :: rate
    real(kind=dp) :: ion_factor, sticking

    if (beta.eq.0.0d0) then
      ion_factor=1.0d0
    else if (beta.eq.1.0d0) then
      ion_factor=1.0d0+16.71d-4/(grain_radius*temperature)
    else
      ion_factor=0.0d0
    end if

    sticking=0.3d0
    rate=alpha*4.57d4*2.4d-22*sqrt(temperature/gamma)*ion_factor*sticking
  end function freeze_out_reaction_rate

  function cosmic_ray_desorption_rate(gamma, cosmic_ray_ionization_rate) result(rate)
    real(kind=dp), intent(in) :: gamma, cosmic_ray_ionization_rate
    real(kind=dp) :: rate
    real(kind=dp) :: flux, yield

    if (gamma.le.1210.0d0) then
      yield=1.0d5
    else
      yield=0.0d0
    end if

    flux=2.06d-3
    rate=flux*cosmic_ray_ionization_rate*2.4d-22*yield
  end function cosmic_ray_desorption_rate

  function photodesorption_yield(temperature) result(yield)
    real(kind=dp), intent(in) :: temperature
    real(kind=dp) :: yield

    if (temperature.lt.50.0d0) then
      yield=3.5d-3
    else if (temperature.lt.85.0d0) then
      yield=4.0d-3
    else if (temperature.lt.100.0d0) then
      yield=5.5d-3
    else
      yield=7.5d-3
    end if
  end function photodesorption_yield

  function thermal_desorption_rate(alpha, gamma, dust_temperature) result(rate)
    real(kind=dp), intent(in) :: alpha, gamma, dust_temperature
    real(kind=dp) :: rate

    rate=sqrt(2.0d0*1.5d15*kb/(pi**2*au)*alpha/gamma)*exp(-(alpha/dust_temperature))
  end function thermal_desorption_rate

  subroutine clamp_reaction_rate(rate, reaction_id, first_reactant)
    real(kind=dp), intent(inout) :: rate
    integer(kind=i4b), intent(in) :: reaction_id
    character(len=*), intent(in) :: first_reactant

    if (rate.lt.0.0d0) then
      print *,'ERROR! Negative RATE for reaction',reaction_id
      stop
    end if
    if (rate.gt.1.0d0 .and. first_reactant(1:1).ne."G") then
      write(10,*)'WARNING! RATE is too large for reaction',reaction_id
      write(10,*)'RATE =',rate
      rate=1.0d0
    end if
    if (rate.lt.1.0d-99) rate=0.0d0
  end subroutine clamp_reaction_rate

end module reaction_rate_kernels_module
