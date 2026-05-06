module heating_rate_kernels_module

  use definitions
  use healpix_types
  use global_module, only : grain_radius, g2d, metallicity, pi

  implicit none

  integer(kind=i4b), parameter :: heating_dust_photoelectric_index = 1
  integer(kind=i4b), parameter :: heating_pah_photoelectric_index = 2
  integer(kind=i4b), parameter :: heating_weingartner_index = 3
  integer(kind=i4b), parameter :: heating_carbon_ionization_index = 4
  integer(kind=i4b), parameter :: heating_h2_formation_index = 5
  integer(kind=i4b), parameter :: heating_h2_photodissociation_index = 6
  integer(kind=i4b), parameter :: heating_fuv_pumping_index = 7
  integer(kind=i4b), parameter :: heating_cosmic_ray_index = 8
  integer(kind=i4b), parameter :: heating_turbulent_index = 9
  integer(kind=i4b), parameter :: heating_chemical_index = 10
  integer(kind=i4b), parameter :: heating_gas_grain_index = 11
  integer(kind=i4b), parameter :: heating_total_index = 12

contains

  function cosmic_ray_heating_rate(cosmic_ray_ionization_rate, density, h2_abundance) result(heating_rate)
    real(kind=dp), intent(in) :: cosmic_ray_ionization_rate, density, h2_abundance
    real(kind=dp) :: heating_rate

    heating_rate=(9.4D0*ev)*(1.3D-17*cosmic_ray_ionization_rate)*density*h2_abundance
  end function cosmic_ray_heating_rate

  function turbulent_heating_rate(turbulent_velocity, density) result(heating_rate)
    real(kind=dp), intent(in) :: turbulent_velocity, density
    real(kind=dp) :: heating_rate
    real(kind=dp) :: turbulent_scale

    turbulent_scale=5.0D0
    heating_rate=3.5D-28*((turbulent_velocity/1.0D5)**3)*(1.0D0/turbulent_scale)*density
  end function turbulent_heating_rate

  function gas_grain_exchange_rate(density, gas_temperature, dust_temperature) result(heating_rate)
    real(kind=dp), intent(in) :: density, gas_temperature, dust_temperature
    real(kind=dp) :: heating_rate
    real(kind=dp) :: accommodation, grain_density, grain_cross_section

    accommodation=0.35D0*exp(-sqrt((dust_temperature+gas_temperature)/5.0D2))+0.1D0
    grain_density=1.998D-12*density*metallicity*100.0D0/g2d
    grain_cross_section=pi*grain_radius**2
    heating_rate=4.003D-12*density*grain_density*grain_cross_section*accommodation*sqrt(gas_temperature) &
        & *(dust_temperature-gas_temperature)
  end function gas_grain_exchange_rate

  subroutine store_heating_rates(heating_rate, dust_photoelectric, pah_photoelectric, weingartner, &
      & carbon_ionization, h2_formation, h2_photodissociation, fuv_pumping, cosmic_ray, turbulent, &
      & chemical, gas_grain, total)
    real(kind=dp), intent(out) :: heating_rate(1:12)
    real(kind=dp), intent(in) :: dust_photoelectric, pah_photoelectric, weingartner
    real(kind=dp), intent(in) :: carbon_ionization, h2_formation, h2_photodissociation
    real(kind=dp), intent(in) :: fuv_pumping, cosmic_ray, turbulent, chemical, gas_grain, total

    heating_rate(heating_dust_photoelectric_index)=dust_photoelectric
    heating_rate(heating_pah_photoelectric_index)=pah_photoelectric
    heating_rate(heating_weingartner_index)=weingartner
    heating_rate(heating_carbon_ionization_index)=carbon_ionization
    heating_rate(heating_h2_formation_index)=h2_formation
    heating_rate(heating_h2_photodissociation_index)=h2_photodissociation
    heating_rate(heating_fuv_pumping_index)=fuv_pumping
    heating_rate(heating_cosmic_ray_index)=cosmic_ray
    heating_rate(heating_turbulent_index)=turbulent
    heating_rate(heating_chemical_index)=chemical
    heating_rate(heating_gas_grain_index)=gas_grain
    heating_rate(heating_total_index)=total
  end subroutine store_heating_rates

end module heating_rate_kernels_module
