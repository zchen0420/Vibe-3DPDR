program test_heating_rate_layout

  use definitions
  use healpix_types
  use heating_rate_kernels_module

  implicit none

  real(kind=dp) :: heating_rate(1:12)
  real(kind=dp) :: expected(1:12)

  expected=(/1.0D0,2.0D0,3.0D0,4.0D0,5.0D0,6.0D0,7.0D0,8.0D0,9.0D0,10.0D0,11.0D0,12.0D0/)

  call store_heating_rates(heating_rate, expected(1), expected(2), expected(3), expected(4), &
      & expected(5), expected(6), expected(7), expected(8), expected(9), expected(10), &
      & expected(11), expected(12))

  call assert_equal('dust photoelectric index', heating_dust_photoelectric_index, 1)
  call assert_equal('pah photoelectric index', heating_pah_photoelectric_index, 2)
  call assert_equal('weingartner index', heating_weingartner_index, 3)
  call assert_equal('carbon ionization index', heating_carbon_ionization_index, 4)
  call assert_equal('h2 formation index', heating_h2_formation_index, 5)
  call assert_equal('h2 photodissociation index', heating_h2_photodissociation_index, 6)
  call assert_equal('fuv pumping index', heating_fuv_pumping_index, 7)
  call assert_equal('cosmic ray index', heating_cosmic_ray_index, 8)
  call assert_equal('turbulent index', heating_turbulent_index, 9)
  call assert_equal('chemical index', heating_chemical_index, 10)
  call assert_equal('gas grain index', heating_gas_grain_index, 11)
  call assert_equal('total index', heating_total_index, 12)
  call assert_vector_equal('heating rate layout', heating_rate, expected)

  write(6,*) 'test_heating_rate_layout: ok'

contains

  subroutine assert_equal(label, actual, expected_value)
    character(len=*), intent(in) :: label
    integer(kind=i4b), intent(in) :: actual, expected_value

    if (actual.ne.expected_value) then
      write(6,*) 'test_heating_rate_layout failed: ', trim(label), actual, expected_value
      stop 1
    endif
  end subroutine assert_equal

  subroutine assert_vector_equal(label, actual, expected_value)
    character(len=*), intent(in) :: label
    real(kind=dp), intent(in) :: actual(1:12), expected_value(1:12)

    if (any(abs(actual-expected_value).gt.1.0D-12)) then
      write(6,*) 'test_heating_rate_layout failed: ', trim(label)
      stop 1
    endif
  end subroutine assert_vector_equal

end program test_heating_rate_layout
