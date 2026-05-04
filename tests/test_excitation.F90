program test_excitation
  use healpix_types
  use excitation_module
  implicit none

  integer, parameter :: nlev = 2
  real(kind=dp) :: energies(nlev)
  real(kind=dp) :: weights(nlev)
  real(kind=dp) :: partition
  real(kind=dp) :: populations(nlev)
  real(kind=dp) :: density
  real(kind=dp) :: temperature

  energies = 0.0D0
  weights = 1.0D0
  density = 10.0D0
  temperature = 100.0D0

  call calculate_partition_function(partition, nlev, energies, weights, temperature)
  call assert_close(partition, 2.0D0, 1.0D-12, 'partition')

  call calculate_lte_populations(nlev, populations, energies, weights, partition, density, temperature)
  call assert_close(populations(1), 5.0D0, 1.0D-12, 'population(1)')
  call assert_close(populations(2), 5.0D0, 1.0D-12, 'population(2)')
  call assert_close(sum(populations), density, 1.0D-12, 'population sum')

  write(6,*) 'test_excitation: ok'

contains

  subroutine assert_close(actual, expected, tolerance, label)
    real(kind=dp), intent(in) :: actual
    real(kind=dp), intent(in) :: expected
    real(kind=dp), intent(in) :: tolerance
    character(len=*), intent(in) :: label

    if (abs(actual - expected) > tolerance) then
      write(6,*) 'test_excitation failed: ', trim(label), actual, expected
      stop 1
    endif
  end subroutine assert_close

end program test_excitation
