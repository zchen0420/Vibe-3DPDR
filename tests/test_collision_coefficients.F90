program test_collision_coefficients
  use collision_coefficients_module, only : calculate_collision_coefficients
  use coolants_module, only : COLLIDER_H, COLLIDER_H2, COLLISION_PARTNER_COUNT, coolant_data
  use definitions, only : dp
  implicit none

  type(coolant_data) :: coolant_table
  real(kind=dp) :: collider_density(1:COLLISION_PARTNER_COUNT)
  real(kind=dp) :: coefficients(2,2)

  coolant_table%nlevels = 2
  coolant_table%ntemperatures = 2
  allocate(coolant_table%collision_temperatures(1:COLLISION_PARTNER_COUNT,1:2))
  allocate(coolant_table%collision_rates(1:COLLISION_PARTNER_COUNT,1:2,1:2,1:2))

  coolant_table%collision_temperatures = 0.0D0
  coolant_table%collision_rates = 0.0D0
  coolant_table%collision_temperatures(COLLIDER_H,1:2) = (/100.0D0, 200.0D0/)
  coolant_table%collision_temperatures(COLLIDER_H2,1:2) = (/100.0D0, 200.0D0/)
  coolant_table%collision_rates(COLLIDER_H,2,1,1:2) = (/1.0D0, 3.0D0/)
  coolant_table%collision_rates(COLLIDER_H2,2,1,1:2) = (/2.0D0, 4.0D0/)

  collider_density = 0.0D0
  collider_density(COLLIDER_H) = 2.0D0
  collider_density(COLLIDER_H2) = 3.0D0

  call calculate_collision_coefficients(coolant_table, 150.0D0, collider_density, coefficients)

  call assert_close(coefficients(2,1), 13.0D0, 1.0D-12, 'interpolated partner sum')
  call assert_close(coefficients(1,2), 0.0D0, 1.0D-12, 'unset transition')

  write(6,*) 'test_collision_coefficients: ok'

contains

  subroutine assert_close(actual, expected, tolerance, label)
    real(kind=dp), intent(in) :: actual
    real(kind=dp), intent(in) :: expected
    real(kind=dp), intent(in) :: tolerance
    character(len=*), intent(in) :: label

    if (abs(actual - expected) > tolerance) then
      write(6,*) 'test_collision_coefficients failed: ', trim(label), actual, expected
      stop 1
    endif
  end subroutine assert_close

end program test_collision_coefficients
