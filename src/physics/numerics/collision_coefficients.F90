module collision_coefficients_module
  use coolants_module, only : collision_partner_count, coolant_data
  use definitions, only : dp
  use healpix_types, only : i4b
  implicit none

contains

  subroutine calculate_collision_coefficients(coolant_table, gas_temperature, collider_density, coefficients)
    type(coolant_data), intent(in) :: coolant_table
    real(kind=dp), intent(in) :: gas_temperature
    real(kind=dp), intent(in) :: collider_density(1:collision_partner_count)
    real(kind=dp), intent(out) :: coefficients(:,:)

    integer(kind=i4b) :: lower_temperature_index
    integer(kind=i4b) :: partner_id
    integer(kind=i4b) :: upper_temperature_index
    real(kind=dp) :: interpolation_fraction

    call assert_collision_dimensions(coolant_table, coefficients)

    coefficients = 0.0d0
    do partner_id=1,collision_partner_count
      if (coolant_table%collision_temperatures(partner_id,1).eq.0.0d0) cycle

      call find_temperature_bracket(coolant_table%collision_temperatures(partner_id,:), gas_temperature, &
          &lower_temperature_index, upper_temperature_index, interpolation_fraction)
      call add_partner_collision_rates(coolant_table, partner_id, lower_temperature_index, &
          &upper_temperature_index, interpolation_fraction, collider_density(partner_id), coefficients)
    end do
  end subroutine calculate_collision_coefficients

  subroutine assert_collision_dimensions(coolant_table, coefficients)
    type(coolant_data), intent(in) :: coolant_table
    real(kind=dp), intent(in) :: coefficients(:,:)

    if (size(coefficients,1).ne.coolant_table%nlevels .or. size(coefficients,2).ne.coolant_table%nlevels) then
      stop 'calculate_collision_coefficients received inconsistent dimensions'
    end if
  end subroutine assert_collision_dimensions

  subroutine find_temperature_bracket(temperatures, gas_temperature, lower_index, upper_index, interpolation_fraction)
    real(kind=dp), intent(in) :: temperatures(:)
    real(kind=dp), intent(in) :: gas_temperature
    integer(kind=i4b), intent(out) :: lower_index
    integer(kind=i4b), intent(out) :: upper_index
    real(kind=dp), intent(out) :: interpolation_fraction

    integer(kind=i4b) :: temperature_index
    integer(kind=i4b) :: temperature_count

    temperature_count = size(temperatures)
    lower_index = 0
    upper_index = 0

    do temperature_index=1,temperature_count
      if (temperatures(temperature_index).gt.gas_temperature) then
        lower_index = temperature_index - 1
        upper_index = temperature_index
        exit
      else if (temperatures(temperature_index).eq.0.0d0) then
        lower_index = temperature_index - 1
        upper_index = temperature_index - 1
        exit
      end if
    end do

    if (upper_index.eq.0) then
      lower_index = temperature_count
      upper_index = temperature_count
    else if (upper_index.eq.1) then
      lower_index = 1
      upper_index = 1
    end if

    if (lower_index.eq.upper_index) then
      interpolation_fraction = 0.0d0
    else
      interpolation_fraction = (gas_temperature-temperatures(lower_index)) &
          &/(temperatures(upper_index)-temperatures(lower_index))
    end if
  end subroutine find_temperature_bracket

  subroutine add_partner_collision_rates(coolant_table, partner_id, lower_temperature_index, upper_temperature_index, &
        &interpolation_fraction, density, coefficients)
    type(coolant_data), intent(in) :: coolant_table
    integer(kind=i4b), intent(in) :: partner_id
    integer(kind=i4b), intent(in) :: lower_temperature_index
    integer(kind=i4b), intent(in) :: upper_temperature_index
    real(kind=dp), intent(in) :: interpolation_fraction
    real(kind=dp), intent(in) :: density
    real(kind=dp), intent(inout) :: coefficients(:,:)

    coefficients = coefficients + density*(coolant_table%collision_rates(partner_id,:,:,lower_temperature_index) &
        &+(coolant_table%collision_rates(partner_id,:,:,upper_temperature_index) &
        &-coolant_table%collision_rates(partner_id,:,:,lower_temperature_index))*interpolation_fraction)
  end subroutine add_partner_collision_rates

end module collision_coefficients_module
