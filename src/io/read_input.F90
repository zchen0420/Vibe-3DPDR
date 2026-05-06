module coolant_input_module
  use coolants_module, only : COLLISION_PARTNER_COUNT, coolant_data
  use healpix_types, only : c, dp, hp, i4b, kb
  implicit none

contains

  subroutine read_lamda_coolant_file(coolant_table)
    type(coolant_data), intent(inout) :: coolant_table

    integer(kind=i4b), parameter :: input_unit = 8
    integer(kind=i4b) :: file_level_count
    integer(kind=i4b) :: level_index
    integer(kind=i4b) :: line_count
    integer(kind=i4b) :: line_id
    integer(kind=i4b) :: line_index
    integer(kind=i4b) :: partner_count
    integer(kind=i4b) :: partner_index
    integer(kind=i4b) :: partner_id
    integer(kind=i4b) :: collision_count
    integer(kind=i4b) :: temperature_count
    integer(kind=i4b) :: i, j
    real(kind=dp) :: einstein_a
    real(kind=dp) :: energy
    real(kind=dp) :: frequency
    real(kind=dp) :: weight

    call clear_coolant_table(coolant_table)

    open(input_unit,file=coolant_table%input_file,status='old')
    read(input_unit,'(////)')
    read(input_unit,*) file_level_count
    if (file_level_count.ne.coolant_table%nlevels) stop 'Incorrect number of energy levels in coolant file'

    read(input_unit,*)
    do level_index=1,file_level_count
      read(input_unit,*) i, energy, weight
      coolant_table%energies(i) = energy*c*hp
      coolant_table%weights(i) = weight
    enddo

    read(input_unit,*)
    read(input_unit,*) line_count
    read(input_unit,*)
    do line_index=1,line_count
      read(input_unit,*) line_id, i, j, einstein_a, frequency
      coolant_table%frequencies(i,j) = frequency*1.0D9
      coolant_table%frequencies(j,i) = coolant_table%frequencies(i,j)
      coolant_table%a_coeffs(i,j) = einstein_a
      coolant_table%b_coeffs(i,j) = coolant_table%a_coeffs(i,j) &
          &/(2.0D0*hp*(coolant_table%frequencies(i,j)**3)/(c**2))
      coolant_table%b_coeffs(j,i) = coolant_table%b_coeffs(i,j) &
          &*(coolant_table%weights(i)/coolant_table%weights(j))
    enddo

    call fill_missing_transition_frequencies(coolant_table)

    read(input_unit,*)
    read(input_unit,*) partner_count
    do partner_index=1,partner_count
      read(input_unit,*)
      read(input_unit,*) partner_id
      call assert_valid_collision_partner(partner_id)

      read(input_unit,*)
      read(input_unit,*) collision_count
      read(input_unit,*)
      read(input_unit,*) temperature_count
      if (temperature_count.gt.coolant_table%ntemperatures) then
        write(6,*) 'ERROR! Too many temperature values (>NTEMP):', temperature_count
        stop
      endif

      call read_collision_partner(input_unit, partner_id, collision_count, temperature_count, coolant_table)
    enddo

    close(input_unit)
    write(6,*) 'Cooling datafile: ',trim(coolant_table%input_file),' read successfully'
  end subroutine read_lamda_coolant_file

  subroutine clear_coolant_table(coolant_table)
    type(coolant_data), intent(inout) :: coolant_table

    coolant_table%energies = 0.0D0
    coolant_table%weights = 0.0D0
    coolant_table%a_coeffs = 0.0D0
    coolant_table%b_coeffs = 0.0D0
    coolant_table%frequencies = 0.0D0
    coolant_table%collision_temperatures = 0.0D0
    coolant_table%collision_rates = 0.0D0
  end subroutine clear_coolant_table

  subroutine fill_missing_transition_frequencies(coolant_table)
    type(coolant_data), intent(inout) :: coolant_table

    integer(kind=i4b) :: i, j
    real(kind=dp) :: calculated_frequency

    do i=1,coolant_table%nlevels
      do j=1,coolant_table%nlevels
        calculated_frequency = abs(coolant_table%energies(i)-coolant_table%energies(j))/hp
        if (coolant_table%frequencies(i,j).ne.0.0D0) then
          if (abs(calculated_frequency-coolant_table%frequencies(i,j))/coolant_table%frequencies(i,j).gt.1.0D-2) then
            write(6,*) 'ERROR! Calculated frequency differs by >1%:'
            write(6,*) calculated_frequency,' Hz vs',coolant_table%frequencies(i,j),' Hz'
            stop
          endif
        else
          coolant_table%frequencies(i,j) = calculated_frequency
        endif
      enddo
    enddo
  end subroutine fill_missing_transition_frequencies

  subroutine assert_valid_collision_partner(partner_id)
    integer(kind=i4b), intent(in) :: partner_id

    if (partner_id.lt.1 .or. partner_id.gt.COLLISION_PARTNER_COUNT) then
      write(6,*) 'ERROR! Unrecognized collision partner ID:', partner_id
      stop
    endif
  end subroutine assert_valid_collision_partner

  subroutine read_collision_partner(input_unit, partner_id, collision_count, temperature_count, coolant_table)
    integer(kind=i4b), intent(in) :: input_unit
    integer(kind=i4b), intent(in) :: partner_id
    integer(kind=i4b), intent(in) :: collision_count
    integer(kind=i4b), intent(in) :: temperature_count
    type(coolant_data), intent(inout) :: coolant_table

    integer(kind=i4b) :: collision_index
    integer(kind=i4b) :: row_id
    integer(kind=i4b) :: lower_level
    integer(kind=i4b) :: upper_level
    integer(kind=i4b) :: temperature_index
    real(kind=dp) :: coefficients(1:coolant_table%ntemperatures)

    read(input_unit,*)
    read(input_unit,*) (coolant_table%collision_temperatures(partner_id,temperature_index), &
        &temperature_index=1,temperature_count)
    read(input_unit,*)

    do collision_index=1,collision_count
      read(input_unit,*) row_id, upper_level, lower_level, &
          &(coefficients(temperature_index),temperature_index=1,temperature_count)

      do temperature_index=1,temperature_count
        coolant_table%collision_rates(partner_id,upper_level,lower_level,temperature_index) = &
            &coefficients(temperature_index)
        call fill_reverse_collision_rate(coolant_table, partner_id, upper_level, lower_level, temperature_index)
      enddo
    enddo
  end subroutine read_collision_partner

  subroutine fill_reverse_collision_rate(coolant_table, partner_id, upper_level, lower_level, temperature_index)
    type(coolant_data), intent(inout) :: coolant_table
    integer(kind=i4b), intent(in) :: partner_id
    integer(kind=i4b), intent(in) :: upper_level
    integer(kind=i4b), intent(in) :: lower_level
    integer(kind=i4b), intent(in) :: temperature_index

    real(kind=dp) :: downward_rate
    real(kind=dp) :: excitation_factor

    downward_rate = coolant_table%collision_rates(partner_id,upper_level,lower_level,temperature_index)
    if (downward_rate.eq.0.0D0) return
    if (coolant_table%collision_rates(partner_id,lower_level,upper_level,temperature_index).ne.0.0D0) return

    excitation_factor = (coolant_table%weights(upper_level)/coolant_table%weights(lower_level)) &
        &*exp(-(coolant_table%energies(upper_level)-coolant_table%energies(lower_level)) &
        &/(kb*coolant_table%collision_temperatures(partner_id,temperature_index)))
    coolant_table%collision_rates(partner_id,lower_level,upper_level,temperature_index) = &
        &downward_rate*excitation_factor
  end subroutine fill_reverse_collision_rate

end module coolant_input_module
