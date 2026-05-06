module thermal_balance_module
  use definitions, only : dp
  use healpix_types, only : i4b
  use heating_rates_module, only : calculate_heating_rates, heating_rate_environment
  use maincode_module, only : chemistry, Fcrit, first_time, grid, level_conv, pdr_ptot, &
      &runtime, Tdiff, thermal, Tmin, Tmax
  use global_module, only : all_heating
  use point_reaction_rates_module, only : calculate_point_reaction_rates
  use reaction_rates_module, only : reaction_rate_indices

  implicit none

contains

  subroutine calculate_heating_and_temperature_updates(dobinarychop, previouschange)
    logical, allocatable, intent(inout) :: dobinarychop(:)
    character(len=1), allocatable, intent(inout) :: previouschange(:)

    integer(kind=i4b) :: point_index
    integer(kind=i4b) :: point_id
    type(heating_rate_environment) :: heating_environment
    type(reaction_rate_indices) :: rate_indices
    real(kind=dp) :: point_heating(1:12)
    real(kind=dp) :: current_temperature

    do point_index=1,pdr_ptot
      point_id=grid%pdr_ids(point_index)
#ifdef THERMALBALANCE
      if (thermal%thermal_converged(point_index)) cycle
#endif

      call calculate_point_reaction_rates(point_index, point_id, rate_indices)
      call set_heating_environment(point_index, point_id, heating_environment)
      call calculate_heating_rates(heating_environment, rate_indices, point_heating)

      all_heating(point_index,:)=point_heating

#ifdef THERMALBALANCE
      thermal%mean_balance(point_index) = all_heating(point_index,12) - thermal%total_cooling_rate(point_index)
      thermal%relative_balance(point_index) = 2.0D0*abs(thermal%mean_balance(point_index))/ &
          &abs(all_heating(point_index,12) + thermal%total_cooling_rate(point_index))

      if (level_conv.and.first_time) then
        current_temperature = thermal%gas_temperature(point_index)

        if (thermal%mean_balance(point_index).eq.0) then
          thermal%low_temperature(point_index) = thermal%gas_temperature(point_index)
          thermal%high_temperature(point_index) = thermal%gas_temperature(point_index)
        else if (thermal%mean_balance(point_index).gt.0) then
          thermal%low_temperature(point_index) = thermal%gas_temperature(point_index)
          thermal%gas_temperature(point_index) = 1.3D0*thermal%gas_temperature(point_index)
          previouschange(point_index) = "H"
        else if (thermal%mean_balance(point_index).lt.0) then
          thermal%high_temperature(point_index) = thermal%gas_temperature(point_index)
          thermal%gas_temperature(point_index) = 0.7D0*thermal%gas_temperature(point_index)
          if (thermal%gas_temperature(point_index).lt.Tmin) thermal%gas_temperature(point_index)=Tmin
          previouschange(point_index) = "C"
        endif

        thermal%previous_gas_temperature(point_index) = current_temperature
        if (thermal%previous_gas_temperature(point_index).lt.Tmin) thermal%previous_gas_temperature(point_index)=Tmin

      else if (level_conv.and..not.first_time) then
        current_temperature = thermal%gas_temperature(point_index)

        if (thermal%relative_balance(point_index).le.Fcrit) thermal%thermal_converged(point_index) = .true.

        if (.not.dobinarychop(point_index)) then
          if (thermal%mean_balance(point_index).gt.0.and.previouschange(point_index).eq."H") then
            thermal%low_temperature(point_index) = thermal%gas_temperature(point_index)
            thermal%gas_temperature(point_index) = 1.3D0*thermal%gas_temperature(point_index)
            thermal%high_temperature(point_index) = thermal%gas_temperature(point_index)
            previouschange(point_index) = "H"
          endif

          if (thermal%mean_balance(point_index).lt.0.and.previouschange(point_index).eq."C") then
            thermal%high_temperature(point_index) = thermal%gas_temperature(point_index)
            thermal%gas_temperature(point_index) = 0.7D0*thermal%gas_temperature(point_index)
            thermal%low_temperature(point_index) = thermal%gas_temperature(point_index)
            previouschange(point_index) = "C"
            if (thermal%gas_temperature(point_index).lt.Tmin) then
              thermal%gas_temperature(point_index)=Tmin
              thermal%low_temperature(point_index)=Tmin
              thermal%high_temperature(point_index)=Tmin
            endif
          endif

          if (thermal%mean_balance(point_index).gt.0.and.previouschange(point_index).eq."C") then
            thermal%gas_temperature(point_index) = (thermal%high_temperature(point_index) + thermal%low_temperature(point_index))/2.0D0
            dobinarychop(point_index)=.true.
          endif

          if (thermal%mean_balance(point_index).lt.0.and.previouschange(point_index).eq."H") then
            thermal%gas_temperature(point_index) = (thermal%high_temperature(point_index) + thermal%low_temperature(point_index))/2.0D0
            dobinarychop(point_index)=.true.
          endif

        else
          if (thermal%mean_balance(point_index).gt.0) then
            thermal%low_temperature(point_index) = thermal%gas_temperature(point_index)
            thermal%gas_temperature(point_index) = (thermal%gas_temperature(point_index) + thermal%high_temperature(point_index)) / 2.0D0
          endif
          if (thermal%mean_balance(point_index).lt.0) then
            thermal%high_temperature(point_index) = thermal%gas_temperature(point_index)
            thermal%gas_temperature(point_index) = (thermal%gas_temperature(point_index) + thermal%low_temperature(point_index)) / 2.0D0
          endif
        endif

#ifdef TEMP_FIX
        if ((abs(thermal%gas_temperature(point_index)-thermal%previous_gas_temperature(point_index)).le.Tdiff).and. &
            &(thermal%relative_balance(point_index).gt.Fcrit)) then
          thermal%thermal_converged(point_index)=.true.
        endif
#endif

        thermal%previous_gas_temperature(point_index) = current_temperature

        if ((current_temperature.lt.Tmin).and.(thermal%mean_balance(point_index).lt.0)) thermal%thermal_converged(point_index)=.true.
        if ((current_temperature.gt.Tmax).and.(thermal%mean_balance(point_index).gt.0)) thermal%thermal_converged(point_index)=.true.

        if (thermal%thermal_converged(point_index)) then
          if (current_temperature.lt.Tmin) then
            thermal%previous_gas_temperature(point_index) = Tmin
            thermal%gas_temperature(point_index) = Tmin
            if (thermal%force_level_minimum(point_index)) then
              thermal%thermal_converged(point_index)=.true.
            else
              thermal%thermal_converged(point_index)=.false.
              thermal%level_population_converged(point_index)=.false.
              current_temperature=Tmin
            endif
            thermal%force_level_minimum(point_index)=.true.
          endif

          if (current_temperature.gt.Tmax) then
            thermal%previous_gas_temperature(point_index) = Tmax
            thermal%gas_temperature(point_index) = Tmax
          endif
        endif
      endif
#endif
    enddo
  end subroutine calculate_heating_and_temperature_updates

  subroutine set_heating_environment(point_index, point_id, environment)
    integer(kind=i4b), intent(in) :: point_index
    integer(kind=i4b), intent(in) :: point_id
    type(heating_rate_environment), intent(out) :: environment

    environment%density = grid%points(point_id)%rho
    environment%gas_temperature = thermal%gas_temperature(point_index)
    environment%dust_temperature = thermal%dust_temperature(point_index)
    environment%uv_field = grid%points(point_id)%UVfield
    environment%turbulent_velocity = runtime%turbulent_velocity
    environment%abundance => grid%points(point_id)%abundance
    environment%reaction_rate => chemistry%rate
  end subroutine set_heating_environment

end module thermal_balance_module
