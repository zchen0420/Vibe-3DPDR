module thermal_balance_module
  use definitions, only : dp
  use healpix_types, only : i4b
  use maincode_module, only : chemistry, Fcrit, first_time, grid, level_conv, nreac, nspec, pdr_ptot, &
      &runtime, Tdiff, thermal, Tmin, Tmax
  use global_module, only : all_heating
  use point_reaction_rates_module, only : calculate_point_reaction_rates, reaction_rate_indices

  implicit none

contains

  subroutine calculate_heating_and_temperature_updates(dobinarychop, previouschange)
    logical, allocatable, intent(inout) :: dobinarychop(:)
    character(len=1), allocatable, intent(inout) :: previouschange(:)

    integer(kind=i4b) :: point_index
    integer(kind=i4b) :: point_id
    type(reaction_rate_indices) :: rate_indices
    real(kind=dp) :: point_heating(1:12)
    real(kind=dp) :: current_temperature

    do point_index=1,pdr_ptot
      point_id=grid%pdr_ids(point_index)
#ifdef THERMALBALANCE
      if (thermal%thermal_converged(point_index)) cycle
#endif

      call calculate_point_reaction_rates(point_index, point_id, rate_indices)
      call calculate_heating_rates(grid%points(point_id)%rho,thermal%gas_temperature(point_index),thermal%dust_temperature(point_index),&
          &grid%points(point_id)%UVfield,runtime%turbulent_velocity,nspec,grid%points(point_id)%abundance(:),nreac,chemistry%rate,point_heating,&
          &rate_indices%grain_surface, rate_indices%h2_photodissociation, &
          &rate_indices%hd_photodissociation, rate_indices%co_photodissociation, &
          &rate_indices%carbon_photoionization, rate_indices%silicon_photoionization)

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

end module thermal_balance_module
