module initial_conditions_module
  use healpix_types
  use maincode_module
  use global_module, only : metallicity

contains

  subroutine initialise_temperatures
    thermal%dust_temperature = dust_temperature

#ifndef GUESS_TEMP
    thermal%gas_temperature = runtime%temperature_guess
    thermal%previous_gas_temperature = runtime%temperature_guess
#ifdef THERMALBALANCE
    thermal%low_temperature = tlow0
    thermal%high_temperature = thigh0
#endif
#endif
  end subroutine initialise_temperatures

  subroutine initialise_particle_abundances
    integer(kind=i4b) :: point_id
    integer(kind=i4b) :: species_index

    do point_id=1,grand_ptot
      allocate(grid%points(point_id)%abundance(1:nspec))
      do species_index=1,nspec
        if (chemistry%network%species(species_index).eq.'H2'.or.chemistry%network%species(species_index).eq.'H'.or.chemistry%network%species(species_index).eq.'He') then
          grid%points(point_id)%abundance(species_index) = chemistry%network%initial_abundance(species_index)
        else
          grid%points(point_id)%abundance(species_index) = chemistry%network%initial_abundance(species_index)*metallicity
        end if
      end do
    end do
  end subroutine initialise_particle_abundances

end module initial_conditions_module
