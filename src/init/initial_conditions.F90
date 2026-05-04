module initial_conditions_module
  use healpix_types
  use maincode_module
  use global_module, only : metallicity

contains

  subroutine initialise_temperatures
    dusttemperature = dust_temperature

#ifndef GUESS_TEMP
    gastemperature = Tguess
    previousgastemperature = Tguess
#ifdef THERMALBALANCE
    Tlow = Tlow0
    Thigh = Thigh0
#endif
#endif
  end subroutine initialise_temperatures

  subroutine initialise_particle_abundances
    integer(kind=i4b) :: species_index

    do p=1,grand_ptot
      allocate(pdr(p)%abundance(1:nspec))
      do species_index=1,nspec
        if (species(species_index).eq.'H2'.OR.species(species_index).eq.'H'.OR.species(species_index).eq.'He') then
          pdr(p)%abundance(species_index) = dummyabundance(species_index)
        else
          pdr(p)%abundance(species_index) = dummyabundance(species_index)*metallicity
        endif
      enddo
    enddo
  end subroutine initialise_particle_abundances

end module initial_conditions_module
