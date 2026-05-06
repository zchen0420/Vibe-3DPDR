module chemistry_io_module
  use maincode_module
  use global_module, only : initialize_species_indices

contains

  subroutine load_chemistry_network
    call initialize_species_indices
    call read_species(nspec, chemistry%network%species, chemistry%network%initial_abundance, chemistry%network%mass)
    call read_rates(nreac,chemistry%network%reactant,chemistry%network%product,chemistry%network%alpha,chemistry%network%beta,chemistry%network%gamma,chemistry%rate,chemistry%network%duplicate,chemistry%network%rtmin,chemistry%network%rtmax)
  end subroutine load_chemistry_network

end module chemistry_io_module
