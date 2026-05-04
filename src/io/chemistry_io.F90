module chemistry_io_module
  use maincode_module

contains

  subroutine load_chemistry_network
    call read_species(nspec, species, dummyabundance, mass)
    call READ_RATES(NREAC,REACTANT,PRODUCT,ALPHA,BETA,GAMMA,rate,DUPLICATE,RTMIN,RTMAX)
  end subroutine load_chemistry_network

end module chemistry_io_module
