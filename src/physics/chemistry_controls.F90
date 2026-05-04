module chemistry_module
  USE ISO_C_BINDING
  use healpix_types

  real(kind=dp),bind(c,name='chemistry_module_mp_relative_abundance_tolerance_') :: relative_abundance_tolerance
  real(kind=dp),bind(c,name='chemistry_module_mp_absolute_abundance_tolerance_') :: absolute_abundance_tolerance
end module chemistry_module
