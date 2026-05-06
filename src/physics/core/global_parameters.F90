module global_module

  USE ISO_C_BINDING
  use definitions
  use healpix_types
  use species_indices_module, only : species_index_map, allocate_species_index_map

  type(species_index_map) :: species_idx
  integer(kind=i4b),bind(c,name='global_module_mp_nelect_')::NELECT

  ! REAL(kind=dp), save :: ZETA=3.85D0,OMEGA=0.42D0,GRAIN_RADIUS=1.0D-5,METALLICITY=1.0D0
  REAL(kind=dp) :: g2d
  REAL(kind=dp) :: metallicity
  REAL(kind=dp) :: omega
  REAL(kind=dp) :: grain_radius
  ! REAL(kind=dp), save :: OMEGA=0.42D0,GRAIN_RADIUS=1.0D-7,METALLICITY=1.0D0
  ! REAL(kind=dp), save :: ZETA=1.0D0,OMEGA=0.42D0,GRAIN_RADIUS=1.0D-5,METALLICITY=1.0D0

  real(kind=dp),allocatable :: all_heating(:,:)

contains

  subroutine initialize_species_indices
    call allocate_species_index_map(species_idx)
  end subroutine initialize_species_indices

end module global_module
