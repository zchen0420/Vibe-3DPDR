module species_indices_module
  use healpix_types, only : i4b
  implicit none

  type :: species_index_map
    integer(kind=i4b), pointer :: NH => null()
    integer(kind=i4b), pointer :: ND => null()
    integer(kind=i4b), pointer :: NH2 => null()
    integer(kind=i4b), pointer :: NHD => null()
    integer(kind=i4b), pointer :: NC => null()
    integer(kind=i4b), pointer :: NCx => null()
    integer(kind=i4b), pointer :: NCO => null()
    integer(kind=i4b), pointer :: NO => null()
    integer(kind=i4b), pointer :: NPROTON => null()
    integer(kind=i4b), pointer :: NH2O => null()
    integer(kind=i4b), pointer :: NHe => null()
    integer(kind=i4b), pointer :: NMG => null()
    integer(kind=i4b), pointer :: NMGx => null()
    integer(kind=i4b), pointer :: NN => null()
    integer(kind=i4b), pointer :: NFE => null()
    integer(kind=i4b), pointer :: NFEx => null()
    integer(kind=i4b), pointer :: NSI => null()
    integer(kind=i4b), pointer :: NSIx => null()
    integer(kind=i4b), pointer :: NCA => null()
    integer(kind=i4b), pointer :: NCAx => null()
    integer(kind=i4b), pointer :: NCAxx => null()
    integer(kind=i4b), pointer :: NS => null()
    integer(kind=i4b), pointer :: NSx => null()
    integer(kind=i4b), pointer :: NCS => null()
    integer(kind=i4b), pointer :: NOSH => null()
    integer(kind=i4b), pointer :: NCL => null()
    integer(kind=i4b), pointer :: NCLx => null()
    integer(kind=i4b), pointer :: NH2x => null()
    integer(kind=i4b), pointer :: NHEx => null()
    integer(kind=i4b), pointer :: NOx => null()
    integer(kind=i4b), pointer :: NNx => null()
    integer(kind=i4b), pointer :: NNA => null()
    integer(kind=i4b), pointer :: NNAx => null()
    integer(kind=i4b), pointer :: NCH => null()
    integer(kind=i4b), pointer :: NCH2 => null()
    integer(kind=i4b), pointer :: NOH => null()
    integer(kind=i4b), pointer :: NO2 => null()
    integer(kind=i4b), pointer :: NH3x => null()
    integer(kind=i4b), pointer :: NH3Ox => null()
    integer(kind=i4b), pointer :: NHCOx => null()
  end type species_index_map

contains

  subroutine allocate_species_index_map(indices)
    type(species_index_map), intent(inout) :: indices

    if (associated(indices%NH)) return
    allocate(indices%NH, indices%ND, indices%NH2, indices%NHD)
    allocate(indices%NC, indices%NCx, indices%NCO, indices%NO, indices%NPROTON)
    allocate(indices%NH2O, indices%NHe, indices%NMG, indices%NMGx, indices%NN)
    allocate(indices%NFE, indices%NFEx, indices%NSI, indices%NSIx)
    allocate(indices%NCA, indices%NCAx, indices%NCAxx)
    allocate(indices%NS, indices%NSx, indices%NCS, indices%NOSH)
    allocate(indices%NCL, indices%NCLx, indices%NH2x, indices%NHEx)
    allocate(indices%NOx, indices%NNx, indices%NNA, indices%NNAx)
    allocate(indices%NCH, indices%NCH2, indices%NOH, indices%NO2)
    allocate(indices%NH3x, indices%NH3Ox, indices%NHCOx)
  end subroutine allocate_species_index_map

end module species_indices_module
