module species_indices_module
  use healpix_types, only : i4b
  implicit none

  type :: species_index_map
    integer(kind=i4b), pointer :: nh => null()
    integer(kind=i4b), pointer :: nd => null()
    integer(kind=i4b), pointer :: nh2 => null()
    integer(kind=i4b), pointer :: nhd => null()
    integer(kind=i4b), pointer :: nc => null()
    integer(kind=i4b), pointer :: ncx => null()
    integer(kind=i4b), pointer :: nco => null()
    integer(kind=i4b), pointer :: no => null()
    integer(kind=i4b), pointer :: nproton => null()
    integer(kind=i4b), pointer :: nh2o => null()
    integer(kind=i4b), pointer :: nhe => null()
    integer(kind=i4b), pointer :: nmg => null()
    integer(kind=i4b), pointer :: nmgx => null()
    integer(kind=i4b), pointer :: nn => null()
    integer(kind=i4b), pointer :: nfe => null()
    integer(kind=i4b), pointer :: nfex => null()
    integer(kind=i4b), pointer :: nsi => null()
    integer(kind=i4b), pointer :: nsix => null()
    integer(kind=i4b), pointer :: nca => null()
    integer(kind=i4b), pointer :: ncax => null()
    integer(kind=i4b), pointer :: ncaxx => null()
    integer(kind=i4b), pointer :: ns => null()
    integer(kind=i4b), pointer :: nsx => null()
    integer(kind=i4b), pointer :: ncs => null()
    integer(kind=i4b), pointer :: nosh => null()
    integer(kind=i4b), pointer :: ncl => null()
    integer(kind=i4b), pointer :: nclx => null()
    integer(kind=i4b), pointer :: nh2x => null()
    integer(kind=i4b), pointer :: nhex => null()
    integer(kind=i4b), pointer :: nox => null()
    integer(kind=i4b), pointer :: nnx => null()
    integer(kind=i4b), pointer :: nna => null()
    integer(kind=i4b), pointer :: nnax => null()
    integer(kind=i4b), pointer :: nch => null()
    integer(kind=i4b), pointer :: nch2 => null()
    integer(kind=i4b), pointer :: noh => null()
    integer(kind=i4b), pointer :: no2 => null()
    integer(kind=i4b), pointer :: nh3x => null()
    integer(kind=i4b), pointer :: nh3ox => null()
    integer(kind=i4b), pointer :: nhcox => null()
  end type species_index_map

contains

  subroutine allocate_species_index_map(indices)
    type(species_index_map), intent(inout) :: indices

    if (associated(indices%nh)) return
    allocate(indices%nh, indices%nd, indices%nh2, indices%nhd)
    allocate(indices%nc, indices%ncx, indices%nco, indices%no, indices%nproton)
    allocate(indices%nh2o, indices%nhe, indices%nmg, indices%nmgx, indices%nn)
    allocate(indices%nfe, indices%nfex, indices%nsi, indices%nsix)
    allocate(indices%nca, indices%ncax, indices%ncaxx)
    allocate(indices%ns, indices%nsx, indices%ncs, indices%nosh)
    allocate(indices%ncl, indices%nclx, indices%nh2x, indices%nhex)
    allocate(indices%nox, indices%nnx, indices%nna, indices%nnax)
    allocate(indices%nch, indices%nch2, indices%noh, indices%no2)
    allocate(indices%nh3x, indices%nh3ox, indices%nhcox)
  end subroutine allocate_species_index_map

end module species_indices_module
