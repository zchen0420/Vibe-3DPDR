module simulation_grid_module
  use definitions, only : dp
  use healpix_types, only : i4b
  use coolants_module, only : point_coolant_state
  implicit none

  type :: pdr_node
    integer(kind=i4b), pointer :: epray(:) => null()
    integer(kind=i4b), pointer :: projected(:,:) => null()
    integer(kind=i4b), pointer :: raytype(:) => null()
    real(kind=dp), pointer :: epoint(:,:,:) => null()
    real(kind=dp), pointer :: columndensity(:) => null()
    real(kind=dp), pointer :: AV(:) => null()
    real(kind=dp), pointer :: rad_surface(:) => null()
    real(kind=dp), pointer :: abundance(:) => null()
    type(point_coolant_state), allocatable :: coolant_state(:)
    real(kind=dp) :: UVfield
    real(kind=dp) :: rho
    real(kind=dp) :: position(1:3)
    integer(kind=i4b) :: etype
#ifdef DUST2
    real(kind=dp) :: dust_temperature
#endif
  end type pdr_node

  type :: simulation_grid
    type(pdr_node), pointer :: points(:) => null()
    integer(kind=i4b), pointer :: pdr_ids(:) => null()
    integer(kind=i4b), pointer :: ion_ids(:) => null()
    integer(kind=i4b), pointer :: dark_ids(:) => null()
  end type simulation_grid

contains

  subroutine allocate_simulation_grid(grid, point_count)
    type(simulation_grid), intent(inout) :: grid
    integer(kind=i4b), intent(in) :: point_count

    allocate(grid%points(1:point_count))
    allocate(grid%pdr_ids(1:point_count))
    allocate(grid%ion_ids(1:point_count))
    allocate(grid%dark_ids(1:point_count))
  end subroutine allocate_simulation_grid

end module simulation_grid_module
