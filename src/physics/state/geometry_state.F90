module geometry_state_module
  use definitions, only : dp
  use healpix_types, only : i4b
  implicit none

  type :: columndens_node
    real(kind=dp), pointer :: columndens_point(:,:) => null()
  end type columndens_node

  type :: geometry_state
    real(kind=dp), pointer :: ray_vectors(:,:) => null()
    real(kind=dp), pointer :: radial_order(:) => null()
    integer(kind=i4b), pointer :: point_order(:) => null()
    type(columndens_node), pointer :: column_density(:) => null()
  end type geometry_state

contains

  subroutine allocate_ray_geometry(geometry, nrays)
    type(geometry_state), intent(inout) :: geometry
    integer(kind=i4b), intent(in) :: nrays

    allocate(geometry%ray_vectors(1:3,0:nrays-1))
  end subroutine allocate_ray_geometry

  subroutine allocate_spatial_order(geometry, pdr_count)
    type(geometry_state), intent(inout) :: geometry
    integer(kind=i4b), intent(in) :: pdr_count

    allocate(geometry%radial_order(0:pdr_count))
    allocate(geometry%point_order(1:pdr_count))
  end subroutine allocate_spatial_order

  subroutine allocate_column_density_storage(geometry, pdr_count)
    type(geometry_state), intent(inout) :: geometry
    integer(kind=i4b), intent(in) :: pdr_count

    allocate(geometry%column_density(0:pdr_count))
  end subroutine allocate_column_density_storage

end module geometry_state_module
