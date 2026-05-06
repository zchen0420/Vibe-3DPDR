module spatial_index_module
  use definitions
  use healpix_types
  use geometry_state_module, only : allocate_spatial_order
  use maincode_module

contains

  subroutine prepare_spatial_index_inputs
    call swap_last_two_pdr_points
    call build_pdr_point_index
  end subroutine prepare_spatial_index_inputs

  subroutine swap_last_two_pdr_points
    integer(kind=i4b) :: point_id, point_index
    real(kind=dp), allocatable :: position_rev(:,:), density_rev(:)

    if (pdr_ptot.lt.2) return

    allocate(position_rev(1:3,1:pdr_ptot))
    allocate(density_rev(1:pdr_ptot))

    do point_index=1,pdr_ptot-2
      point_id=grid%pdr_ids(point_index)
      position_rev(1:3,point_index)=grid%points(point_id)%position
      density_rev(point_index)=grid%points(point_id)%rho
    end do

    position_rev(1:3,pdr_ptot-1)=grid%points(grid%pdr_ids(pdr_ptot))%position
    density_rev(pdr_ptot-1)=grid%points(grid%pdr_ids(pdr_ptot))%rho

    position_rev(1:3,pdr_ptot)=grid%points(grid%pdr_ids(pdr_ptot-1))%position
    density_rev(pdr_ptot)=grid%points(grid%pdr_ids(pdr_ptot-1))%rho

    do point_index=1,pdr_ptot
      point_id=grid%pdr_ids(point_index)
      grid%points(point_id)%position=position_rev(1:3,point_index)
      grid%points(point_id)%rho=density_rev(point_index)
    end do

    deallocate(position_rev)
    deallocate(density_rev)
  end subroutine swap_last_two_pdr_points

  subroutine build_pdr_point_index
    integer(kind=i4b) :: point_id, point_index

    call allocate_spatial_order(geometry, pdr_ptot)

    do point_index=1,pdr_ptot
      point_id=grid%pdr_ids(point_index)
      geometry%radial_order(point_index) = sqrt(sum(grid%points(point_id)%position**2))
      geometry%point_order(point_index) = point_id
    end do
  end subroutine build_pdr_point_index

end module spatial_index_module
