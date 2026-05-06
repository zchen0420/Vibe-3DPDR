module columns_module
  use definitions, only : dp
  use healpix_types, only : i4b
  use ray_path_module, only : projected_point_id, ray_step_length, set_ray_origin_projection
  implicit none

contains

  real(kind=dp) function column_increment(left_density, right_density, left_abundance, right_abundance, step_length_pc)
  real(kind=dp), intent(in) :: left_density
  real(kind=dp), intent(in) :: right_density
  real(kind=dp), intent(in) :: left_abundance
  real(kind=dp), intent(in) :: right_abundance
  real(kind=dp), intent(in) :: step_length_pc

  column_increment = step_length_pc*(left_density*left_abundance + right_density*right_abundance)/2.0d0
end function column_increment

subroutine calc_columndens(allocate_storage)
  use maincode_module

  logical, intent(in) :: allocate_storage
  integer(kind=i4b) :: point_id
  integer(kind=i4b) :: pdr_index
  integer(kind=i4b) :: ray_index
  integer(kind=i4b) :: eval_index
  integer(kind=i4b) :: species_index
  real(kind=dp) :: ray_step

  if (dark_ptot.gt.0) then
    point_id=grid%dark_ids(1)
    if (allocate_storage) allocate(geometry%column_density(0)%columndens_point(0:nrays-1,1:nspec))
    geometry%column_density(0)%columndens_point = 0.0d0
    call calculate_point_columns(point_id, 0)
  end if

  do pdr_index=1,pdr_ptot
    point_id=grid%pdr_ids(pdr_index)
    call set_ray_origin_projection(grid%points(point_id), point_id)
#ifdef THERMALBALANCE
    if (thermal%thermal_converged(pdr_index)) cycle
#endif
    if (allocate_storage) allocate(geometry%column_density(pdr_index)%columndens_point(0:nrays-1,1:nspec))
    geometry%column_density(pdr_index)%columndens_point = 0.0d0
    call calculate_point_columns(point_id, pdr_index)
  end do

contains

  subroutine calculate_point_columns(point_id, column_index)
    integer(kind=i4b), intent(in) :: point_id
    integer(kind=i4b), intent(in) :: column_index
    integer(kind=i4b) :: left_point_id
    integer(kind=i4b) :: right_point_id

    do ray_index=0,nrays-1
      if (grid%points(point_id)%epray(ray_index).eq.0) cycle
      do eval_index=1,grid%points(point_id)%epray(ray_index)
        ray_step = ray_step_length(grid%points(point_id), ray_index, eval_index)
        left_point_id = projected_point_id(grid%points(point_id), ray_index, eval_index-1)
        right_point_id = projected_point_id(grid%points(point_id), ray_index, eval_index)

        do species_index=1,nspec
          geometry%column_density(column_index)%columndens_point(ray_index,species_index) = &
              &geometry%column_density(column_index)%columndens_point(ray_index,species_index) + &
              &column_increment(grid%points(left_point_id)%rho, grid%points(right_point_id)%rho, &
              &grid%points(left_point_id)%abundance(species_index), &
              &grid%points(right_point_id)%abundance(species_index), ray_step*pc)
        end do
      end do
    end do
  end subroutine calculate_point_columns

end subroutine calc_columndens

end module columns_module
