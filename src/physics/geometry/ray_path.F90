module ray_path_module
  use definitions, only : dp
  use healpix_types, only : i4b
  use simulation_grid_module, only : pdr_node, simulation_grid
  implicit none

contains

  real(kind=dp) function distance_between_points(left_position, right_position)
    real(kind=dp), intent(in) :: left_position(1:3)
    real(kind=dp), intent(in) :: right_position(1:3)

    distance_between_points = sqrt(sum((left_position(1:3)-right_position(1:3))**2))
  end function distance_between_points

  real(kind=dp) function ray_step_length(point, ray_index, eval_index)
    type(pdr_node), intent(in) :: point
    integer(kind=i4b), intent(in) :: ray_index
    integer(kind=i4b), intent(in) :: eval_index

    ray_step_length = distance_between_points(point%epoint(1:3,ray_index,eval_index-1), &
        & point%epoint(1:3,ray_index,eval_index))
  end function ray_step_length

  real(kind=dp) function ray_distance_from_origin(point, ray_index, eval_index)
    type(pdr_node), intent(in) :: point
    integer(kind=i4b), intent(in) :: ray_index
    integer(kind=i4b), intent(in) :: eval_index

    ray_distance_from_origin = distance_between_points(point%epoint(1:3,ray_index,0), &
        & point%epoint(1:3,ray_index,eval_index))
  end function ray_distance_from_origin

  real(kind=dp) function adaptive_step_length(point, ray_index, eval_index)
    type(pdr_node), intent(in) :: point
    integer(kind=i4b), intent(in) :: ray_index
    integer(kind=i4b), intent(in) :: eval_index

    adaptive_step_length = ray_distance_from_origin(point, ray_index, eval_index) - &
        & ray_distance_from_origin(point, ray_index, eval_index-1)
  end function adaptive_step_length

  integer(kind=i4b) function projected_point_id(point, ray_index, eval_index)
    type(pdr_node), intent(in) :: point
    integer(kind=i4b), intent(in) :: ray_index
    integer(kind=i4b), intent(in) :: eval_index

    projected_point_id = point%projected(ray_index,eval_index)
  end function projected_point_id

  subroutine set_ray_origin_projection(point, point_id)
    type(pdr_node), intent(inout) :: point
    integer(kind=i4b), intent(in) :: point_id

    point%projected(:,0) = point_id
  end subroutine set_ray_origin_projection

  subroutine initialize_ray_origin(point, point_id, source_origin)
    type(pdr_node), intent(inout) :: point
    integer(kind=i4b), intent(in) :: point_id
    real(kind=dp), intent(in) :: source_origin(1:3)
    integer(kind=i4b) :: ray_index

    point%epray = 0
    call set_ray_origin_projection(point, point_id)
    do ray_index=lbound(point%epoint,2),ubound(point%epoint,2)
      point%epoint(1:3,ray_index,0) = source_origin(1:3)
    enddo
  end subroutine initialize_ray_origin

  subroutine store_projected_point(grid, source_point_id, source_origin, projected_position, ray_index, &
      & target_point_id, max_eval_points, kill_current_ray)
    type(simulation_grid), intent(inout) :: grid
    integer(kind=i4b), intent(in) :: source_point_id
    real(kind=dp), intent(in) :: source_origin(1:3)
    real(kind=dp), intent(in) :: projected_position(1:3)
    integer(kind=i4b), intent(in) :: ray_index
    integer(kind=i4b), intent(in) :: target_point_id
    integer(kind=i4b), intent(in) :: max_eval_points
    logical, intent(inout) :: kill_current_ray
    integer(kind=i4b) :: projected_eval_index

    grid%points(source_point_id)%epray(ray_index) = grid%points(source_point_id)%epray(ray_index)+1
    projected_eval_index = grid%points(source_point_id)%epray(ray_index)
    if (projected_eval_index.gt.max_eval_points) stop 'Increase maxpoints!'

    grid%points(source_point_id)%epoint(1:3,ray_index,projected_eval_index) = &
        & projected_position(1:3)+source_origin(1:3)
    grid%points(source_point_id)%projected(ray_index,projected_eval_index) = target_point_id
    if (grid%points(target_point_id)%etype.eq.2) kill_current_ray=.true.
  end subroutine store_projected_point

  subroutine update_minimum_adaptive_step(point, adaptive_minimum)
    type(pdr_node), intent(in) :: point
    real(kind=dp), intent(inout) :: adaptive_minimum
    integer(kind=i4b) :: eval_index
    integer(kind=i4b) :: ray_index
    real(kind=dp) :: step_length

    do ray_index=lbound(point%epray,1),ubound(point%epray,1)
      if (point%epray(ray_index).le.0) cycle
      do eval_index=1,point%epray(ray_index)
        step_length = adaptive_step_length(point, ray_index, eval_index)
        if (step_length.lt.0.0_dp) stop 'found negative adaptive step!'
        if (step_length.lt.adaptive_minimum) adaptive_minimum = step_length
      enddo
    enddo
  end subroutine update_minimum_adaptive_step

end module ray_path_module
