subroutine evaluation_points
  !calculation of evaluation points
  !T.Bisbas
  use definitions
  use healpix_types
  use healpix_module
  use maincode_module
  use ray_path_module, only : distance_between_points, initialize_ray_origin, projected_point_id, &
      & store_projected_point, update_minimum_adaptive_step

  double precision::adaptivemin
  integer(kind=i4b) :: pdr_index, point_index, point_id
  integer(kind=i4b) :: total_evaluation_points

  if (dark_ptot.gt.0) then
    write(6,*) 'Creating evaluation points for the Dark Molecular element'
    call build_evaluation_points_for_source(grid%dark_ids(1))
  endif

  write(6,*) 'Proceeding for the grid%points (SERIAL)...'
  do pdr_index=1,pdr_ptot
    call build_evaluation_points_for_source(grid%pdr_ids(pdr_index))
  enddo

  total_evaluation_points=0
  do point_index=1,pdr_ptot
    point_id=grid%pdr_ids(point_index)
    total_evaluation_points = total_evaluation_points + sum(grid%points(point_id)%epray(:))
  enddo
  if (dark_ptot.gt.0) then
    total_evaluation_points = total_evaluation_points + sum(grid%points(grid%dark_ids(1))%epray(:))
  endif
  write(6,*) 'No. evaluation points:',total_evaluation_points
  write(6,*) 'Done!';write(6,*) ''

  write(6,*) 'Checking for negative steps...'
  adaptivemin=100.0D0
  do point_index=1,pdr_ptot
    call update_minimum_adaptive_step(grid%points(grid%pdr_ids(point_index)), adaptivemin)
  enddo
  if (dark_ptot.gt.0) then
    call update_minimum_adaptive_step(grid%points(grid%dark_ids(1)), adaptivemin)
  endif

  write(6,*) 'No negative steps found'
  write(6,*) 'Minimum adaptive step = ',adaptivemin

  write(6,*) 'Assigning raytypes'
  do point_index=1,pdr_ptot
    call assign_raytypes(grid%pdr_ids(point_index))
  enddo

  if (dark_ptot.gt.0) then
    call assign_raytypes(grid%dark_ids(1))
  endif

  return

contains

  subroutine build_evaluation_points_for_source(source_point_id)
    integer(kind=i4b), intent(in) :: source_point_id
    logical :: killray(0:nrays-1)
    integer(kind=i4b) :: candidate_point_id
    integer(kind=i4b) :: pixel_index
    integer(kind=i4b) :: sorted_point_count
    integer(kind=i4b) :: sorted_point_index
    integer(kind=i4b) :: sorted_point_total
    integer(kind=i4b), allocatable :: sorted_point_ids(:)
    real(kind=dp) :: evaluation_point(1:3,0:nrays-1)
    real(kind=dp) :: line_of_sight_angle
    real(kind=dp) :: max_origin_distance
    real(kind=dp) :: relative_vector(1:3)
    real(kind=dp) :: source_origin(1:3)
    real(kind=dp), allocatable :: sorted_distance(:)

    killray=.false.
    source_origin(1:3) = grid%points(source_point_id)%position
    allocate(sorted_distance(0:grand_ptot-1)) !needs one extra place for sorting in heapsort
    allocate(sorted_point_ids(1:grand_ptot-1)) !-1 to avoid overlapping source_origin & current point

    call initialize_ray_origin(grid%points(source_point_id), source_point_id, source_origin)

    sorted_point_count=0
    do candidate_point_id=1,grand_ptot
      if (candidate_point_id.eq.source_point_id) cycle
      sorted_point_count=sorted_point_count+1
      sorted_distance(sorted_point_count)=distance_between_points(grid%points(candidate_point_id)%position, source_origin)
      sorted_point_ids(sorted_point_count)=candidate_point_id
    enddo

    sorted_point_total=sorted_point_count
    call validate_sorted_point_count(sorted_point_total)
    call heapsort(sorted_point_total,sorted_point_ids,sorted_distance)

    max_origin_distance=sorted_distance(sorted_point_total)
    evaluation_point=0.0D0

    do sorted_point_index=1,sorted_point_total
      relative_vector=grid%points(sorted_point_ids(sorted_point_index))%position-source_origin
      call ray_index_for_relative_vector(relative_vector, pixel_index)
      if (killray(pixel_index)) cycle

      call project_candidate_on_ray(relative_vector,evaluation_point(1:3,pixel_index), &
          & max_origin_distance,pixel_index,line_of_sight_angle)
      if (line_of_sight_angle.le.theta_crit) then
        call store_projected_point(grid, source_point_id, source_origin, evaluation_point(1:3,pixel_index), &
            & pixel_index, sorted_point_ids(sorted_point_index), maxpoints, killray(pixel_index))
      endif
    enddo

    deallocate(sorted_distance)
    deallocate(sorted_point_ids)
  end subroutine build_evaluation_points_for_source

  subroutine validate_sorted_point_count(sorted_point_total)
    integer(kind=i4b), intent(in) :: sorted_point_total

    if (sorted_point_total.ne.(grand_ptot-1)) then
      write(6,*) 'sorted_point_total = ',sorted_point_total,' grand_ptot-1 = ',grand_ptot-1
      stop 'sorted_point_total is not equal to grand_ptot-1 !!'
    endif
  end subroutine validate_sorted_point_count

  subroutine ray_index_for_relative_vector(relative_vector, pixel_index)
    real(kind=dp), intent(in) :: relative_vector(1:3)
    integer(kind=i4b), intent(out) :: pixel_index
    real(kind=dp) :: ray_phi
    real(kind=dp) :: ray_theta

    call vec2ang(relative_vector,ray_theta,ray_phi)
    call ang2pix_nest_id(nside,ray_theta,ray_phi,pixel_index)
  end subroutine ray_index_for_relative_vector

  subroutine assign_raytypes(point_id)
    integer(kind=i4b), intent(in) :: point_id
    integer(kind=i4b) :: ray_index
    integer(kind=i4b) :: terminal_point_id

    if (.not.associated(grid%points(point_id)%raytype)) allocate(grid%points(point_id)%raytype(0:nrays-1))
    do ray_index=0,nrays-1
      if (grid%points(point_id)%epray(ray_index).gt.0) then
        terminal_point_id = projected_point_id(grid%points(point_id), ray_index, grid%points(point_id)%epray(ray_index))
        grid%points(point_id)%raytype(ray_index) = -grid%points(terminal_point_id)%etype
      else
        grid%points(point_id)%raytype(ray_index) = -grid%points(point_id)%etype
      endif
    enddo
  end subroutine assign_raytypes

  subroutine project_candidate_on_ray(relative_vector,evaluation_point_on_ray,max_origin_distance, &
      & pixel_index,line_of_sight_angle)
    real(kind=dp), intent(in) :: relative_vector(1:3)
    real(kind=dp), intent(inout) :: evaluation_point_on_ray(1:3)
    real(kind=dp), intent(in) :: max_origin_distance
    integer(kind=i4b), intent(in) :: pixel_index
    real(kind=dp), intent(out) :: line_of_sight_angle
    real(kind=dp) :: ray_endpoint(1:3)

    ray_endpoint(1:3) = 1.1_DP*max_origin_distance*geometry%ray_vectors(1:3,pixel_index)
    line_of_sight_angle=acos(dot_product(relative_vector(1:3)-evaluation_point_on_ray(1:3), &
        & ray_endpoint(1:3)-evaluation_point_on_ray(1:3))/ &
        & (distance_between_points(relative_vector(1:3), evaluation_point_on_ray(1:3)) * &
        & distance_between_points(ray_endpoint(1:3), evaluation_point_on_ray(1:3))))

    if (ray_endpoint(3).ne.0.0_dp) then
      evaluation_point_on_ray(3) = (ray_endpoint(1)*ray_endpoint(3)*relative_vector(1) + &
          & ray_endpoint(2)*ray_endpoint(3)*relative_vector(2) + &
          & (ray_endpoint(3)**2)*relative_vector(3))/(ray_endpoint(1)**2+ray_endpoint(2)**2+ &
          & ray_endpoint(3)**2)
      evaluation_point_on_ray(2) = evaluation_point_on_ray(3)*ray_endpoint(2)/ray_endpoint(3)
      evaluation_point_on_ray(1) = evaluation_point_on_ray(3)*ray_endpoint(1)/ray_endpoint(3)
    else
      if (ray_endpoint(1).eq.0.0_dp) then
        evaluation_point_on_ray(1) = 0.0_dp
        evaluation_point_on_ray(2) = relative_vector(2)
        evaluation_point_on_ray(3) = 0.0_dp
      else if (ray_endpoint(2).eq.0.0_dp) then
        evaluation_point_on_ray(1) = relative_vector(1)
        evaluation_point_on_ray(2) = 0.0_dp
        evaluation_point_on_ray(3) = 0.0_dp
      else
        evaluation_point_on_ray(3) = 0.0_dp
        evaluation_point_on_ray(1) = ((ray_endpoint(1)**2)*relative_vector(1) + &
            & ray_endpoint(1)*ray_endpoint(2)*relative_vector(2))/(ray_endpoint(1)**2 + ray_endpoint(2)**2)
        evaluation_point_on_ray(2) = evaluation_point_on_ray(1) * ray_endpoint(2)/ray_endpoint(1)
      endif
    endif
  end subroutine project_candidate_on_ray

end subroutine evaluation_points
