module radiation_module
  use definitions, only : dp
  use healpix_types, only : i4b
  use ray_path_module, only : projected_point_id, ray_step_length, set_ray_origin_projection
  implicit none

contains

  integer(kind=i4b) function incident_ray_index(field_vector, ray_vectors, ray_count)
  integer(kind=i4b), intent(in) :: ray_count
  real(kind=dp), intent(in) :: field_vector(:)
  real(kind=dp), intent(in) :: ray_vectors(1:3,0:ray_count-1)

  integer(kind=i4b) :: ray_index
  real(kind=dp) :: best_projection
  real(kind=dp) :: ray_projection

  incident_ray_index = 0
  best_projection = 0.0d0

  do ray_index=0,ray_count-1
    ray_projection = -ray_dot_product(field_vector, ray_vectors(:,ray_index))
    if (ray_projection.ge.best_projection) then
      best_projection = ray_projection
      incident_ray_index = ray_index
    end if
  end do
end function incident_ray_index

real(kind=dp) function ray_dot_product(field_vector, ray_vector)
real(kind=dp), intent(in) :: field_vector(:)
real(kind=dp), intent(in) :: ray_vector(1:3)
integer(kind=i4b) :: component_index

ray_dot_product = 0.0d0
do component_index=1,min(size(field_vector),3)
  ray_dot_product = ray_dot_product + field_vector(component_index)*ray_vector(component_index)
end do
end function ray_dot_product

real(kind=dp) function attenuated_uv(ray_surface, visual_extinction, uv_factor)
real(kind=dp), intent(in) :: ray_surface
real(kind=dp), intent(in) :: visual_extinction
real(kind=dp), intent(in) :: uv_factor

attenuated_uv = ray_surface*exp(-visual_extinction*uv_factor)
if (attenuated_uv.lt.1.0d-50) attenuated_uv = 0.0d0
end function attenuated_uv

subroutine calc_uvfield
  use maincode_module
#ifdef THERMALBALANCE
  use maincode_module, only : tmin_array, tmax_array
#endif

  integer(kind=i4b) :: grid_id
  integer(kind=i4b) :: pdr_index
  integer(kind=i4b) :: ray_index
  integer(kind=i4b) :: eval_index
  integer(kind=i4b) :: selected_ray_index
  real(kind=dp) :: ray_step
  real(kind=dp) :: guessed_temperature

  do grid_id=1,grand_ptot
    grid%points(grid_id)%uvfield = 0.0d0
  end do

  do pdr_index=1,pdr_ptot
    grid_id=grid%pdr_ids(pdr_index)
    call calculate_point_uvfield(grid_id)
  end do

#ifdef THERMALBALANCE
#ifdef GUESS_TEMP
  allocate(tmin_array(0:pdr_ptot))
  allocate(tmax_array(0:pdr_ptot))
  do pdr_index=1,pdr_ptot
    guessed_temperature = 10.0d0*(1.0d0+(2.*grid%points(grid%pdr_ids(pdr_index))%uvfield)**(1.0d0/3.0d0))
    thermal%gas_temperature(pdr_index) = guessed_temperature
    thermal%previous_gas_temperature(pdr_index) = guessed_temperature

    thermal%low_temperature(pdr_index) = guessed_temperature/2.0d0
    thermal%high_temperature(pdr_index) = guessed_temperature*1.5d0

    if (thermal%low_temperature(pdr_index).lt.tmin)  thermal%low_temperature(pdr_index)  = tmin
    if (thermal%high_temperature(pdr_index).gt.tmax) thermal%high_temperature(pdr_index) = tmax

    tmin_array(pdr_index) = thermal%low_temperature(pdr_index)/3.0d0
    tmax_array(pdr_index) = thermal%high_temperature(pdr_index)*2.0d0
    if (tmin_array(pdr_index).lt.tmin) tmin_array(pdr_index) = tmin
    if (tmax_array(pdr_index).gt.tmax) tmax_array(pdr_index) = tmax
  end do
#endif
#endif

  if (dark_ptot.gt.0) then
    grid_id=grid%dark_ids(1)
    call calculate_point_uvfield(grid_id)
    thermal%gas_temperature(0)=tmin
    thermal%previous_gas_temperature(0)=tmin
#ifdef THERMALBALANCE
    thermal%low_temperature(0)=tmin
    thermal%high_temperature(0)=tmin
#endif
  end if

contains

  subroutine calculate_point_uvfield(point_id)
    integer(kind=i4b), intent(in) :: point_id

    grid%points(point_id)%rad_surface = 0.0d0
    grid%points(point_id)%columndensity = 0.0d0
    grid%points(point_id)%av = 0.0d0
    call set_ray_origin_projection(grid%points(point_id), point_id)

    do ray_index=0,nrays-1
      if (grid%points(point_id)%epray(ray_index).gt.0) then
        do eval_index=1,grid%points(point_id)%epray(ray_index)
          ray_step = ray_step_length(grid%points(point_id), ray_index, eval_index)
          grid%points(point_id)%columndensity(ray_index) = grid%points(point_id)%columndensity(ray_index) + &
              &((grid%points(projected_point_id(grid%points(point_id), ray_index, eval_index-1))%rho + &
              &grid%points(projected_point_id(grid%points(point_id), ray_index, eval_index))%rho)/2.)*ray_step*pc
        end do
      end if
      grid%points(point_id)%av(ray_index) = grid%points(point_id)%columndensity(ray_index)*runtime%av_scale
    end do

    if (fieldchoice.eq.'UNI') selected_ray_index = incident_ray_index(gext(:), geometry%ray_vectors, nrays)

    do ray_index=0,nrays-1
      if (fieldchoice.eq.'ISO') then
        grid%points(point_id)%rad_surface(ray_index) = gext(1)/real(nrays,kind=dp)
      else if (fieldchoice.eq.'UNI') then
        if (ray_index.eq.selected_ray_index) then
          grid%points(point_id)%rad_surface(ray_index) = -ray_dot_product(gext(:),geometry%ray_vectors(:,ray_index))
        end if
      end if
    end do

    do ray_index=0,nrays-1
      grid%points(point_id)%uvfield = grid%points(point_id)%uvfield + &
          &attenuated_uv(grid%points(point_id)%rad_surface(ray_index), grid%points(point_id)%av(ray_index), runtime%uv_scale)
    end do
  end subroutine calculate_point_uvfield

end subroutine calc_uvfield

end module radiation_module
