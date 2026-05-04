module radiation_module
  use definitions, only : dp
  use healpix_types, only : i4b
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
    best_projection = 0.0D0

    do ray_index=0,ray_count-1
      ray_projection = -ray_dot_product(field_vector, ray_vectors(:,ray_index))
      if (ray_projection.ge.best_projection) then
        best_projection = ray_projection
        incident_ray_index = ray_index
      endif
    enddo
  end function incident_ray_index

  real(kind=dp) function ray_dot_product(field_vector, ray_vector)
    real(kind=dp), intent(in) :: field_vector(:)
    real(kind=dp), intent(in) :: ray_vector(1:3)
    integer(kind=i4b) :: component_index

    ray_dot_product = 0.0D0
    do component_index=1,min(size(field_vector),3)
      ray_dot_product = ray_dot_product + field_vector(component_index)*ray_vector(component_index)
    enddo
  end function ray_dot_product

  real(kind=dp) function attenuated_uv(ray_surface, visual_extinction, uv_factor)
    real(kind=dp), intent(in) :: ray_surface
    real(kind=dp), intent(in) :: visual_extinction
    real(kind=dp), intent(in) :: uv_factor

    attenuated_uv = ray_surface*exp(-visual_extinction*uv_factor)
    if (attenuated_uv.lt.1.0D-50) attenuated_uv = 0.0D0
  end function attenuated_uv

  subroutine calc_uvfield
    use maincode_module
#ifdef THERMALBALANCE
    use maincode_module, only : Tmin_array, Tmax_array
#endif

    integer(kind=i4b) :: grid_id
    integer(kind=i4b) :: pdr_index
    integer(kind=i4b) :: ray_index
    integer(kind=i4b) :: eval_index
    integer(kind=i4b) :: selected_ray_index
    real(kind=dp) :: ray_step

    do grid_id=1,grand_ptot
      pdr(grid_id)%UVfield = 0.0D0
    enddo

    do pdr_index=1,pdr_ptot
      grid_id=IDlist_pdr(pdr_index)
      call calculate_point_uvfield(grid_id)
    enddo

#ifdef THERMALBALANCE
#ifdef GUESS_TEMP
    allocate(Tmin_array(0:pdr_ptot))
    allocate(Tmax_array(0:pdr_ptot))
    do pdr_index=1,pdr_ptot
      Tguess = 10.0D0*(1.0D0+(2.*pdr(IDlist_pdr(pdr_index))%UVfield)**(1.0D0/3.0D0))
      gastemperature(pdr_index) = Tguess
      previousgastemperature(pdr_index) = Tguess

      Tlow(pdr_index) = Tguess/2.0D0
      Thigh(pdr_index) = Tguess*1.5D0

      if (Tlow(pdr_index).lt.Tmin)  Tlow(pdr_index)  = Tmin
      if (Thigh(pdr_index).gt.Tmax) Thigh(pdr_index) = Tmax

      Tmin_array(pdr_index) = Tlow(pdr_index)/3.0D0
      Tmax_array(pdr_index) = Thigh(pdr_index)*2.0D0
      if (Tmin_array(pdr_index).lt.Tmin) Tmin_array(pdr_index) = Tmin
      if (Tmax_array(pdr_index).gt.Tmax) Tmax_array(pdr_index) = Tmax
    enddo
#endif
#endif

    if (dark_ptot.gt.0) then
      grid_id=IDlist_dark(1)
      call calculate_point_uvfield(grid_id)
      gastemperature(0)=Tmin
      previousgastemperature(0)=Tmin
#ifdef THERMALBALANCE
      Tlow(0)=Tmin
      Thigh(0)=Tmin
#endif
    endif

  contains

    subroutine calculate_point_uvfield(point_id)
      integer(kind=i4b), intent(in) :: point_id

      pdr(point_id)%rad_surface = 0.0D0
      pdr(point_id)%columndensity = 0.0D0
      pdr(point_id)%AV = 0.0D0
      pdr(point_id)%projected(:,0)=point_id

      do ray_index=0,nrays-1
        if (pdr(point_id)%epray(ray_index).gt.0) then
          do eval_index=1,pdr(point_id)%epray(ray_index)
            ray_step = sqrt((pdr(point_id)%epoint(1,ray_index,eval_index-1)-&
                &pdr(point_id)%epoint(1,ray_index,eval_index))**2 + &
                &(pdr(point_id)%epoint(2,ray_index,eval_index-1)-&
                &pdr(point_id)%epoint(2,ray_index,eval_index))**2 + &
                &(pdr(point_id)%epoint(3,ray_index,eval_index-1)-&
                &pdr(point_id)%epoint(3,ray_index,eval_index))**2)
            pdr(point_id)%columndensity(ray_index) = pdr(point_id)%columndensity(ray_index) + &
                &((pdr(int(pdr(point_id)%projected(ray_index,eval_index-1)))%rho + &
                &pdr(int(pdr(point_id)%projected(ray_index,eval_index)))%rho)/2.)*ray_step*pc
          enddo
        endif
        pdr(point_id)%AV(ray_index) = pdr(point_id)%columndensity(ray_index)*AV_fac
      enddo

      if (fieldchoice.eq.'UNI') selected_ray_index = incident_ray_index(Gext(:), vectors, nrays)

      do ray_index=0,nrays-1
        if (fieldchoice.eq.'ISO') then
          pdr(point_id)%rad_surface(ray_index) = Gext(1)/real(nrays,kind=dp)
        else if (fieldchoice.eq.'UNI') then
          if (ray_index.eq.selected_ray_index) then
            pdr(point_id)%rad_surface(ray_index) = -ray_dot_product(Gext(:),vectors(:,ray_index))
          endif
        endif
      enddo

      do ray_index=0,nrays-1
        pdr(point_id)%UVfield = pdr(point_id)%UVfield + &
            &attenuated_uv(pdr(point_id)%rad_surface(ray_index), pdr(point_id)%AV(ray_index), UV_fac)
      enddo
    end subroutine calculate_point_uvfield

  end subroutine calc_uvfield

end module radiation_module
