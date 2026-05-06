module geometry_setup_module
  use healpix_types
  use healpix_module
  use geometry_state_module, only : allocate_ray_geometry
  use maincode_module

contains

  subroutine initialise_pdr_order
    call heapsort(pdr_ptot,geometry%point_order,geometry%radial_order)
  end subroutine initialise_pdr_order

  subroutine initialise_healpix_geometry
    nside=2**level
    nrays=12*nside**2
    ns_max=8192

    allocate(vector(1:3))
    allocate(vertex(1:3,1:4))
    call allocate_ray_geometry(geometry, nrays)

    call mk_xy2pix
    call build_healpix_vectors
  end subroutine initialise_healpix_geometry

  subroutine build_healpix_vectors
    integer(kind=i4b) :: ray_number, pixel_index

    write(6,*) 'Building HEALPix geometry%ray_vectors...'
    open(unit=77,file='HEALPix_vectors.dat',status='replace')
    do ray_number=1,nrays
      pixel_index=ray_number-1 !pixel_index is the ID of a HEALPix ray. Runs with values 0:nrays-1
      call pix2vec_nest(nside,pixel_index,pix2x,pix2y,vector,vertex)
      geometry%ray_vectors(1:3,pixel_index)=vector(1:3) !Store in memory
      write(77,'(3ES11.3,I7)') geometry%ray_vectors(1:3,pixel_index),pixel_index
    end do
    close(77)

    write(6,*) 'Done!';write(6,*) ''
  end subroutine build_healpix_vectors

  subroutine build_evaluation_geometry
    real :: evaluation_time

    write(6,*) 'Building evaluation points...'
    call evaluation_points

#ifndef PSEUDO_1D
    call trim_boundary_evaluation_points
#endif

    call cpu_time(evaluation_time)
    write(6,*) 'Time = ',evaluation_time,' seconds'

    call update_maxpoints_from_evaluation_geometry
  end subroutine build_evaluation_geometry

  subroutine trim_boundary_evaluation_points
    integer(kind=i4b) :: point_id, point_index
    integer(kind=i4b) :: ray_index

    write(6,*) 'subtracting 1 evaluation point in raytype(j)=-2'
    do point_index=1,pdr_ptot
      point_id=grid%pdr_ids(point_index)
      do ray_index=0,nrays-1
        if (grid%points(point_id)%epray(ray_index).eq.0) cycle
        if (grid%points(point_id)%raytype(ray_index).eq.-2) grid%points(point_id)%epray(ray_index)=grid%points(point_id)%epray(ray_index)-1
      end do
    end do
  end subroutine trim_boundary_evaluation_points

  subroutine update_maxpoints_from_evaluation_geometry
    integer(kind=i4b) :: point_id, point_index
    integer(kind=i4b) :: newmaxpoints

    maxpoints = 0
    do point_index=1,pdr_ptot
      point_id=grid%pdr_ids(point_index)
      newmaxpoints = maxval(grid%points(point_id)%epray)
      if (newmaxpoints.gt.maxpoints) maxpoints = newmaxpoints
    end do

    if (dark_ptot.gt.0) then
      point_id=grid%dark_ids(1)
      newmaxpoints = maxval(grid%points(point_id)%epray)
      if (newmaxpoints.gt.maxpoints) maxpoints = newmaxpoints
    end if

    write(6,*) '';write(6,*) 'new maxpoints = ',maxpoints
  end subroutine update_maxpoints_from_evaluation_geometry

end module geometry_setup_module
