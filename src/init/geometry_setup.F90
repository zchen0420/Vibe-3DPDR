module geometry_setup_module
  use healpix_types
  use healpix_module
  use maincode_module

contains

  subroutine initialise_pdr_order
    call heapsort(pdr_ptot,rrb,rra)
  end subroutine initialise_pdr_order

  subroutine initialise_healpix_geometry
    nside=2**level
    nrays=12*nside**2
    ns_max=8192

    allocate(vector(1:3))
    allocate(vertex(1:3,1:4))
    allocate(vectors(1:3,0:nrays-1))

    call mk_xy2pix
    call build_healpix_vectors
  end subroutine initialise_healpix_geometry

  subroutine build_healpix_vectors
    write(6,*) 'Building HEALPix vectors...'
    open(unit=77,file='HEALPix_vectors.dat',status='replace')
    do i=1,nrays
      ipix=i-1 !ipix is the ID of a HEALPix ray. Runs with values 0:nrays-1
      call pix2vec_nest(nside,ipix,pix2x,pix2y,vector,vertex)
      vectors(1:3,ipix)=vector(1:3) !Store in memory
      write(77,'(3ES11.3,I7)') vectors(1:3,ipix),ipix
    enddo
    close(77)

    write(6,*) 'Done!';write(6,*) ''
  end subroutine build_healpix_vectors

  subroutine build_evaluation_geometry
    write(6,*) 'Building evaluation points...'
    call evaluation_points

#ifndef PSEUDO_1D
    call trim_boundary_evaluation_points
#endif

    call cpu_time(time_evalpoints)
    write(6,*) 'Time = ',time_evalpoints,' seconds'

    call update_maxpoints_from_evaluation_geometry
  end subroutine build_evaluation_geometry

  subroutine trim_boundary_evaluation_points
    integer(kind=i4b) :: ray_index

    write(6,*) 'subtracting 1 evaluation point in raytype(j)=-2'
    do pp=1,pdr_ptot
      p=IDlist_pdr(pp)
      do ray_index=0,nrays-1
        if (pdr(p)%epray(ray_index).eq.0) cycle
        if (pdr(p)%raytype(ray_index).eq.-2) pdr(p)%epray(ray_index)=pdr(p)%epray(ray_index)-1
      enddo
    enddo
  end subroutine trim_boundary_evaluation_points

  subroutine update_maxpoints_from_evaluation_geometry
    integer(kind=i4b) :: newmaxpoints

    maxpoints = 0
    do pp=1,pdr_ptot
      p=IDlist_pdr(pp)
      newmaxpoints = maxval(pdr(p)%epray)
      if (newmaxpoints.gt.maxpoints) maxpoints = newmaxpoints
    enddo

    if (dark_ptot.gt.0) then
      p=IDlist_dark(1)
      newmaxpoints = maxval(pdr(p)%epray)
      if (newmaxpoints.gt.maxpoints) maxpoints = newmaxpoints
    endif

    write(6,*) '';write(6,*) 'new maxpoints = ',maxpoints
  end subroutine update_maxpoints_from_evaluation_geometry

end module geometry_setup_module
