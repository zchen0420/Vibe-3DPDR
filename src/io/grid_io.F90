module grid_io_module
  use definitions
  use healpix_types
  use maincode_module

contains

  subroutine read_initial_grid
    real(kind=dp) :: maximum_density, minimum_density

    write(6,*) 'Reading initial conditions file'

    call count_grid_points
    call allocate_grid_storage
    call classify_grid_points(maximum_density, minimum_density)
    call print_grid_classification_summary(maximum_density, minimum_density)
  end subroutine read_initial_grid

  subroutine count_grid_points
    integer :: read_status

    open(unit=2,file=input,status='old')

    i=0
    do
      read(2,*,iostat=read_status) points
      if (read_status.ne.0) exit
      i=i+1
    enddo

    grand_ptot=i
    write(6,*) 'Total elements: ',grand_ptot

    close(2)
    write(6,*) ''
  end subroutine count_grid_points

  subroutine allocate_grid_storage
    allocate(pdr(1:grand_ptot))
    allocate(IDlist_pdr(1:grand_ptot))
    allocate(IDlist_ion(1:grand_ptot))
    allocate(IDlist_dark(1:grand_ptot))
  end subroutine allocate_grid_storage

  subroutine classify_grid_points(maximum_density, minimum_density)
    real(kind=dp), intent(out) :: maximum_density, minimum_density

    pdr_ptot=0
    ion_ptot=0
    dark_ptot=0
    maximum_density=0.0D0
    minimum_density=1.0D10

    open(unit=2,file=input,status='old')

    do p=1,grand_ptot
      read(2,*) xpos,ypos,zpos,denst
      if (denst.le.rho_min) then
        ion_ptot = ion_ptot + 1
        pdr(p)%etype = 2 !IONIZED
        pdr(p)%x=xpos
        pdr(p)%y=ypos
        pdr(p)%z=zpos
        pdr(p)%rho=denst
        IDlist_ion(ion_ptot)=p
      endif
      if ((denst.gt.rho_min).AND.(denst.le.rho_max)) then
        pdr_ptot = pdr_ptot + 1
        pdr(p)%etype = 1 !PDR
        pdr(p)%x=xpos
        pdr(p)%y=ypos
        pdr(p)%z=zpos
        pdr(p)%rho=denst
        if (denst.gt.maximum_density) maximum_density=denst
        if (denst.lt.minimum_density) minimum_density=denst
        IDlist_pdr(pdr_ptot)=p
      endif
      if (denst.gt.rho_max) then
        dark_ptot = dark_ptot + 1
        pdr(p)%etype = 3 !DARK MOLECULAR
        pdr(p)%x=xpos
        pdr(p)%y=ypos
        pdr(p)%z=zpos
        pdr(p)%rho=denst
        IDlist_dark(dark_ptot)=p
      endif
    enddo

    close(2)
  end subroutine classify_grid_points

  subroutine print_grid_classification_summary(maximum_density, minimum_density)
    real(kind=dp), intent(in) :: maximum_density, minimum_density

    write(6,*) 'PDR elements       = ',pdr_ptot
    write(6,*) 'IONIZED elements   = ',ion_ptot
    write(6,*) 'MOLECULAR elements = ',dark_ptot
    write(6,*) 'Maximum PDR density = ',maximum_density
    write(6,*) 'Minimum PDR density = ',minimum_density
    write(6,*) 'Density used in DMR = ',2.0D0*rho_max
  end subroutine print_grid_classification_summary

end module grid_io_module
