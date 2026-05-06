module grid_io_module
  use definitions
  use healpix_types
  use simulation_grid_module, only : allocate_simulation_grid
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
    integer(kind=i4b) :: record_count
    integer :: read_status
    real(kind=dp) :: dummy_value

    open(unit=2,file=runtime%input_file,status='old')

    record_count=0
    do
      read(2,*,iostat=read_status) dummy_value
      if (read_status.ne.0) exit
      record_count=record_count+1
    end do

    grand_ptot=record_count
    write(6,*) 'Total elements: ',grand_ptot

    close(2)
    write(6,*) ''
  end subroutine count_grid_points

  subroutine allocate_grid_storage
    call allocate_simulation_grid(grid, grand_ptot)
  end subroutine allocate_grid_storage

  subroutine classify_grid_points(maximum_density, minimum_density)
    real(kind=dp), intent(out) :: maximum_density, minimum_density
    integer(kind=i4b) :: point_id
    real(kind=dp) :: x_position, y_position, z_position, density_value

    pdr_ptot=0
    ion_ptot=0
    dark_ptot=0
    maximum_density=0.0d0
    minimum_density=1.0d10

    open(unit=2,file=runtime%input_file,status='old')

    do point_id=1,grand_ptot
      read(2,*) x_position,y_position,z_position,density_value
      if (density_value.le.rho_min) then
        ion_ptot = ion_ptot + 1
        grid%points(point_id)%etype = 2 !IONIZED
        grid%points(point_id)%position = (/x_position, y_position, z_position/)
        grid%points(point_id)%rho=density_value
        grid%ion_ids(ion_ptot)=point_id
      end if
      if ((density_value.gt.rho_min).and.(density_value.le.rho_max)) then
        pdr_ptot = pdr_ptot + 1
        grid%points(point_id)%etype = 1 !grid%points
        grid%points(point_id)%position = (/x_position, y_position, z_position/)
        grid%points(point_id)%rho=density_value
        if (density_value.gt.maximum_density) maximum_density=density_value
        if (density_value.lt.minimum_density) minimum_density=density_value
        grid%pdr_ids(pdr_ptot)=point_id
      end if
      if (density_value.gt.rho_max) then
        dark_ptot = dark_ptot + 1
        grid%points(point_id)%etype = 3 !DARK MOLECULAR
        grid%points(point_id)%position = (/x_position, y_position, z_position/)
        grid%points(point_id)%rho=density_value
        grid%dark_ids(dark_ptot)=point_id
      end if
    end do

    close(2)
  end subroutine classify_grid_points

  subroutine print_grid_classification_summary(maximum_density, minimum_density)
    real(kind=dp), intent(in) :: maximum_density, minimum_density

    write(6,*) 'PDR elements       = ',pdr_ptot
    write(6,*) 'IONIZED elements   = ',ion_ptot
    write(6,*) 'MOLECULAR elements = ',dark_ptot
    write(6,*) 'Maximum PDR density = ',maximum_density
    write(6,*) 'Minimum PDR density = ',minimum_density
    write(6,*) 'Density used in DMR = ',2.0d0*rho_max
  end subroutine print_grid_classification_summary

end module grid_io_module
