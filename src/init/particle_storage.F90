module particle_storage_module
  use healpix_types
  use coolants_module, only : coolant_count, coolant_cii, coolant_ci, coolant_oi, coolant_c12o
  use maincode_module

contains

  subroutine allocate_particle_radiation_state
    write(6,*) 'allocating memory...'

    call allocate_level_populations
    call allocate_radiation_arrays
    call allocate_dark_region_radiation_arrays
    call allocate_solution_arrays

    write(6,*) 'Memory OK';write(6,*) ''
  end subroutine allocate_particle_radiation_state

  subroutine allocate_level_populations
    integer(kind=i4b) :: point_id

    do point_id=1,grand_ptot
      allocate(grid%points(point_id)%coolant_state(1:coolant_count))
      call allocate_point_coolant_population(point_id, coolant_cii)
      call allocate_point_coolant_population(point_id, coolant_ci)
      call allocate_point_coolant_population(point_id, coolant_oi)
      call allocate_point_coolant_population(point_id, coolant_c12o)
    end do
  end subroutine allocate_level_populations

  subroutine allocate_radiation_arrays
    integer(kind=i4b) :: point_id, point_index

    do point_index=1,pdr_ptot
      point_id=grid%pdr_ids(point_index)
      call allocate_point_radiation_arrays(point_id)
    end do
  end subroutine allocate_radiation_arrays

  subroutine allocate_dark_region_radiation_arrays
    if (dark_ptot.gt.0) then
      call allocate_point_radiation_arrays(grid%dark_ids(1))
    end if
  end subroutine allocate_dark_region_radiation_arrays

  subroutine allocate_point_radiation_arrays(point_id)
    integer, intent(in) :: point_id

    allocate(grid%points(point_id)%epray(0:nrays-1))
    allocate(grid%points(point_id)%epoint(1:3,0:nrays-1,0:maxpoints))
    allocate(grid%points(point_id)%projected(0:nrays-1,0:maxpoints))
    allocate(grid%points(point_id)%columndensity(0:nrays-1))
    allocate(grid%points(point_id)%av(0:nrays-1))
    allocate(grid%points(point_id)%rad_surface(0:nrays-1))
    allocate(grid%points(point_id)%raytype(0:nrays-1))
    call allocate_point_coolant_radiation(point_id, coolant_cii)
    call allocate_point_coolant_radiation(point_id, coolant_ci)
    call allocate_point_coolant_radiation(point_id, coolant_oi)
    call allocate_point_coolant_radiation(point_id, coolant_c12o)
  end subroutine allocate_point_radiation_arrays

  subroutine allocate_point_coolant_population(point_id, coolant_id)
    integer(kind=i4b), intent(in) :: point_id
    integer(kind=i4b), intent(in) :: coolant_id
    integer(kind=i4b) :: level_count

    level_count = coolant(coolant_id)%nlevels
    allocate(grid%points(point_id)%coolant_state(coolant_id)%population(1:level_count))
  end subroutine allocate_point_coolant_population

  subroutine allocate_point_coolant_radiation(point_id, coolant_id)
    integer(kind=i4b), intent(in) :: point_id
    integer(kind=i4b), intent(in) :: coolant_id
    integer(kind=i4b) :: level_count

    level_count = coolant(coolant_id)%nlevels
    allocate(grid%points(point_id)%coolant_state(coolant_id)%line(1:level_count,1:level_count))
    allocate(grid%points(point_id)%coolant_state(coolant_id)%optical_depth(1:level_count,1:level_count,0:nrays-1))
  end subroutine allocate_point_coolant_radiation

  subroutine allocate_solution_arrays
    integer(kind=i4b) :: coolant_id

    if (.not.allocated(coolant_iteration)) allocate(coolant_iteration(1:coolant_count))
    do coolant_id = 1, coolant_count
      allocate(coolant_iteration(coolant_id)%solution(0:pdr_ptot,1:coolant(coolant_id)%nlevels))
    end do
  end subroutine allocate_solution_arrays

end module particle_storage_module
