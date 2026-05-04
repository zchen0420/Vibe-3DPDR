module particle_storage_module
  use healpix_types
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
    do p=1,grand_ptot
      allocate(pdr(p)%CII_pop(1:CII_nlev))
      allocate(pdr(p)%CI_pop(1:CI_nlev))
      allocate(pdr(p)%OI_pop(1:OI_nlev))
      allocate(pdr(p)%C12O_pop(1:C12O_nlev))
    enddo
  end subroutine allocate_level_populations

  subroutine allocate_radiation_arrays
    do pp=1,pdr_ptot
      p=IDlist_pdr(pp)
      call allocate_point_radiation_arrays(p)
    enddo
  end subroutine allocate_radiation_arrays

  subroutine allocate_dark_region_radiation_arrays
    if (dark_ptot.gt.0) then
      call allocate_point_radiation_arrays(IDlist_dark(1))
    endif
  end subroutine allocate_dark_region_radiation_arrays

  subroutine allocate_point_radiation_arrays(point_id)
    integer, intent(in) :: point_id

    allocate(pdr(point_id)%epray(0:nrays-1))
    allocate(pdr(point_id)%epoint(1:3,0:nrays-1,0:maxpoints))
    allocate(pdr(point_id)%projected(0:nrays-1,0:maxpoints))
    allocate(pdr(point_id)%columndensity(0:nrays-1))
    allocate(pdr(point_id)%AV(0:nrays-1))
    allocate(pdr(point_id)%rad_surface(0:nrays-1))
    allocate(pdr(point_id)%CII_line(1:CII_nlev,1:CII_nlev))
    allocate(pdr(point_id)%CI_line(1:CI_nlev,1:CI_nlev))
    allocate(pdr(point_id)%OI_line(1:OI_nlev,1:OI_nlev))
    allocate(pdr(point_id)%C12O_line(1:C12O_nlev,1:C12O_nlev))
    allocate(pdr(point_id)%raytype(0:nrays-1))
    allocate(pdr(point_id)%CII_optdepth(1:CII_nlev,1:CII_nlev,0:nrays-1))
    allocate(pdr(point_id)%CI_optdepth(1:CI_nlev,1:CI_nlev,0:nrays-1))
    allocate(pdr(point_id)%OI_optdepth(1:OI_nlev,1:OI_nlev,0:nrays-1))
    allocate(pdr(point_id)%C12O_optdepth(1:C12O_nlev,1:C12O_nlev,0:nrays-1))
  end subroutine allocate_point_radiation_arrays

  subroutine allocate_solution_arrays
    allocate(CII_solution(0:pdr_ptot,1:CII_nlev))
    allocate(CI_solution(0:pdr_ptot,1:CI_nlev))
    allocate(OI_solution(0:pdr_ptot,1:OI_nlev))
    allocate(C12O_solution(0:pdr_ptot,1:C12O_nlev))
  end subroutine allocate_solution_arrays

end module particle_storage_module
