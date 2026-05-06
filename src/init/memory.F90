module memory_module
  use definitions, only : dp
  use healpix_types, only : i4b
  use chemistry_network_module, only : allocate_chemistry_state
  use thermal_state_module, only : allocate_temperature_state, allocate_thermal_balance_state
  use coolants_module, only : collision_partner_count, coolant_count, coolant_cii, coolant_ci, coolant_oi, &
      &coolant_c12o, coolant_data

contains

  subroutine allocations

    use maincode_module
    use uclpdr_module, only : sco_grid
    use global_module, only : all_heating
    implicit none
    integer(kind=i4b) :: coolant_id

    !load SCO_GRID data [UCL_PDR]
    sco_grid(1:8,1) = (/0.000d+00,-1.408d-02,-1.099d-01,-4.400d-01,-1.154d+00,-1.888d+00,-2.760d+00,-4.001d+00/)
    sco_grid(1:8,2) = (/-8.539d-02,-1.015d-01,-2.104d-01,-5.608d-01,-1.272d+00,-1.973d+00,-2.818d+00,-4.055d+00/)
    sco_grid(1:8,3) = (/-1.451d-01,-1.612d-01,-2.708d-01,-6.273d-01,-1.355d+00,-2.057d+00,-2.902d+00,-4.122d+00/)
    sco_grid(1:8,4) = (/-4.559d-01,-4.666d-01,-5.432d-01,-8.665d-01,-1.602d+00,-2.303d+00,-3.146d+00,-4.421d+00/)
    sco_grid(1:8,5) = (/-1.303d+00,-1.312d+00,-1.367d+00,-1.676d+00,-2.305d+00,-3.034d+00,-3.758d+00,-5.077d+00/)
    sco_grid(1:8,6) = (/-3.883d+00,-3.888d+00,-3.936d+00,-4.197d+00,-4.739d+00,-5.165d+00,-5.441d+00,-6.446d+00/)

    call allocate_coolant_storage(coolant_cii)
    call allocate_coolant_storage(coolant_ci)
    call allocate_coolant_storage(coolant_oi)
    call allocate_coolant_storage(coolant_c12o)

    call allocate_chemistry_state(chemistry, nspec, nreac)

    !allocations start from 0 to cope with the ONE dark molecular element
    if (.not.allocated(coolant_iteration)) allocate(coolant_iteration(1:coolant_count))
    call allocate_temperature_state(thermal, pdr_ptot)
    do coolant_id = 1, coolant_count
      allocate(coolant_iteration(coolant_id)%cooling_rate(0:pdr_ptot))
    end do

    allocate(all_heating(0:pdr_ptot,1:12))

#ifdef THERMALBALANCE
    call allocate_thermal_balance_state(thermal, pdr_ptot)
#endif

    return
  end subroutine allocations

  subroutine allocate_coolant_storage(coolant_id)
    use maincode_module
    integer(kind=i4b), intent(in) :: coolant_id

    if (.not.allocated(coolant)) stop 'Coolant metadata not initialized'
    if (coolant_id.lt.coolant_cii .or. coolant_id.gt.coolant_c12o) then
      stop 'Invalid coolant storage id'
    end if

    call allocate_coolant_arrays(coolant(coolant_id))
  end subroutine allocate_coolant_storage

  subroutine allocate_coolant_arrays(coolant_table)
    type(coolant_data), intent(inout) :: coolant_table

    allocate(coolant_table%energies(1:coolant_table%nlevels))
    allocate(coolant_table%weights(1:coolant_table%nlevels))
    allocate(coolant_table%a_coeffs(1:coolant_table%nlevels,1:coolant_table%nlevels))
    allocate(coolant_table%b_coeffs(1:coolant_table%nlevels,1:coolant_table%nlevels))
    allocate(coolant_table%frequencies(1:coolant_table%nlevels,1:coolant_table%nlevels))
    allocate(coolant_table%collision_temperatures(1:collision_partner_count,1:coolant_table%ntemperatures))
    allocate(coolant_table%collision_rates(1:collision_partner_count,1:coolant_table%nlevels, &
        &1:coolant_table%nlevels,1:coolant_table%ntemperatures))
  end subroutine allocate_coolant_arrays

end module memory_module
