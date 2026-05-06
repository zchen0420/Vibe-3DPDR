module memory_module
  use definitions, only : dp
  use healpix_types, only : i4b
  use chemistry_network_module, only : allocate_chemistry_state
  use thermal_state_module, only : allocate_temperature_state, allocate_thermal_balance_state
  use coolants_module, only : COLLISION_PARTNER_COUNT, COOLANT_COUNT, COOLANT_CII, COOLANT_CI, COOLANT_OI, &
      &COOLANT_C12O, coolant_data

contains

  subroutine allocations

    use maincode_module
    use uclpdr_module, only : SCO_GRID
    use global_module, only : all_heating
    implicit none
    integer(kind=i4b) :: coolant_id

    !load SCO_GRID data [UCL_PDR]
    SCO_GRID(1:8,1) = (/0.000D+00,-1.408D-02,-1.099D-01,-4.400D-01,-1.154D+00,-1.888D+00,-2.760D+00,-4.001D+00/)
    SCO_GRID(1:8,2) = (/-8.539D-02,-1.015D-01,-2.104D-01,-5.608D-01,-1.272D+00,-1.973D+00,-2.818D+00,-4.055D+00/)
    SCO_GRID(1:8,3) = (/-1.451D-01,-1.612D-01,-2.708D-01,-6.273D-01,-1.355D+00,-2.057D+00,-2.902D+00,-4.122D+00/)
    SCO_GRID(1:8,4) = (/-4.559D-01,-4.666D-01,-5.432D-01,-8.665D-01,-1.602D+00,-2.303D+00,-3.146D+00,-4.421D+00/)
    SCO_GRID(1:8,5) = (/-1.303D+00,-1.312D+00,-1.367D+00,-1.676D+00,-2.305D+00,-3.034D+00,-3.758D+00,-5.077D+00/)
    SCO_GRID(1:8,6) = (/-3.883D+00,-3.888D+00,-3.936D+00,-4.197D+00,-4.739D+00,-5.165D+00,-5.441D+00,-6.446D+00/)

    call allocate_coolant_storage(COOLANT_CII)
    call allocate_coolant_storage(COOLANT_CI)
    call allocate_coolant_storage(COOLANT_OI)
    call allocate_coolant_storage(COOLANT_C12O)

    call allocate_chemistry_state(chemistry, nspec, nreac)

    !allocations start from 0 to cope with the ONE dark molecular element
    if (.not.allocated(coolant_iteration)) allocate(coolant_iteration(1:COOLANT_COUNT))
    call allocate_temperature_state(thermal, pdr_ptot)
    do coolant_id = 1, COOLANT_COUNT
      allocate(coolant_iteration(coolant_id)%cooling_rate(0:pdr_ptot))
    enddo

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
    if (coolant_id.lt.COOLANT_CII .or. coolant_id.gt.COOLANT_C12O) then
      stop 'Invalid coolant storage id'
    endif

    call allocate_coolant_arrays(coolant(coolant_id))
  end subroutine allocate_coolant_storage

  subroutine allocate_coolant_arrays(coolant_table)
    type(coolant_data), intent(inout) :: coolant_table

    allocate(coolant_table%energies(1:coolant_table%nlevels))
    allocate(coolant_table%weights(1:coolant_table%nlevels))
    allocate(coolant_table%a_coeffs(1:coolant_table%nlevels,1:coolant_table%nlevels))
    allocate(coolant_table%b_coeffs(1:coolant_table%nlevels,1:coolant_table%nlevels))
    allocate(coolant_table%frequencies(1:coolant_table%nlevels,1:coolant_table%nlevels))
    allocate(coolant_table%collision_temperatures(1:COLLISION_PARTNER_COUNT,1:coolant_table%ntemperatures))
    allocate(coolant_table%collision_rates(1:COLLISION_PARTNER_COUNT,1:coolant_table%nlevels, &
        &1:coolant_table%nlevels,1:coolant_table%ntemperatures))
  end subroutine allocate_coolant_arrays

end module memory_module
