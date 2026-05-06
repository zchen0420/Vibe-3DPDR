module convergence_module
  use definitions, only : dp
  use healpix_types, only : i4b
  implicit none

contains

  subroutine set_lte_populations(point_id, gas_temperature, gas_density, partition_functions)
    use maincode_module
    use coolants_module, only : COOLANT_CII, COOLANT_CI, COOLANT_OI, COOLANT_C12O
    use global_module, only : species_idx
    use excitation_module, only : calculate_partition_function, calculate_lte_populations

    integer(kind=i4b), intent(in) :: point_id
    real(kind=dp), intent(in) :: gas_temperature
    real(kind=dp), intent(in) :: gas_density
    real(kind=dp), intent(out), optional :: partition_functions(4)
    real(kind=dp) :: cii_partition_function
    real(kind=dp) :: ci_partition_function
    real(kind=dp) :: oi_partition_function
    real(kind=dp) :: c12o_partition_function

    call calculate_partition_function(cii_partition_function,coolant(COOLANT_CII)%nlevels,&
        &coolant(COOLANT_CII)%energies,coolant(COOLANT_CII)%weights,gas_temperature)
    call calculate_partition_function(ci_partition_function,coolant(COOLANT_CI)%nlevels,&
        &coolant(COOLANT_CI)%energies,coolant(COOLANT_CI)%weights,gas_temperature)
    call calculate_partition_function(oi_partition_function,coolant(COOLANT_OI)%nlevels,&
        &coolant(COOLANT_OI)%energies,coolant(COOLANT_OI)%weights,gas_temperature)
    call calculate_partition_function(c12o_partition_function,coolant(COOLANT_C12O)%nlevels,&
        &coolant(COOLANT_C12O)%energies,coolant(COOLANT_C12O)%weights,gas_temperature)

    if (present(partition_functions)) then
      partition_functions(1) = cii_partition_function
      partition_functions(2) = ci_partition_function
      partition_functions(3) = oi_partition_function
      partition_functions(4) = c12o_partition_function
    endif

    call calculate_lte_populations(coolant(COOLANT_CII)%nlevels,&
        &grid%points(point_id)%coolant_state(COOLANT_CII)%population,coolant(COOLANT_CII)%energies,&
        &coolant(COOLANT_CII)%weights,cii_partition_function,grid%points(point_id)%ABUNDANCE(species_idx%NCx)*gas_density,gas_temperature)
    call calculate_lte_populations(coolant(COOLANT_CI)%nlevels,&
        &grid%points(point_id)%coolant_state(COOLANT_CI)%population,coolant(COOLANT_CI)%energies,&
        &coolant(COOLANT_CI)%weights,ci_partition_function,grid%points(point_id)%ABUNDANCE(species_idx%NC)*gas_density,gas_temperature)
    call calculate_lte_populations(coolant(COOLANT_OI)%nlevels,&
        &grid%points(point_id)%coolant_state(COOLANT_OI)%population,coolant(COOLANT_OI)%energies,&
        &coolant(COOLANT_OI)%weights,oi_partition_function,grid%points(point_id)%ABUNDANCE(species_idx%NO)*gas_density,gas_temperature)
    call calculate_lte_populations(coolant(COOLANT_C12O)%nlevels,&
        &grid%points(point_id)%coolant_state(COOLANT_C12O)%population,coolant(COOLANT_C12O)%energies,&
        &coolant(COOLANT_C12O)%weights,c12o_partition_function,grid%points(point_id)%ABUNDANCE(species_idx%NCO)*gas_density,gas_temperature)
  end subroutine set_lte_populations

  logical function population_has_converged(solution, current_population, abundance, relch_values, nlevels)
    integer(kind=i4b), intent(in) :: nlevels
    real(kind=dp), intent(in) :: solution(1:nlevels)
    real(kind=dp), intent(in) :: current_population(1:nlevels)
    real(kind=dp), intent(in) :: abundance
    real(kind=dp), intent(out) :: relch_values(1:nlevels)

    integer(kind=i4b) :: level_index
    real(kind=dp) :: relative_change

    population_has_converged = .true.
    relch_values = 0.0D0

    do level_index=1,nlevels
      if (solution(level_index).ge.abundance*1.0D-10) then
        if (solution(level_index).eq.0.0D0 .and. current_population(level_index).eq.0.0D0) then
          relative_change = 0.0D0
        else
          relative_change = 2.0D0*abs((solution(level_index)-current_population(level_index))&
              &/(solution(level_index)+current_population(level_index)))
        endif

        relch_values(level_index) = relative_change
        if (relative_change.gt.1.0D-2) population_has_converged = .false.
      endif
    enddo
  end function population_has_converged

  subroutine print_convergence_count(label, percentage, converged_count, total_count)
    character(len=*), intent(in) :: label
    real(kind=dp), intent(in) :: percentage
    integer(kind=i4b), intent(in) :: converged_count
    integer(kind=i4b), intent(in) :: total_count

    if (trim(label).eq.'Level populations') then
      write(*,'(" ",A," are ",F5.1,"% thermal%thermal_converged.")') trim(label), percentage
    else
      write(*,'(" ",A," is ",F5.1,"% thermal%thermal_converged.")') trim(label), percentage
    endif
    write(*,'(" [",I6,"/",I6,"]")') converged_count,total_count
  end subroutine print_convergence_count

end module convergence_module
