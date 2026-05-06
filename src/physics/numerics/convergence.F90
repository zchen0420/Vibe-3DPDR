module convergence_module
  use definitions, only : dp
  use healpix_types, only : i4b
  implicit none

contains

  subroutine set_lte_populations(point_id, gas_temperature, gas_density, partition_functions)
    use maincode_module
    use coolants_module, only : coolant_cii, coolant_ci, coolant_oi, coolant_c12o
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

    call calculate_partition_function(cii_partition_function,coolant(coolant_cii)%nlevels,&
        &coolant(coolant_cii)%energies,coolant(coolant_cii)%weights,gas_temperature)
    call calculate_partition_function(ci_partition_function,coolant(coolant_ci)%nlevels,&
        &coolant(coolant_ci)%energies,coolant(coolant_ci)%weights,gas_temperature)
    call calculate_partition_function(oi_partition_function,coolant(coolant_oi)%nlevels,&
        &coolant(coolant_oi)%energies,coolant(coolant_oi)%weights,gas_temperature)
    call calculate_partition_function(c12o_partition_function,coolant(coolant_c12o)%nlevels,&
        &coolant(coolant_c12o)%energies,coolant(coolant_c12o)%weights,gas_temperature)

    if (present(partition_functions)) then
      partition_functions(1) = cii_partition_function
      partition_functions(2) = ci_partition_function
      partition_functions(3) = oi_partition_function
      partition_functions(4) = c12o_partition_function
    end if

    call calculate_lte_populations(coolant(coolant_cii)%nlevels,&
        &grid%points(point_id)%coolant_state(coolant_cii)%population,coolant(coolant_cii)%energies,&
        &coolant(coolant_cii)%weights,cii_partition_function,grid%points(point_id)%abundance(species_idx%ncx)*gas_density,gas_temperature)
    call calculate_lte_populations(coolant(coolant_ci)%nlevels,&
        &grid%points(point_id)%coolant_state(coolant_ci)%population,coolant(coolant_ci)%energies,&
        &coolant(coolant_ci)%weights,ci_partition_function,grid%points(point_id)%abundance(species_idx%nc)*gas_density,gas_temperature)
    call calculate_lte_populations(coolant(coolant_oi)%nlevels,&
        &grid%points(point_id)%coolant_state(coolant_oi)%population,coolant(coolant_oi)%energies,&
        &coolant(coolant_oi)%weights,oi_partition_function,grid%points(point_id)%abundance(species_idx%no)*gas_density,gas_temperature)
    call calculate_lte_populations(coolant(coolant_c12o)%nlevels,&
        &grid%points(point_id)%coolant_state(coolant_c12o)%population,coolant(coolant_c12o)%energies,&
        &coolant(coolant_c12o)%weights,c12o_partition_function,grid%points(point_id)%abundance(species_idx%nco)*gas_density,gas_temperature)
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
  relch_values = 0.0d0

  do level_index=1,nlevels
    if (solution(level_index).ge.abundance*1.0d-10) then
      if (solution(level_index).eq.0.0d0 .and. current_population(level_index).eq.0.0d0) then
        relative_change = 0.0d0
      else
        relative_change = 2.0d0*abs((solution(level_index)-current_population(level_index))&
            &/(solution(level_index)+current_population(level_index)))
      end if

      relch_values(level_index) = relative_change
      if (relative_change.gt.1.0d-2) population_has_converged = .false.
    end if
  end do
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
  end if
  write(*,'(" [",I6,"/",I6,"]")') converged_count,total_count
end subroutine print_convergence_count

end module convergence_module
