module convergence_module
  use definitions, only : dp
  use healpix_types, only : i4b
  implicit none

contains

  subroutine set_lte_populations(point_id, gas_temperature, gas_density)
    use maincode_module
    use global_module, only : NCx, NC, NO, NCO
    use excitation_module, only : calculate_partition_function, calculate_lte_populations

    integer(kind=i4b), intent(in) :: point_id
    real(kind=dp), intent(in) :: gas_temperature
    real(kind=dp), intent(in) :: gas_density

    call calculate_partition_function(CII_Z_FUNCTION,CII_NLEV,CII_ENERGIES,CII_WEIGHTS,gas_temperature)
    call calculate_partition_function(CI_Z_FUNCTION,CI_NLEV,CI_ENERGIES,CI_WEIGHTS,gas_temperature)
    call calculate_partition_function(OI_Z_FUNCTION,OI_NLEV,OI_ENERGIES,OI_WEIGHTS,gas_temperature)
    call calculate_partition_function(C12O_Z_FUNCTION,C12O_NLEV,C12O_ENERGIES,C12O_WEIGHTS,gas_temperature)

    call calculate_lte_populations(CII_NLEV,pdr(point_id)%CII_POP,CII_ENERGIES,&
        &CII_WEIGHTS,CII_Z_FUNCTION,pdr(point_id)%ABUNDANCE(NCx)*gas_density,gas_temperature)
    call calculate_lte_populations(CI_NLEV,pdr(point_id)%CI_POP,CI_ENERGIES,&
        &CI_WEIGHTS,CI_Z_FUNCTION,pdr(point_id)%ABUNDANCE(NC)*gas_density,gas_temperature)
    call calculate_lte_populations(OI_NLEV,pdr(point_id)%OI_POP,OI_ENERGIES,&
        &OI_WEIGHTS,OI_Z_FUNCTION,pdr(point_id)%ABUNDANCE(NO)*gas_density,gas_temperature)
    call calculate_lte_populations(C12O_NLEV,pdr(point_id)%C12O_POP,C12O_ENERGIES,&
        &C12O_WEIGHTS,C12O_Z_FUNCTION,pdr(point_id)%ABUNDANCE(NCO)*gas_density,gas_temperature)
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
      write(*,'(" ",A," are ",F5.1,"% converged.")') trim(label), percentage
    else
      write(*,'(" ",A," is ",F5.1,"% converged.")') trim(label), percentage
    endif
    write(*,'(" [",I6,"/",I6,"]")') converged_count,total_count
  end subroutine print_convergence_count

end module convergence_module
