module evolution_setup_module
  use definitions
  use chemistry_module
  use maincode_module
  use global_module
  use columns_module, only : calc_columndens
  use dark_region_module, only : dark_molecular_region
  use iteration_chemistry_module, only : run_initial_chemistry_iterations
  use level_population_solver_module, only : set_initial_lte_populations
  use radiation_module, only : calc_uvfield

  implicit none

contains

  subroutine initialise_evolution_state(dobinarychop, previouschange)
    logical, allocatable, intent(inout) :: dobinarychop(:)
    character(len=1), allocatable, intent(inout) :: previouschange(:)

    write(6,*) 'Calculating UV field ...'

    call calc_uvfield

#ifdef DUST2
    call calculate_dust_temperatures
    do pp=1,pdr_ptot
      p=IDlist_pdr(pp)
      dusttemperature(pp)=pdr(p)%dust_t
    enddo
#endif

#ifdef THERMALBALANCE
    allocate(converged(0:pdr_ptot))
    allocate(doleveltmin(0:pdr_ptot))
    doleveltmin=.false.
    converged=.false.
    level_conv=.false.
    first_time=.true.

    if (allocated(dobinarychop)) deallocate(dobinarychop)
    if (allocated(previouschange)) deallocate(previouschange)
    allocate(dobinarychop(0:pdr_ptot))
    dobinarychop=.false.
    allocate(previouschange(0:pdr_ptot))
#endif

    allocate(column(0:pdr_ptot))
    write(6,*) 'Calculating column densities...'
    referee=0
    call calc_columndens
    referee=1

    start_time = 0.0D0

    if (dark_ptot.gt.0) call dark_molecular_region

    iteration = 0
  end subroutine initialise_evolution_state

  subroutine prepare_initial_lte_populations
    real(kind=dp) :: start_time_lte
    real(kind=dp) :: end_time_lte

    write(6,*) ''
    write(6,*) 'Calculating LTE level populations...'
    call cpu_time(t3b)

    relative_abundance_tolerance = 1.0D-8
    absolute_abundance_tolerance = 1.0D-30

    call cpu_time(start_time_lte)
    call run_initial_chemistry_iterations
    call set_initial_lte_populations
    call cpu_time(end_time_lte)

    write(6,*) 'Time for LTE level populations (SERIAL)= ', &
        &(end_time_lte-start_time_lte),' seconds'

    call cpu_time(t3)
    write(6,*) 'Total time = ',t3,' seconds'
    write(6,*) ''
  end subroutine prepare_initial_lte_populations

end module evolution_setup_module
