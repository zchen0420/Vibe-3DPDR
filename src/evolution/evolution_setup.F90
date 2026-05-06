module evolution_setup_module
  use definitions, only : dp
  use healpix_types, only : i4b
  use chemistry_module, only : absolute_abundance_tolerance, relative_abundance_tolerance
  use geometry_state_module, only : allocate_column_density_storage
  use thermal_state_module, only : allocate_thermal_convergence_state
  use maincode_module, only : dark_ptot, first_time, geometry, grid, iteration, level_conv, pdr_ptot, start_time, thermal
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
    integer(kind=i4b) :: point_id, point_index

    write(6,*) 'Calculating UV field ...'

    call calc_uvfield

#ifdef DUST2
    call calculate_dust_temperatures
    do point_index=1,pdr_ptot
      point_id=grid%pdr_ids(point_index)
      thermal%dust_temperature(point_index)=grid%points(point_id)%dust_temperature
    end do
#endif

#ifdef THERMALBALANCE
    call allocate_thermal_convergence_state(thermal, pdr_ptot)
    level_conv=.false.
    first_time=.true.

    if (allocated(dobinarychop)) deallocate(dobinarychop)
    if (allocated(previouschange)) deallocate(previouschange)
    allocate(dobinarychop(0:pdr_ptot))
    dobinarychop=.false.
    allocate(previouschange(0:pdr_ptot))
#endif

    call allocate_column_density_storage(geometry, pdr_ptot)
    write(6,*) 'Calculating geometry%column_density densities...'
    call calc_columndens(.true.)

    start_time = 0.0d0

    if (dark_ptot.gt.0) call dark_molecular_region

    iteration = 0
  end subroutine initialise_evolution_state

  subroutine prepare_initial_lte_populations
    real(kind=dp) :: start_time_lte
    real(kind=dp) :: end_time_lte
    real :: total_time

    write(6,*) ''
    write(6,*) 'Calculating LTE level populations...'

    relative_abundance_tolerance = 1.0d-8
    absolute_abundance_tolerance = 1.0d-30

    call cpu_time(start_time_lte)
    call run_initial_chemistry_iterations
    call set_initial_lte_populations
    call cpu_time(end_time_lte)

    write(6,*) 'Time for LTE level populations (SERIAL)= ', &
        &(end_time_lte-start_time_lte),' seconds'

    call cpu_time(total_time)
    write(6,*) 'Total time = ',total_time,' seconds'
    write(6,*) ''
  end subroutine prepare_initial_lte_populations

end module evolution_setup_module
