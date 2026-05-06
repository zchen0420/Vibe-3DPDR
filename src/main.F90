program threedpdr

  use maincode_module
  use chemistry_io_module, only : load_chemistry_network
  use coolant_io_module, only : load_coolant_data
  use evolution_setup_module, only : initialise_evolution_state, prepare_initial_lte_populations
  use geometry_setup_module, only : initialise_pdr_order, initialise_healpix_geometry, build_evaluation_geometry
  use grid_io_module, only : read_initial_grid
  use initial_conditions_module, only : initialise_temperatures, initialise_particle_abundances
  use iteration_chemistry_module, only : refresh_chemistry_after_temperature_change
  use iteration_convergence_module, only : initialise_level_convergence_state, evaluate_iteration_convergence, &
      &report_iteration_convergence, update_population_densities, write_iteration_timing
  use level_population_solver_module, only : refresh_lte_populations_for_unconverged_points, &
      &solve_level_populations_lvg
  use memory_module, only : allocations
  use output_module, only : write_final_outputs, write_simulation_finished
  use particle_storage_module, only : allocate_particle_radiation_state
  use runtime_config_module, only : load_runtime_configuration
  use runtime_state_module, only : initialize_runtime_config_state
  use spatial_index_module, only : prepare_spatial_index_inputs
  use thermal_balance_module, only : calculate_heating_and_temperature_updates

  implicit none

  character(len=128) :: config_file
  real :: simulation_start_time

  call initialize_runtime_config_state(runtime)
  call cpu_time(simulation_start_time)
  call read_simulation_files
  call initialise_simulation
  call evolve_simulation

contains

  subroutine read_simulation_files
    call load_runtime_configuration(config_file)
    call read_initial_grid
    call prepare_spatial_index_inputs
  end subroutine read_simulation_files

  subroutine initialise_simulation
    call allocations
    call load_coolant_data
    call load_chemistry_network
    call initialise_temperatures
    call initialise_particle_abundances
    call initialise_pdr_order
    call initialise_healpix_geometry
    call allocate_particle_radiation_state
    call build_evaluation_geometry
  end subroutine initialise_simulation

  subroutine evolve_simulation
    logical :: stop_iterations
    logical :: atomic_coolants_converged
    logical, allocatable :: dobinarychop(:)
    character(len=1), allocatable :: previouschange(:)
    real :: iteration_start_time

    call initialise_evolution_state(dobinarychop, previouschange)
    call prepare_initial_lte_populations
    call initialise_level_convergence_state

    write(6,*) 'Begin iterations...'

    levpop_iteration=0
    atomic_coolants_converged=.false.

    do iteration=1,runtime%total_iterations
      write(6,*) ''
      write(6,*) 'Iteration ',iteration
      call cpu_time(iteration_start_time)

      levpop_iteration=levpop_iteration+1
      write(6,*) 'Level population iteration ',levpop_iteration

      if (iteration.gt.1.and.levpop_iteration.eq.1) then
        call refresh_chemistry_after_temperature_change
        call refresh_lte_populations_for_unconverged_points
      else
        call solve_level_populations_lvg(atomic_coolants_converged)
      endif

      call calculate_heating_and_temperature_updates(dobinarychop, previouschange)
      call evaluate_iteration_convergence(stop_iterations)
      if (stop_iterations) exit

      call report_iteration_convergence(atomic_coolants_converged)
      call update_population_densities
      call write_iteration_timing(iteration_start_time)
    enddo

    call write_final_outputs(config_file)
    call write_simulation_finished(simulation_start_time)
  end subroutine evolve_simulation

end program threedpdr
