module level_population_solver_module
  use definitions, only : dp
  use healpix_types, only : i4b
  use collision_coefficients_module, only : calculate_collision_coefficients
  use coolants_module, only : collider_electron, collider_h, collider_h2, collider_he, collider_ortho_h2, &
      &collider_para_h2, collider_proton, collision_partner_count, coolant_count, coolant_cii, coolant_ci, &
      &coolant_oi, coolant_c12o
  use escape_probability_module, only : calculate_lvg_transition_rates, lvg_local_conditions
  use maincode_module, only : coolant, coolant_iteration, grid, levpop_iteration, maxpoints, nrays, pdr_ptot, runtime, thermal
  use global_module, only : metallicity, nelect, species_idx
  use convergence_module, only : set_lte_populations
  use level_population_system_module, only : solve_statistical_equilibrium
  use ray_path_module, only : projected_point_id, set_ray_origin_projection

  implicit none

  type coolant_work_item
    integer(kind=i4b) :: coolant_id
    integer(kind=i4b) :: species_abundance_index
    integer(kind=i4b) :: level_count
    real(kind=dp), allocatable :: collision_coefficients(:,:)
    real(kind=dp), allocatable :: transition(:,:)
    real(kind=dp), allocatable :: line(:,:)
    real(kind=dp), allocatable :: optical_depth(:,:,:)
    real(kind=dp), allocatable :: beta(:,:,:)
    real(kind=dp), allocatable :: solution(:)
    real(kind=dp), allocatable :: evalpop(:,:,:)
  end type coolant_work_item

contains

  subroutine set_initial_lte_populations
    integer(kind=i4b) :: point_id, point_index
    real(kind=dp) :: partition_functions(4)

    do point_index=1,pdr_ptot
      point_id=grid%pdr_ids(point_index)
      call set_lte_populations(point_id, thermal%gas_temperature(point_index), grid%points(point_id)%rho, partition_functions)
#ifndef GUESS_TEMP
      if (point_index.eq.1) then
        write(6,*) ''
        write(6,*) 'Z(CII)  = ',partition_functions(1)
        write(6,*) 'Z(CI)   = ',partition_functions(2)
        write(6,*) 'Z(OI)   = ',partition_functions(3)
        write(6,*) 'Z(C12O) = ',partition_functions(4)
        write(6,*) ''
      end if
#endif
    end do
  end subroutine set_initial_lte_populations

  subroutine refresh_lte_populations_for_unconverged_points
    integer(kind=i4b) :: point_index
    integer(kind=i4b) :: point_id

    do point_index=1,pdr_ptot
      point_id=grid%pdr_ids(point_index)
#ifdef THERMALBALANCE
      if (thermal%level_population_converged(point_index).or.thermal%thermal_converged(point_index)) cycle
#else
      if (thermal%level_population_converged(point_index)) cycle
#endif
      call set_lte_populations(point_id, thermal%gas_temperature(point_index), grid%points(point_id)%rho)
    end do
  end subroutine refresh_lte_populations_for_unconverged_points

  subroutine initialize_coolant_work(work_item, coolant_id, species_abundance_index)
    type(coolant_work_item), intent(out) :: work_item
    integer(kind=i4b), intent(in) :: coolant_id
    integer(kind=i4b), intent(in) :: species_abundance_index

    work_item%coolant_id = coolant_id
    work_item%species_abundance_index = species_abundance_index
    work_item%level_count = coolant(coolant_id)%nlevels
  end subroutine initialize_coolant_work

  subroutine allocate_coolant_workspace(work_item)
    type(coolant_work_item), intent(inout) :: work_item

    allocate(work_item%collision_coefficients(1:work_item%level_count,1:work_item%level_count))
    allocate(work_item%transition(1:work_item%level_count,1:work_item%level_count))
    allocate(work_item%line(1:work_item%level_count,1:work_item%level_count))
    allocate(work_item%optical_depth(1:work_item%level_count,1:work_item%level_count,0:nrays-1))
    allocate(work_item%beta(1:work_item%level_count,1:work_item%level_count,0:nrays-1))
    allocate(work_item%solution(1:work_item%level_count))
    allocate(work_item%evalpop(0:nrays-1,0:maxpoints,1:work_item%level_count))

    work_item%evalpop=0.0d0
  end subroutine allocate_coolant_workspace

  subroutine deallocate_coolant_workspace(work_item)
    type(coolant_work_item), intent(inout) :: work_item

    deallocate(work_item%collision_coefficients)
    deallocate(work_item%transition)
    deallocate(work_item%line)
    deallocate(work_item%optical_depth)
    deallocate(work_item%beta)
    deallocate(work_item%solution)
    deallocate(work_item%evalpop)
  end subroutine deallocate_coolant_workspace

  subroutine solve_level_populations_lvg(atomic_coolants_converged)
    logical, intent(in) :: atomic_coolants_converged
    integer(kind=i4b) :: point_index
    integer(kind=i4b) :: point_id

    integer(kind=i4b) :: ray_index
    integer(kind=i4b) :: eval_index
    integer(kind=i4b) :: level_index
    integer(kind=i4b) :: coolant_id
    type(coolant_work_item) :: coolant_work(coolant_count)

    call initialize_coolant_work(coolant_work(coolant_cii),coolant_cii,species_idx%ncx)
    call initialize_coolant_work(coolant_work(coolant_ci),coolant_ci,species_idx%nc)
    call initialize_coolant_work(coolant_work(coolant_oi),coolant_oi,species_idx%no)
    call initialize_coolant_work(coolant_work(coolant_c12o),coolant_c12o,species_idx%nco)

    do point_index=1,pdr_ptot
      point_id=grid%pdr_ids(point_index)

#ifdef THERMALBALANCE
      if (thermal%level_population_converged(point_index).or.thermal%thermal_converged(point_index)) cycle
#else
      if (thermal%level_population_converged(point_index)) cycle
#endif

      call set_ray_origin_projection(grid%points(point_id), point_id)

      call allocate_lvg_workspace
      call populate_evaluation_populations(point_id)

      call solve_coolant_population(coolant_work(coolant_cii), point_index, point_id, .true.)
      call solve_coolant_population(coolant_work(coolant_ci), point_index, point_id, .true.)
      call solve_coolant_population(coolant_work(coolant_oi), point_index, point_id, .true.)
      call solve_coolant_population(coolant_work(coolant_c12o), point_index, point_id, atomic_coolants_converged)

      call deallocate_lvg_workspace
    end do

    thermal%total_cooling_rate=0.0d0
    do coolant_id = 1, coolant_count
      thermal%total_cooling_rate = thermal%total_cooling_rate + coolant_iteration(coolant_id)%cooling_rate
    end do

  contains

    subroutine allocate_lvg_workspace
      do coolant_id=1,coolant_count
        call allocate_coolant_workspace(coolant_work(coolant_id))
      end do
    end subroutine allocate_lvg_workspace

    subroutine populate_evaluation_populations(source_point_id)
      integer(kind=i4b), intent(in) :: source_point_id
      integer(kind=i4b) :: eval_point_id

      do ray_index=0,nrays-1
        do eval_index=0,grid%points(source_point_id)%epray(ray_index)
          eval_point_id=projected_point_id(grid%points(source_point_id), ray_index, eval_index)

          do coolant_id=1,coolant_count
            do level_index=1,coolant_work(coolant_id)%level_count
              coolant_work(coolant_id)%evalpop(ray_index,eval_index,level_index)=&
                  &grid%points(eval_point_id)%coolant_state(coolant_id)%population(level_index)
            end do
          end do
        end do
      end do
    end subroutine populate_evaluation_populations

    subroutine solve_coolant_population(work_item, target_point_index, target_point_id, relaxation_enabled)
      type(coolant_work_item), intent(inout) :: work_item
      integer(kind=i4b), intent(in) :: target_point_index
      integer(kind=i4b), intent(in) :: target_point_id
      logical, intent(in) :: relaxation_enabled
      real(kind=dp) :: species_density

      species_density = abundance_density(target_point_id, work_item%species_abundance_index)

      call update_collision_rates(work_item, target_point_index, target_point_id)
      call update_lvg_transition_rates(work_item, target_point_index, target_point_id)
      call solve_statistical_equilibrium(work_item%transition, species_density, work_item%solution)

      coolant_iteration(work_item%coolant_id)%solution(target_point_index,:)=work_item%solution
      call apply_population_relaxation(work_item%coolant_id, target_point_index, target_point_id, relaxation_enabled)
    end subroutine solve_coolant_population

    real(kind=dp) function abundance_density(point_id, abundance_index)
    integer(kind=i4b), intent(in) :: point_id
    integer(kind=i4b), intent(in) :: abundance_index

    abundance_density = grid%points(point_id)%abundance(abundance_index)*grid%points(point_id)%rho
  end function abundance_density

  subroutine update_collision_rates(work_item, target_point_index, target_point_id)
    type(coolant_work_item), intent(inout) :: work_item
    integer(kind=i4b), intent(in) :: target_point_index
    integer(kind=i4b), intent(in) :: target_point_id
    real(kind=dp) :: collider_density(1:collision_partner_count)

    call calculate_collider_densities(thermal%gas_temperature(target_point_index), target_point_id, collider_density)
    call calculate_collision_coefficients(coolant(work_item%coolant_id), thermal%gas_temperature(target_point_index), &
        &collider_density, work_item%collision_coefficients)
  end subroutine update_collision_rates

  subroutine calculate_collider_densities(gas_temperature, point_id, collider_density)
    real(kind=dp), intent(in) :: gas_temperature
    integer(kind=i4b), intent(in) :: point_id
    real(kind=dp), intent(out) :: collider_density(1:collision_partner_count)

    real(kind=dp) :: h2_density
    real(kind=dp) :: ortho_h2_fraction
    real(kind=dp) :: para_h2_fraction

    collider_density = 0.0d0
    h2_density = abundance_density(point_id, species_idx%nh2)

    if (h2_density.gt.0.0d0) then
      para_h2_fraction = 1.0d0/(1.0d0+9.0d0*exp(-170.5d0/gas_temperature))
      ortho_h2_fraction = 1.0d0 - para_h2_fraction
    else
      para_h2_fraction = 0.0d0
      ortho_h2_fraction = 0.0d0
    end if

    collider_density(collider_h2) = h2_density
    collider_density(collider_para_h2) = h2_density*para_h2_fraction
    collider_density(collider_ortho_h2) = h2_density*ortho_h2_fraction
    collider_density(collider_electron) = abundance_density(point_id, nelect)
    collider_density(collider_h) = abundance_density(point_id, species_idx%nh)
    collider_density(collider_he) = abundance_density(point_id, species_idx%nhe)
    collider_density(collider_proton) = abundance_density(point_id, species_idx%nproton)
  end subroutine calculate_collider_densities

  subroutine update_lvg_transition_rates(work_item, target_point_index, target_point_id)
    type(coolant_work_item), intent(inout) :: work_item
    integer(kind=i4b), intent(in) :: target_point_index
    integer(kind=i4b), intent(in) :: target_point_id
    type(lvg_local_conditions) :: conditions

    conditions%gas_temperature = thermal%gas_temperature(target_point_index)
    conditions%dust_temperature = thermal%dust_temperature(target_point_index)
    conditions%turbulent_velocity = runtime%turbulent_velocity
    conditions%gas_density = grid%points(target_point_id)%rho
    conditions%metallicity = metallicity

    call calculate_lvg_transition_rates(coolant(work_item%coolant_id), conditions, &
        &grid%points(target_point_id)%epray, &
        &grid%points(target_point_id)%coolant_state(work_item%coolant_id)%population, &
        &grid%points(target_point_id)%epoint, work_item%evalpop, work_item%collision_coefficients, &
        &work_item%transition, coolant_iteration(work_item%coolant_id)%cooling_rate(target_point_index), &
        &work_item%line, work_item%optical_depth, work_item%beta)

    grid%points(target_point_id)%coolant_state(work_item%coolant_id)%line = work_item%line
    grid%points(target_point_id)%coolant_state(work_item%coolant_id)%optical_depth = work_item%optical_depth
  end subroutine update_lvg_transition_rates

  subroutine apply_population_relaxation(coolant_id, target_point_index, target_point_id, relaxation_enabled)
    integer(kind=i4b), intent(in) :: coolant_id
    integer(kind=i4b), intent(in) :: target_point_index
    integer(kind=i4b), intent(in) :: target_point_id
    logical, intent(in) :: relaxation_enabled

#ifdef CO_FIX
    if (.not.relaxation_enabled) return

    if (levpop_iteration.ge.120) then
      coolant_iteration(coolant_id)%solution(target_point_index,:) = &
          &grid%points(target_point_id)%coolant_state(coolant_id)%population
    else if (levpop_iteration.ge.75) then
      coolant_iteration(coolant_id)%solution(target_point_index,:) = &
          &0.5d0*(coolant_iteration(coolant_id)%solution(target_point_index,:) + &
          &grid%points(target_point_id)%coolant_state(coolant_id)%population)
    end if
#endif
  end subroutine apply_population_relaxation

  subroutine deallocate_lvg_workspace
    do coolant_id=1,coolant_count
      call deallocate_coolant_workspace(coolant_work(coolant_id))
    end do
  end subroutine deallocate_lvg_workspace

end subroutine solve_level_populations_lvg

end module level_population_solver_module
