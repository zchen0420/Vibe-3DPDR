module level_population_solver_module
  use definitions, only : dp
  use healpix_types, only : i4b
  use coolants_module, only : COOLANT_COUNT, COOLANT_CII, COOLANT_CI, COOLANT_OI, COOLANT_C12O
  use maincode_module, only : coolant, coolant_iteration, grid, levpop_iteration, maxpoints, nrays, pdr_ptot, runtime, thermal
  use global_module, only : metallicity, NELECT, species_idx
  use convergence_module, only : set_lte_populations
  use level_population_system_module, only : solve_statistical_equilibrium
  use ray_path_module, only : projected_point_id, set_ray_origin_projection

  implicit none

  type coolant_work_item
    integer(kind=i4b) :: coolant_id
    integer(kind=i4b) :: species_abundance_index
    integer(kind=i4b) :: level_count
    integer(kind=i4b) :: temperature_count
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
      endif
#endif
    enddo
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
    enddo
  end subroutine refresh_lte_populations_for_unconverged_points

  subroutine initialize_coolant_work(work_item, coolant_id, species_abundance_index)
    type(coolant_work_item), intent(out) :: work_item
    integer(kind=i4b), intent(in) :: coolant_id
    integer(kind=i4b), intent(in) :: species_abundance_index

    work_item%coolant_id = coolant_id
    work_item%species_abundance_index = species_abundance_index
    work_item%level_count = coolant(coolant_id)%nlevels
    work_item%temperature_count = coolant(coolant_id)%ntemperatures
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

    work_item%evalpop=0.0D0
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
    type(coolant_work_item) :: coolant_work(COOLANT_COUNT)

    call initialize_coolant_work(coolant_work(COOLANT_CII),COOLANT_CII,species_idx%NCx)
    call initialize_coolant_work(coolant_work(COOLANT_CI),COOLANT_CI,species_idx%NC)
    call initialize_coolant_work(coolant_work(COOLANT_OI),COOLANT_OI,species_idx%NO)
    call initialize_coolant_work(coolant_work(COOLANT_C12O),COOLANT_C12O,species_idx%NCO)

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

      call solve_coolant_population(coolant_work(COOLANT_CII), point_index, point_id, .true.)
      call solve_coolant_population(coolant_work(COOLANT_CI), point_index, point_id, .true.)
      call solve_coolant_population(coolant_work(COOLANT_OI), point_index, point_id, .true.)
      call solve_coolant_population(coolant_work(COOLANT_C12O), point_index, point_id, atomic_coolants_converged)

      call deallocate_lvg_workspace
    enddo

    thermal%total_cooling_rate=0.0D0
    do coolant_id = 1, COOLANT_COUNT
      thermal%total_cooling_rate = thermal%total_cooling_rate + coolant_iteration(coolant_id)%cooling_rate
    enddo

  contains

    subroutine allocate_lvg_workspace
      do coolant_id=1,COOLANT_COUNT
        call allocate_coolant_workspace(coolant_work(coolant_id))
      enddo
    end subroutine allocate_lvg_workspace

    subroutine populate_evaluation_populations(source_point_id)
      integer(kind=i4b), intent(in) :: source_point_id
      integer(kind=i4b) :: eval_point_id

      do ray_index=0,nrays-1
        do eval_index=0,grid%points(source_point_id)%epray(ray_index)
          eval_point_id=projected_point_id(grid%points(source_point_id), ray_index, eval_index)

          do coolant_id=1,COOLANT_COUNT
            do level_index=1,coolant_work(coolant_id)%level_count
              coolant_work(coolant_id)%evalpop(ray_index,eval_index,level_index)=&
                  &grid%points(eval_point_id)%coolant_state(coolant_id)%population(level_index)
            enddo
          enddo
        enddo
      enddo
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

      call find_collision_coefficients(work_item%temperature_count,work_item%level_count, &
          & thermal%gas_temperature(target_point_index),coolant(work_item%coolant_id)%temperatures, &
          & coolant(work_item%coolant_id)%h,coolant(work_item%coolant_id)%hp, &
          & coolant(work_item%coolant_id)%el,coolant(work_item%coolant_id)%he, &
          & coolant(work_item%coolant_id)%h2,coolant(work_item%coolant_id)%ph2, &
          & coolant(work_item%coolant_id)%oh2,work_item%collision_coefficients, &
          & abundance_density(target_point_id, species_idx%NH), &
          & abundance_density(target_point_id, species_idx%NPROTON), &
          & abundance_density(target_point_id, NELECT), &
          & abundance_density(target_point_id, species_idx%NHe), &
          & abundance_density(target_point_id, species_idx%NH2), &
          & work_item%coolant_id)
    end subroutine update_collision_rates

    subroutine update_lvg_transition_rates(work_item, target_point_index, target_point_id)
      type(coolant_work_item), intent(inout) :: work_item
      integer(kind=i4b), intent(in) :: target_point_index
      integer(kind=i4b), intent(in) :: target_point_id

      call escape_probability(work_item%transition, thermal%dust_temperature(target_point_index), &
          & nrays, work_item%level_count, coolant(work_item%coolant_id)%a_coeffs, &
          & coolant(work_item%coolant_id)%b_coeffs, work_item%collision_coefficients, &
          & coolant(work_item%coolant_id)%frequencies, work_item%evalpop, maxpoints, &
          & thermal%gas_temperature(target_point_index), runtime%turbulent_velocity, &
          & grid%points(target_point_id)%epray, &
          & grid%points(target_point_id)%coolant_state(work_item%coolant_id)%population, &
          & grid%points(target_point_id)%epoint, coolant(work_item%coolant_id)%weights, &
          & coolant_iteration(work_item%coolant_id)%cooling_rate(target_point_index), &
          & work_item%line, work_item%optical_depth, work_item%coolant_id, &
          & grid%points(target_point_id)%rho, metallicity, work_item%beta)

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
            &0.5D0*(coolant_iteration(coolant_id)%solution(target_point_index,:) + &
            &grid%points(target_point_id)%coolant_state(coolant_id)%population)
      endif
#endif
    end subroutine apply_population_relaxation

    subroutine deallocate_lvg_workspace
      do coolant_id=1,COOLANT_COUNT
        call deallocate_coolant_workspace(coolant_work(coolant_id))
      enddo
    end subroutine deallocate_lvg_workspace

  end subroutine solve_level_populations_lvg

end module level_population_solver_module
