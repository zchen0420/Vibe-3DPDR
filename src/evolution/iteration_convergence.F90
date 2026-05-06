module iteration_convergence_module
  use definitions, only : dp
  use healpix_types, only : i4b
  use coolants_module, only : COOLANT_COUNT, COOLANT_CII, COOLANT_CI, COOLANT_OI, COOLANT_C12O, coolant_label
  use thermal_state_module, only : allocate_level_convergence_state
  use maincode_module, only : coolant, coolant_iteration, first_time, grid, level_conv, levpop_iteration, pdr_ptot, thermal
  use global_module, only : species_idx
  use convergence_module, only : population_has_converged, print_convergence_count

  implicit none

contains

  subroutine initialise_level_convergence_state
    integer(kind=i4b) :: coolant_id

    if (.not.allocated(coolant_iteration)) allocate(coolant_iteration(1:COOLANT_COUNT))
    call allocate_level_convergence_state(thermal, pdr_ptot)

    do coolant_id = 1, COOLANT_COUNT
      allocate(coolant_iteration(coolant_id)%converged(0:pdr_ptot))
      coolant_iteration(coolant_id)%converged = .false.
      allocate(coolant_iteration(coolant_id)%relative_change(0:pdr_ptot,1:coolant(coolant_id)%nlevels))
    enddo
  end subroutine initialise_level_convergence_state

  subroutine evaluate_iteration_convergence(stop_iterations)
    logical, intent(out) :: stop_iterations
    integer(kind=i4b) :: converged_count
    integer(kind=i4b) :: coolant_id
    integer(kind=i4b) :: point_index
    integer(kind=i4b) :: point_id
    logical :: populations_converged
#ifdef THERMALBALANCE
    real(kind=dp) :: thermal_percentage
#endif

    stop_iterations=.false.

    write(6,*) 'Checking for convergence...'

#ifdef THERMALBALANCE
    if (level_conv.and.first_time) first_time=.false.

    converged_count=0
    do point_index=1,pdr_ptot
      if (thermal%thermal_converged(point_index)) converged_count=converged_count+1
    enddo

    if (level_conv) then
      write(6,*) 'Resetting [level_conv=.false.]'
      level_conv=.false.
    endif

    thermal_percentage = 100.D0*real(converged_count,kind=dp)/real(pdr_ptot,kind=dp)
    write(*,'(" Thermal balance is ",F5.1,"% thermal%thermal_converged.")') thermal_percentage
    write(*,'(" [",I6,"/",I6,"]")') converged_count,pdr_ptot

    if (converged_count.eq.pdr_ptot) then
      write(6,*) '#### thermal%thermal_converged through thermal balance ####'
      stop_iterations=.true.
      return
    endif
#endif

    populations_converged=.true.
    thermal%level_population_converged=.false.

    do point_index=1,pdr_ptot
      point_id=grid%pdr_ids(point_index)

      do coolant_id = 1, COOLANT_COUNT
        coolant_iteration(coolant_id)%converged(point_index)=.true.
      enddo

      call update_coolant_convergence(COOLANT_CII, point_index, point_id, species_idx%NCx)
      call update_coolant_convergence(COOLANT_CI, point_index, point_id, species_idx%NC)
      call update_coolant_convergence(COOLANT_OI, point_index, point_id, species_idx%NO)
      call update_coolant_convergence(COOLANT_C12O, point_index, point_id, species_idx%NCO)

      if (.not.point_coolants_converged(point_index)) then
        populations_converged=.false.
      endif

      if (point_coolants_converged(point_index)) then
        thermal%level_population_converged(point_index) = .true.
      endif
    enddo

    if (.not.populations_converged) return

    write(6,*) '#### thermal%thermal_converged through level populations ####'
    levpop_iteration=0

#ifdef THERMALBALANCE
    write(6,*) 'Enabling thermal balance routine in next iteration'
    level_conv=.true.
#else
    stop_iterations=.true.
#endif
  end subroutine evaluate_iteration_convergence

  subroutine report_iteration_convergence(atomic_coolants_converged)
    logical, intent(out) :: atomic_coolants_converged
    integer(kind=i4b) :: coolant_id
    integer(kind=i4b) :: level_count
    integer(kind=i4b) :: point_index
    integer(kind=i4b) :: coolant_converged_count(1:COOLANT_COUNT)
    real(kind=dp) :: level_population_percentage
    real(kind=dp) :: coolant_percentage(1:COOLANT_COUNT)

    level_count=0
    coolant_converged_count=0

    do point_index=1,pdr_ptot
      if (thermal%level_population_converged(point_index)) level_count=level_count+1
      do coolant_id = 1, COOLANT_COUNT
        if (coolant_iteration(coolant_id)%converged(point_index)) then
          coolant_converged_count(coolant_id)=coolant_converged_count(coolant_id)+1
        endif
      enddo
    enddo

    level_population_percentage = 100.D0*real(level_count,kind=dp)/real(pdr_ptot,kind=dp)
    do coolant_id = 1, COOLANT_COUNT
      coolant_percentage(coolant_id) = &
          &int(100.D0*real(coolant_converged_count(coolant_id),kind=dp)/real(pdr_ptot,kind=dp),kind=i4b)
    enddo
    atomic_coolants_converged = coolant_converged_count(COOLANT_CII).eq.pdr_ptot.and. &
        &coolant_converged_count(COOLANT_CI).eq.pdr_ptot.and.coolant_converged_count(COOLANT_OI).eq.pdr_ptot

    do coolant_id = 1, COOLANT_COUNT
      call print_convergence_count(coolant_label(coolant_id), coolant_percentage(coolant_id), &
          &coolant_converged_count(coolant_id), pdr_ptot)
    enddo
    call print_convergence_count('Level populations', level_population_percentage, level_count, pdr_ptot)

#ifdef THERMALBALANCE
    if (int(level_population_percentage,kind=i4b).ge.100) then
      thermal%level_population_converged = .false.
      write(6,*) 'Resetting [thermal%level_population_converged=.false.] array'
    endif
#endif
  end subroutine report_iteration_convergence

  subroutine update_population_densities
    integer(kind=i4b) :: point_index
    integer(kind=i4b) :: point_id
    integer(kind=i4b) :: coolant_id
    integer(kind=i4b) :: level_index

    write(6,*) 'Updating population densities...'

    do point_index=1,pdr_ptot
      point_id=grid%pdr_ids(point_index)

      do coolant_id = 1, COOLANT_COUNT
        do level_index=1,coolant(coolant_id)%nlevels
          grid%points(point_id)%coolant_state(coolant_id)%population(level_index) = &
              &coolant_iteration(coolant_id)%solution(point_index,level_index)
        enddo
      enddo
    enddo
  end subroutine update_population_densities

  subroutine update_coolant_convergence(coolant_id, point_index, point_id, abundance_index)
    integer(kind=i4b), intent(in) :: coolant_id
    integer(kind=i4b), intent(in) :: point_index
    integer(kind=i4b), intent(in) :: point_id
    integer(kind=i4b), intent(in) :: abundance_index
    integer(kind=i4b) :: level_count

    level_count = coolant(coolant_id)%nlevels
    coolant_iteration(coolant_id)%converged(point_index) = &
        &population_has_converged(coolant_iteration(coolant_id)%solution(point_index,1:level_count), &
        &grid%points(point_id)%coolant_state(coolant_id)%population, grid%points(point_id)%abundance(abundance_index), &
        &coolant_iteration(coolant_id)%relative_change(point_index,1:level_count), level_count)
  end subroutine update_coolant_convergence

  logical function point_coolants_converged(point_index)
    integer(kind=i4b), intent(in) :: point_index
    integer(kind=i4b) :: coolant_id

    point_coolants_converged = .true.
    do coolant_id = 1, COOLANT_COUNT
      if (.not.coolant_iteration(coolant_id)%converged(point_index)) then
        point_coolants_converged = .false.
        return
      endif
    enddo
  end function point_coolants_converged

  subroutine write_iteration_timing(iteration_start_time)
    real, intent(in) :: iteration_start_time
    real :: iteration_end_time

    call cpu_time(iteration_end_time)
    write(6,*) 'Iteration time = ',iteration_end_time-iteration_start_time,' seconds.'
    write(6,*) 'Total time = ',iteration_end_time,' seconds.'
  end subroutine write_iteration_timing

end module iteration_convergence_module
