module iteration_convergence_module
  use definitions
  use maincode_module
  use global_module
  use convergence_module, only : population_has_converged, print_convergence_count

  implicit none

contains

  subroutine initialise_level_convergence_state
    allocate(level_converged(0:pdr_ptot))
    level_converged=.false.
    allocate(CII_conv(0:pdr_ptot))
    allocate(CI_conv(0:pdr_ptot))
    allocate(OI_conv(0:pdr_ptot))
    allocate(C12O_conv(0:pdr_ptot))
    CII_conv=.false.
    CI_conv=.false.
    OI_conv=.false.
    C12O_conv=.false.
    allocate(CII_RELCH(0:pdr_ptot,1:CII_NLEV))
    allocate(CI_RELCH(0:pdr_ptot,1:CI_NLEV))
    allocate(OI_RELCH(0:pdr_ptot,1:OI_NLEV))
    allocate(C12O_RELCH(0:pdr_ptot,1:C12O_NLEV))
  end subroutine initialise_level_convergence_state

  subroutine evaluate_iteration_convergence(stop_iterations)
    logical, intent(out) :: stop_iterations
    integer(kind=i4b) :: converged_count
    integer(kind=i4b) :: point_index
    integer(kind=i4b) :: point_id

    stop_iterations=.false.

    write(6,*) 'Checking for convergence...'

#ifdef THERMALBALANCE
    if (level_conv.and.first_time) first_time=.false.

    converged_count=0
    do point_index=1,pdr_ptot
      if (converged(point_index)) converged_count=converged_count+1
    enddo

    if (level_conv) then
      write(6,*) 'Resetting [level_conv=.false.]'
      level_conv=.false.
    endif

    thermal_percentage = 100.D0*real(converged_count,kind=dp)/real(pdr_ptot,kind=dp)
    write(*,'(" Thermal balance is ",F5.1,"% converged.")') thermal_percentage
    write(*,'(" [",I6,"/",I6,"]")') converged_count,pdr_ptot

    if (converged_count.eq.pdr_ptot) then
      write(6,*) '#### Converged through thermal balance ####'
      stop_iterations=.true.
      return
    endif
#endif

    RELCH_conv=.true.
    level_converged=.false.

    do point_index=1,pdr_ptot
      point_id=IDlist_pdr(point_index)

      CII_conv(point_index)=.true.
      CI_conv(point_index)=.true.
      OI_conv(point_index)=.true.
      C12O_conv(point_index)=.true.

      CII_conv(point_index) = population_has_converged(CII_SOLUTION(point_index,1:CII_NLEV), &
          &pdr(point_id)%CII_pop, pdr(point_id)%abundance(NCx), CII_RELCH(point_index,1:CII_NLEV), CII_NLEV)
      CI_conv(point_index) = population_has_converged(CI_SOLUTION(point_index,1:CI_NLEV), &
          &pdr(point_id)%CI_pop, pdr(point_id)%abundance(NC), CI_RELCH(point_index,1:CI_NLEV), CI_NLEV)
      OI_conv(point_index) = population_has_converged(OI_SOLUTION(point_index,1:OI_NLEV), &
          &pdr(point_id)%OI_pop, pdr(point_id)%abundance(NO), OI_RELCH(point_index,1:OI_NLEV), OI_NLEV)
      C12O_conv(point_index) = population_has_converged(C12O_SOLUTION(point_index,1:C12O_NLEV), &
          &pdr(point_id)%C12O_pop, pdr(point_id)%abundance(NCO), C12O_RELCH(point_index,1:C12O_NLEV), &
          &C12O_NLEV)

      if (.not.(CII_conv(point_index).and.CI_conv(point_index).and. &
          &OI_conv(point_index).and.C12O_conv(point_index))) then
        RELCH_conv=.false.
      endif

      if (CII_conv(point_index).and.CI_conv(point_index).and. &
          &OI_conv(point_index).and.C12O_conv(point_index)) then
        level_converged(point_index) = .true.
      endif
    enddo

    if (.not.relch_conv) return

    write(6,*) '#### Converged through level populations ####'
    levpop_iteration=0

#ifdef THERMALBALANCE
    write(6,*) 'Enabling thermal balance routine in next iteration'
    level_conv=.true.
#else
    stop_iterations=.true.
#endif
  end subroutine evaluate_iteration_convergence

  subroutine report_iteration_convergence
    integer(kind=i4b) :: level_count
    integer(kind=i4b) :: point_index

    level_count=0
    CII_i=0
    CI_i=0
    OI_i=0
    C12O_i=0

    do point_index=1,pdr_ptot
      if (level_converged(point_index)) level_count=level_count+1
      if (CII_conv(point_index)) CII_i=CII_i+1
      if (CI_conv(point_index)) CI_i=CI_i+1
      if (OI_conv(point_index)) OI_i=OI_i+1
      if (C12O_conv(point_index)) C12O_i=C12O_i+1
    enddo

    levpop_percentage  = 100.D0*real(level_count,kind=dp)/real(pdr_ptot,kind=dp)
    CII_percentage = int(100.D0*real(CII_i,kind=dp)/real(pdr_ptot,kind=dp),kind=i4b)
    CI_percentage = int(100.D0*real(CI_i,kind=dp)/real(pdr_ptot,kind=dp),kind=i4b)
    OI_percentage = int(100.D0*real(OI_i,kind=dp)/real(pdr_ptot,kind=dp),kind=i4b)
    C12O_percentage = int(100.D0*real(C12O_i,kind=dp)/real(pdr_ptot,kind=dp),kind=i4b)

    call print_convergence_count('CII', CII_percentage, CII_i, pdr_ptot)
    call print_convergence_count('CI', CI_percentage, CI_i, pdr_ptot)
    call print_convergence_count('OI', OI_percentage, OI_i, pdr_ptot)
    call print_convergence_count('CO', C12O_percentage, C12O_i, pdr_ptot)
    call print_convergence_count('Level populations', levpop_percentage, level_count, pdr_ptot)

#ifdef THERMALBALANCE
    if (int(levpop_percentage,kind=i4b).ge.100) then
      level_converged = .false.
      write(6,*) 'Resetting [level_converged=.false.] array'
    endif
#endif
  end subroutine report_iteration_convergence

  subroutine update_population_densities
    integer(kind=i4b) :: point_index
    integer(kind=i4b) :: point_id
    integer(kind=i4b) :: level_index

    write(6,*) 'Updating population densities...'

    do point_index=1,pdr_ptot
      point_id=IDlist_pdr(point_index)

      do level_index=1,CII_NLEV
        pdr(point_id)%CII_pop(level_index) = CII_solution(point_index,level_index)
      enddo
      do level_index=1,CI_NLEV
        pdr(point_id)%CI_pop(level_index) = CI_solution(point_index,level_index)
      enddo
      do level_index=1,OI_NLEV
        pdr(point_id)%OI_pop(level_index) = OI_solution(point_index,level_index)
      enddo
      do level_index=1,C12O_NLEV
        pdr(point_id)%C12O_pop(level_index) = C12O_solution(point_index,level_index)
      enddo
    enddo
  end subroutine update_population_densities

  subroutine write_iteration_timing
    call cpu_time(t3)
    write(6,*) 'Iteration time = ',t3-t3b,' seconds.'
    write(6,*) 'Total time = ',t3,' seconds.'
  end subroutine write_iteration_timing

end module iteration_convergence_module
