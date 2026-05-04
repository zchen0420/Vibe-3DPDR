module iteration_chemistry_module
  use definitions
  use maincode_module
  use global_module
  use columns_module, only : calc_columndens

  implicit none

contains

  subroutine run_initial_chemistry_iterations
    integer(kind=i4b) :: NRGR,NRH2,NRHD,NRCO,NRCI,NRSI
    real(kind=dp), allocatable :: dummy_abundance(:,:)
    real(kind=dp), allocatable :: dummy_density(:)
    real(kind=dp), allocatable :: dummy_rate(:,:)
    real(kind=dp), allocatable :: dummy_temperature(:)

    allocate(dummy_rate(1:nreac,1:pdr_ptot))
    allocate(dummy_abundance(1:nspec,1:pdr_ptot))
    allocate(dummy_density(1:pdr_ptot))
    allocate(dummy_temperature(1:pdr_ptot))

    do ii=1,CHEMITERATIONS
      write(6,*) 'Chemical iteration ',ii

      do pp=1,pdr_ptot
        p=IDlist_pdr(pp)
        call calculate_reaction_rates(gastemperature(pp),dusttemperature(pp),nrays,pdr(p)%rad_surface(0:nrays-1),&
            &pdr(p)%AV(0:nrays-1),column(pp)%columndens_point(0:nrays-1,1:nspec),&
            &nreac, reactant, product, alpha, beta, gamma, rate, rtmin, rtmax, duplicate, nspec,&
            &NRGR,NRH2,NRHD,NRCO,NRCI,NRSI)
        dummy_rate(:,pp) = rate
        dummy_abundance(:,pp) = pdr(p)%abundance
        dummy_density(pp) = pdr(p)%rho
        dummy_temperature(pp) = gastemperature(pp)
      enddo

      call calculate_abundances(dummy_abundance,dummy_rate,dummy_density,dummy_temperature,pdr_ptot,nspec,nreac)

      do pp=1,pdr_ptot
        p=IDlist_pdr(pp)
        pdr(p)%abundance = dummy_abundance(:,pp)
      enddo

      call calc_columndens
    enddo

    deallocate(dummy_rate)
    deallocate(dummy_abundance)
    deallocate(dummy_density)
    deallocate(dummy_temperature)
  end subroutine run_initial_chemistry_iterations

  subroutine refresh_chemistry_after_temperature_change
    integer(kind=i4b) :: NRGR,NRH2,NRHD,NRCO,NRCI,NRSI
    real(kind=dp), allocatable :: dummy_abundance(:,:)
    real(kind=dp), allocatable :: dummy_density(:)
    real(kind=dp), allocatable :: dummy_rate(:,:)
    real(kind=dp), allocatable :: dummy_temperature(:)

    allocate(dummy_rate(1:nreac,1:pdr_ptot))
    allocate(dummy_abundance(1:nspec,1:pdr_ptot))
    allocate(dummy_density(1:pdr_ptot))
    allocate(dummy_temperature(1:pdr_ptot))

    do ii=1,3
      write(6,*) 'Chemical iteration ',ii

      do pp=1,pdr_ptot
        p=IDlist_pdr(pp)
        call calculate_reaction_rates(gastemperature(pp),dusttemperature(pp),nrays,pdr(p)%rad_surface(0:nrays-1),&
            &pdr(p)%AV(0:nrays-1),column(pp)%columndens_point(0:nrays-1,1:nspec),&
            &nreac, reactant, product, alpha, beta, gamma, rate, rtmin, rtmax, duplicate, nspec,&
            &NRGR,NRH2,NRHD,NRCO,NRCI,NRSI)
        dummy_rate(:,pp) = rate
        dummy_abundance(:,pp) = pdr(p)%abundance
        dummy_density(pp) = pdr(p)%rho
        dummy_temperature(pp) = gastemperature(pp)
      enddo

      call calculate_abundances(dummy_abundance,dummy_rate,dummy_density,dummy_temperature,pdr_ptot,nspec,nreac)

      do pp=1,pdr_ptot
#ifdef THERMALBALANCE
        if (converged(pp)) cycle
#endif
        p=IDlist_pdr(pp)
        pdr(p)%abundance = dummy_abundance(:,pp)
      enddo

      call calc_columndens
    enddo

    deallocate(dummy_rate)
    deallocate(dummy_abundance)
    deallocate(dummy_density)
    deallocate(dummy_temperature)
  end subroutine refresh_chemistry_after_temperature_change

end module iteration_chemistry_module
