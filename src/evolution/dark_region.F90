module dark_region_module
  use definitions, only : dp
  use healpix_types, only : i4b
  implicit none

contains

  subroutine dark_molecular_region
    use maincode_module
    use global_module
    use columns_module, only : calc_columndens
    use convergence_module, only : set_lte_populations

    integer(kind=i4b) :: chemistry_iteration
    integer(kind=i4b) :: representative_point_id
    integer(kind=i4b) :: dark_point_index
    integer(kind=i4b) :: dark_point_id
    integer(kind=i4b) :: grain_rate_index
    integer(kind=i4b) :: h2_rate_index
    integer(kind=i4b) :: hd_rate_index
    integer(kind=i4b) :: co_rate_index
    integer(kind=i4b) :: ci_rate_index
    integer(kind=i4b) :: si_rate_index
    real(kind=dp), allocatable :: local_rate(:,:)
    real(kind=dp), allocatable :: local_abundance(:,:)
    real(kind=dp), allocatable :: local_density(:)
    real(kind=dp), allocatable :: local_temperature(:)

    write(6,*) ''
    write(6,*) '*** Dark Molecular Region ***'
    write(6,*) 'Calculating LTE level populations...'
    n_H = 2.0D0*rho_max
    write(6,*) 'Density = ',n_H

    representative_point_id=IDlist_dark(1)

    do chemistry_iteration=1,CHEMITERATIONS
      write(6,*) 'Chemical iteration ',chemistry_iteration

      allocate(local_rate(1:nreac,1))
      call CALCULATE_REACTION_RATES(gastemperature(0),dusttemperature(0),nrays,&
          &pdr(representative_point_id)%rad_surface(0:nrays-1),&
          &pdr(representative_point_id)%AV(0:nrays-1),column(0)%columndens_point(0:nrays-1,1:nspec),&
          &nreac, reactant, product, alpha, beta, gamma, rate, rtmin, rtmax, duplicate, nspec,&
          &grain_rate_index,h2_rate_index,hd_rate_index,co_rate_index,ci_rate_index,si_rate_index)
      local_rate(:,1) = rate

      allocate(local_abundance(1:nspec,1))
      allocate(local_density(1))
      allocate(local_temperature(1))

      local_abundance(:,1) = pdr(representative_point_id)%abundance
      local_density(1) = pdr(representative_point_id)%rho
      local_temperature(1) = gastemperature(0)

      call CALCULATE_ABUNDANCES(local_abundance,local_rate,local_density,local_temperature,1,NSPEC,NREAC)
      pdr(representative_point_id)%abundance = local_abundance(:,1)

      deallocate(local_abundance)
      deallocate(local_density)
      deallocate(local_temperature)
      deallocate(local_rate)

      call calc_columndens
    enddo

    call set_lte_populations(representative_point_id, gastemperature(0), n_H)

    write(6,*) 'Assigning properties to all dark particles'
    do dark_point_index=2,dark_ptot
      dark_point_id=IDlist_dark(dark_point_index)
      pdr(dark_point_id)%cii_pop = pdr(representative_point_id)%cii_pop
      pdr(dark_point_id)%ci_pop = pdr(representative_point_id)%ci_pop
      pdr(dark_point_id)%oi_pop = pdr(representative_point_id)%oi_pop
      pdr(dark_point_id)%c12o_pop = pdr(representative_point_id)%c12o_pop
      pdr(dark_point_id)%abundance = pdr(representative_point_id)%abundance
    enddo
    write(6,*) 'Done! Proceeding with PDR...'
  end subroutine dark_molecular_region

end module dark_region_module
