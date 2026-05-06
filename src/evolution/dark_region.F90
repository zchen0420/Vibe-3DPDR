module dark_region_module
  use definitions, only : dp
  use healpix_types, only : i4b
  use point_reaction_rates_module, only : calculate_point_reaction_rates
  implicit none

contains

  subroutine dark_molecular_region
    use maincode_module, only : chemistry, dark_ptot, grid, nreac, nspec, rho_max, runtime, thermal
    use coolants_module, only : coolant_count
    use columns_module, only : calc_columndens
    use convergence_module, only : set_lte_populations

    integer(kind=i4b) :: chemistry_iteration
    integer(kind=i4b) :: representative_point_id
    integer(kind=i4b) :: dark_point_index
    integer(kind=i4b) :: dark_point_id
    integer(kind=i4b) :: coolant_id
    real(kind=dp) :: molecular_density
    real(kind=dp), allocatable :: local_rate(:,:)
    real(kind=dp), allocatable :: local_abundance(:,:)
    real(kind=dp), allocatable :: local_density(:)
    real(kind=dp), allocatable :: local_temperature(:)

    write(6,*) ''
    write(6,*) '*** Dark Molecular Region ***'
    write(6,*) 'Calculating LTE level populations...'
    molecular_density = 2.0d0*rho_max
    write(6,*) 'Density = ',molecular_density

    representative_point_id=grid%dark_ids(1)

    do chemistry_iteration=1,runtime%chemistry_iterations
      write(6,*) 'Chemical iteration ',chemistry_iteration

      allocate(local_rate(1:nreac,1))
      call calculate_point_reaction_rates(0, representative_point_id)
      local_rate(:,1) = chemistry%rate

      allocate(local_abundance(1:nspec,1))
      allocate(local_density(1))
      allocate(local_temperature(1))

      local_abundance(:,1) = grid%points(representative_point_id)%abundance
      local_density(1) = grid%points(representative_point_id)%rho
      local_temperature(1) = thermal%gas_temperature(0)

      call calculate_abundances(local_abundance,local_rate,local_density,local_temperature,1,nspec,nreac)
      grid%points(representative_point_id)%abundance = local_abundance(:,1)

      deallocate(local_abundance)
      deallocate(local_density)
      deallocate(local_temperature)
      deallocate(local_rate)

      call calc_columndens(.false.)
    end do

    call set_lte_populations(representative_point_id, thermal%gas_temperature(0), molecular_density)

    write(6,*) 'Assigning properties to all dark particles'
    do dark_point_index=2,dark_ptot
      dark_point_id=grid%dark_ids(dark_point_index)
      do coolant_id=1,coolant_count
        grid%points(dark_point_id)%coolant_state(coolant_id)%population = &
            &grid%points(representative_point_id)%coolant_state(coolant_id)%population
      end do
      grid%points(dark_point_id)%abundance = grid%points(representative_point_id)%abundance
    end do
    write(6,*) 'Done! Proceeding with grid%points...'
  end subroutine dark_molecular_region

end module dark_region_module
