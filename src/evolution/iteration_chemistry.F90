module iteration_chemistry_module
  use definitions, only : dp
  use healpix_types, only : i4b
  use maincode_module, only : chemistry, grid, nreac, nspec, pdr_ptot, runtime, thermal
  use columns_module, only : calc_columndens
  use point_reaction_rates_module, only : calculate_point_reaction_rates

  implicit none

contains

  subroutine run_initial_chemistry_iterations
    integer(kind=i4b) :: chemistry_iteration, point_id, point_index
    real(kind=dp), allocatable :: dummy_abundance(:,:)
    real(kind=dp), allocatable :: dummy_density(:)
    real(kind=dp), allocatable :: dummy_rate(:,:)
    real(kind=dp), allocatable :: dummy_temperature(:)

    allocate(dummy_rate(1:nreac,1:pdr_ptot))
    allocate(dummy_abundance(1:nspec,1:pdr_ptot))
    allocate(dummy_density(1:pdr_ptot))
    allocate(dummy_temperature(1:pdr_ptot))

    do chemistry_iteration=1,runtime%chemistry_iterations
      write(6,*) 'Chemical iteration ',chemistry_iteration

      do point_index=1,pdr_ptot
        point_id=grid%pdr_ids(point_index)
        call calculate_point_reaction_rates(point_index, point_id)
        dummy_rate(:,point_index) = chemistry%rate
        dummy_abundance(:,point_index) = grid%points(point_id)%abundance
        dummy_density(point_index) = grid%points(point_id)%rho
        dummy_temperature(point_index) = thermal%gas_temperature(point_index)
      end do

      call calculate_abundances(dummy_abundance,dummy_rate,dummy_density,dummy_temperature,pdr_ptot,nspec,nreac)

      do point_index=1,pdr_ptot
        point_id=grid%pdr_ids(point_index)
        grid%points(point_id)%abundance = dummy_abundance(:,point_index)
      end do

      call calc_columndens(.false.)
    end do

    deallocate(dummy_rate)
    deallocate(dummy_abundance)
    deallocate(dummy_density)
    deallocate(dummy_temperature)
  end subroutine run_initial_chemistry_iterations

  subroutine refresh_chemistry_after_temperature_change
    integer(kind=i4b) :: chemistry_iteration, point_id, point_index
    real(kind=dp), allocatable :: dummy_abundance(:,:)
    real(kind=dp), allocatable :: dummy_density(:)
    real(kind=dp), allocatable :: dummy_rate(:,:)
    real(kind=dp), allocatable :: dummy_temperature(:)

    allocate(dummy_rate(1:nreac,1:pdr_ptot))
    allocate(dummy_abundance(1:nspec,1:pdr_ptot))
    allocate(dummy_density(1:pdr_ptot))
    allocate(dummy_temperature(1:pdr_ptot))

    do chemistry_iteration=1,3
      write(6,*) 'Chemical iteration ',chemistry_iteration

      do point_index=1,pdr_ptot
        point_id=grid%pdr_ids(point_index)
        call calculate_point_reaction_rates(point_index, point_id)
        dummy_rate(:,point_index) = chemistry%rate
        dummy_abundance(:,point_index) = grid%points(point_id)%abundance
        dummy_density(point_index) = grid%points(point_id)%rho
        dummy_temperature(point_index) = thermal%gas_temperature(point_index)
      end do

      call calculate_abundances(dummy_abundance,dummy_rate,dummy_density,dummy_temperature,pdr_ptot,nspec,nreac)

      do point_index=1,pdr_ptot
#ifdef THERMALBALANCE
        if (thermal%thermal_converged(point_index)) cycle
#endif
        point_id=grid%pdr_ids(point_index)
        grid%points(point_id)%abundance = dummy_abundance(:,point_index)
      end do

      call calc_columndens(.false.)
    end do

    deallocate(dummy_rate)
    deallocate(dummy_abundance)
    deallocate(dummy_density)
    deallocate(dummy_temperature)
  end subroutine refresh_chemistry_after_temperature_change

end module iteration_chemistry_module
