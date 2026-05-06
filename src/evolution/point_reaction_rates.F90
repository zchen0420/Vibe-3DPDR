module point_reaction_rates_module
  use definitions, only : dp
  use healpix_types, only : i4b
  use maincode_module, only : chemistry, geometry, grid, nreac, nrays, nspec, thermal
  implicit none

  type :: reaction_rate_indices
    integer(kind=i4b) :: grain_surface
    integer(kind=i4b) :: h2_photodissociation
    integer(kind=i4b) :: hd_photodissociation
    integer(kind=i4b) :: co_photodissociation
    integer(kind=i4b) :: carbon_photoionization
    integer(kind=i4b) :: silicon_photoionization
  end type reaction_rate_indices

contains

  subroutine calculate_point_reaction_rates(point_index, point_id, indices)
    integer(kind=i4b), intent(in) :: point_index
    integer(kind=i4b), intent(in) :: point_id
    type(reaction_rate_indices), intent(out), optional :: indices
    type(reaction_rate_indices) :: local_indices

    call calculate_reaction_rates(thermal%gas_temperature(point_index),thermal%dust_temperature(point_index), &
        &nrays,grid%points(point_id)%rad_surface(0:nrays-1),grid%points(point_id)%AV(0:nrays-1), &
        &geometry%column_density(point_index)%columndens_point(0:nrays-1,1:nspec), &
        &nreac, chemistry%network%reactant, chemistry%network%product, chemistry%network%alpha, &
        &chemistry%network%beta, chemistry%network%gamma, chemistry%rate, chemistry%network%rtmin, &
        &chemistry%network%rtmax, chemistry%network%duplicate, nspec, &
        &local_indices%grain_surface, local_indices%h2_photodissociation, &
        &local_indices%hd_photodissociation, local_indices%co_photodissociation, &
        &local_indices%carbon_photoionization, local_indices%silicon_photoionization)

    if (present(indices)) indices = local_indices
  end subroutine calculate_point_reaction_rates

end module point_reaction_rates_module
