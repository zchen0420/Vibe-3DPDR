module point_reaction_rates_module
  use healpix_types, only : i4b
  use maincode_module, only : chemistry, geometry, grid, thermal
  use reaction_rates_module, only : reaction_rate_environment, reaction_rate_indices, calculate_reaction_rates
  implicit none

contains

  subroutine calculate_point_reaction_rates(point_index, point_id, indices)
    integer(kind=i4b), intent(in) :: point_index
    integer(kind=i4b), intent(in) :: point_id
    type(reaction_rate_indices), intent(out), optional :: indices
    type(reaction_rate_environment) :: environment
    type(reaction_rate_indices) :: local_indices

    environment%gas_temperature = thermal%gas_temperature(point_index)
    environment%dust_temperature = thermal%dust_temperature(point_index)
    environment%radiation_surface => grid%points(point_id)%rad_surface
    environment%visual_extinction => grid%points(point_id)%av
    environment%column_density => geometry%column_density(point_index)%columndens_point

    call calculate_reaction_rates(environment, chemistry%network, chemistry%rate, local_indices)

    if (present(indices)) indices = local_indices
  end subroutine calculate_point_reaction_rates

end module point_reaction_rates_module
