module thermal_state_module
  use definitions, only : dp
  use healpix_types, only : i4b
  implicit none

  type :: thermal_state
    real(kind=dp), pointer :: dust_temperature(:) => null()
    real(kind=dp), pointer :: gas_temperature(:) => null()
    real(kind=dp), pointer :: previous_gas_temperature(:) => null()
    real(kind=dp), pointer :: total_cooling_rate(:) => null()
    real(kind=dp), pointer :: mean_balance(:) => null()
    real(kind=dp), pointer :: relative_balance(:) => null()
    real(kind=dp), pointer :: low_temperature(:) => null()
    real(kind=dp), pointer :: high_temperature(:) => null()
    logical, pointer :: thermal_converged(:) => null()
    logical, pointer :: force_level_minimum(:) => null()
    logical, pointer :: level_population_converged(:) => null()
  end type thermal_state

contains

  subroutine allocate_temperature_state(state, pdr_count)
    type(thermal_state), intent(inout) :: state
    integer(kind=i4b), intent(in) :: pdr_count

    allocate(state%dust_temperature(0:pdr_count))
    allocate(state%gas_temperature(0:pdr_count))
    allocate(state%previous_gas_temperature(0:pdr_count))
    allocate(state%total_cooling_rate(0:pdr_count))
  end subroutine allocate_temperature_state

  subroutine allocate_thermal_balance_state(state, pdr_count)
    type(thermal_state), intent(inout) :: state
    integer(kind=i4b), intent(in) :: pdr_count

    allocate(state%relative_balance(0:pdr_count))
    allocate(state%mean_balance(0:pdr_count))
    allocate(state%low_temperature(0:pdr_count))
    allocate(state%high_temperature(0:pdr_count))
    state%mean_balance = 0.0D0
  end subroutine allocate_thermal_balance_state

  subroutine allocate_thermal_convergence_state(state, pdr_count)
    type(thermal_state), intent(inout) :: state
    integer(kind=i4b), intent(in) :: pdr_count

    allocate(state%thermal_converged(0:pdr_count))
    allocate(state%force_level_minimum(0:pdr_count))
    state%thermal_converged = .false.
    state%force_level_minimum = .false.
  end subroutine allocate_thermal_convergence_state

  subroutine allocate_level_convergence_state(state, pdr_count)
    type(thermal_state), intent(inout) :: state
    integer(kind=i4b), intent(in) :: pdr_count

    allocate(state%level_population_converged(0:pdr_count))
    state%level_population_converged = .false.
  end subroutine allocate_level_convergence_state

end module thermal_state_module
