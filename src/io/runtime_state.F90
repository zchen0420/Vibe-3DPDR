module runtime_state_module
  use definitions, only : dp
  use healpix_types, only : i4b
  implicit none

  type :: runtime_config_state
    character(len=128), pointer :: input_file => null()
    character(len=128), pointer :: output_prefix => null()
    integer(kind=i4b), pointer :: total_iterations => null()
    integer(kind=i4b), pointer :: output_interval => null()
    integer, pointer :: chemistry_iterations => null()
    real(kind=dp), pointer :: temperature_guess => null()
    real(kind=dp), pointer :: turbulent_velocity => null()
    real(kind=dp), pointer :: input_turbulent_velocity => null()
    real(kind=dp), pointer :: cosmic_ray_ionization_rate => null()
    real(kind=dp), pointer :: av_scale => null()
    real(kind=dp), pointer :: uv_scale => null()
  end type runtime_config_state

contains

  subroutine initialize_runtime_config_state(config)
    type(runtime_config_state), intent(inout) :: config

    allocate(config%input_file)
    allocate(config%output_prefix)
    allocate(config%total_iterations)
    allocate(config%output_interval)
    allocate(config%chemistry_iterations)
    allocate(config%temperature_guess)
    allocate(config%turbulent_velocity)
    allocate(config%input_turbulent_velocity)
    allocate(config%cosmic_ray_ionization_rate)
    allocate(config%av_scale)
    allocate(config%uv_scale)
  end subroutine initialize_runtime_config_state

end module runtime_state_module
