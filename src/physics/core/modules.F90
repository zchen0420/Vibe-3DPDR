module maincode_module
  use iso_c_binding
  use definitions
  use healpix_types
  use chemistry_network_module, only : chemistry_state
  use geometry_state_module, only : geometry_state
  use runtime_state_module, only : runtime_config_state
  use simulation_grid_module, only : pdr_node, simulation_grid
  use thermal_state_module, only : thermal_state
  use coolants_module, only : coolant_data, point_coolant_state, coolant_iteration_state

  integer, parameter :: x_axis = 1
  integer, parameter :: y_axis = 2
  integer, parameter :: z_axis = 3

  integer(kind=i4b) :: level ! current level
  integer(kind=i4b) :: nrays ! no. of rays on current level
  integer(kind=i4b) :: nside ! refer to healpix manual
  integer(kind=i4b) :: iteration
  integer(kind=i4b) :: nspec, nreac
  integer(kind=i4b) :: maxpoints

  type(runtime_config_state) :: runtime
  real(kind=dp) :: theta_crit ! critical theta angle to produce evaluation point
  real(kind=dp) :: tlow0
  real(kind=dp) :: thigh0
  real(kind=dp) :: tmin
  real(kind=dp) :: tmax
  real(kind=dp) :: fcrit
  real(kind=dp) :: tdiff
  real(kind=dp) :: dust_temperature
  !  real(kind=dp) :: h

  character(len=3) :: fieldchoice
  real(kind=dp) :: gext(1)

  type(geometry_state) :: geometry
  type(coolant_data), allocatable :: coolant(:)

  type(chemistry_state) :: chemistry

  real(kind=dp),bind(c,name='maincode_module_mp_start_time_'):: start_time
  real(kind=dp),bind(c,name='maincode_module_mp_end_time_'):: end_time

  type(thermal_state) :: thermal

  type(simulation_grid) :: grid
  type(coolant_iteration_state), allocatable :: coolant_iteration(:)

  integer(kind=i4b)::levpop_iteration
  real(kind=dp) :: rho_min
  real(kind=dp) :: rho_max
#ifdef GUESS_TEMP
  real(kind=dp),allocatable :: tmin_array(:)
  real(kind=dp),allocatable :: tmax_array(:)
#endif

  logical :: level_conv,first_time
  integer(kind=i4b) :: grand_ptot
  integer(kind=i4b) :: pdr_ptot,ion_ptot,dark_ptot
end module maincode_module
