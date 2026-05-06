MODULE maincode_module
  USE ISO_C_BINDING
  use definitions
  use healpix_types
  use chemistry_network_module, only : chemistry_state
  use geometry_state_module, only : geometry_state
  use runtime_state_module, only : runtime_config_state
  use simulation_grid_module, only : pdr_node, simulation_grid
  use thermal_state_module, only : thermal_state
  use coolants_module, only : coolant_data, point_coolant_state, coolant_iteration_state

  integer, parameter :: X_AXIS = 1
  integer, parameter :: Y_AXIS = 2
  integer, parameter :: Z_AXIS = 3

  integer(kind=I4B) :: level        ! current level
  integer(kind=I4B) :: nrays        ! no. of rays on current level
  integer(kind=I4B) :: nside        ! refer to healpix manual
  integer(kind=i4b) :: iteration
  integer(kind=i4b) :: NSPEC, NREAC
  integer(kind=i4b) :: maxpoints

  type(runtime_config_state) :: runtime
  real(kind=dp) :: theta_crit  ! critical theta angle to produce evaluation point
  real(kind=dp) :: Tlow0
  real(kind=dp) :: Thigh0
  real(kind=dp) :: Tmin
  real(kind=dp) :: Tmax
  real(kind=dp) :: Fcrit
  real(kind=dp) :: Tdiff
  real(kind=dp) :: dust_temperature
  !  real(kind=dp) :: h

  character(len=3) :: fieldchoice
  real(kind=dp) :: Gext(1)

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
  real(kind=dp),allocatable :: Tmin_array(:)
  real(kind=dp),allocatable :: Tmax_array(:)
#endif

  logical :: level_conv,first_time
  integer(kind=i4b) :: grand_ptot
  integer(kind=i4b) :: pdr_ptot,ion_ptot,dark_ptot
END MODULE maincode_module
