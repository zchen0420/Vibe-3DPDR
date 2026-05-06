module run_config_module

  use definitions
  implicit none

  type :: run_config
    character(len=128) :: input
    character(len=128) :: output
    integer :: level
    integer :: chemiterations
    integer :: itertot
    real(kind=dp) :: theta_crit
    real(kind=dp) :: rho_min
    real(kind=dp) :: rho_max
    real(kind=dp) :: zeta
    real(kind=dp) :: v_turb_inp
    real(kind=dp) :: dust_temperature
    real(kind=dp) :: end_time
    real(kind=dp) :: g2d
    real(kind=dp) :: metallicity
    real(kind=dp) :: omega
    real(kind=dp) :: grain_radius
    character(len=128) :: c12oinput
    character(len=128) :: ciiinput
    character(len=128) :: ciinput
    character(len=128) :: oiinput
    real(kind=dp) :: tguess
    real(kind=dp) :: tlow0
    real(kind=dp) :: thigh0
    real(kind=dp) :: tmin
    real(kind=dp) :: tmax
    real(kind=dp) :: fcrit
    real(kind=dp) :: tdiff
    character(len=3) :: fieldchoice
    real(kind=dp) :: gext(1)
  end type run_config

contains

  subroutine read_run_config(config_file, config)
    character(len=*), intent(in) :: config_file
    type(run_config), intent(out) :: config
    integer, parameter :: params_unit = 12

    open(unit=params_unit, file=config_file, status='old', action='read')

    call skip_lines(params_unit, 3)
    read(params_unit,*) config%input
    read(params_unit,*) config%output
    read(params_unit,*) config%level
    read(params_unit,*) config%theta_crit
    read(params_unit,*) config%chemiterations
    read(params_unit,*) config%itertot
    read(params_unit,*) config%rho_min
    read(params_unit,*) config%rho_max
    read(params_unit,*) config%zeta
    read(params_unit,*) config%v_turb_inp
    read(params_unit,*) config%dust_temperature
    read(params_unit,*) config%end_time
    read(params_unit,*) config%g2d
    read(params_unit,*) config%metallicity
    read(params_unit,*) config%omega
    read(params_unit,*) config%grain_radius

    call skip_lines(params_unit, 3)
    read(params_unit,*) config%c12oinput
    read(params_unit,*) config%ciiinput
    read(params_unit,*) config%ciinput
    read(params_unit,*) config%oiinput

    call skip_lines(params_unit, 3)
    read(params_unit,*) config%tguess
    read(params_unit,*) config%tlow0
    read(params_unit,*) config%thigh0
    read(params_unit,*) config%tmin
    read(params_unit,*) config%tmax
    read(params_unit,*) config%fcrit
    read(params_unit,*) config%tdiff

    call skip_lines(params_unit, 3)
    read(params_unit,*) config%fieldchoice
    read(params_unit,*) config%gext

    close(params_unit)
  end subroutine read_run_config

  subroutine skip_lines(unit, nlines)
    integer, intent(in) :: unit
    integer, intent(in) :: nlines
    integer :: line

    do line = 1, nlines
      read(unit,*)
    end do
  end subroutine skip_lines

  integer function count_data_records(data_file)
  character(len=*), intent(in) :: data_file
  integer, parameter :: data_unit = 91
  integer :: ios
  real(kind=dp) :: dummy

  count_data_records = 0
  open(unit=data_unit, file=data_file, status='old', action='read')
  do
    read(data_unit,*,iostat=ios) dummy
    if (ios /= 0) exit
    count_data_records = count_data_records + 1
  end do
  close(data_unit)
end function count_data_records

end module run_config_module


subroutine readparams(config_file)

  !T.Bisbas, T.Bell
  use definitions
  use healpix_types
  use run_config_module
  use coolants_module, only : coolant_count, coolant_c12o, coolant_cii, coolant_ci, coolant_oi, &
      &coolant_default_nlevels, coolant_default_ntemps
  use maincode_module, only : runtime, level, theta_crit, rho_min, rho_max,&
      & fieldchoice, gext, nspec, nreac, maxpoints, &
      & tlow0, thigh0, tmin, tmax, fcrit, tdiff, dust_temperature,&
      & end_time, coolant
  use global_module, only : g2d, metallicity, omega, grain_radius
  implicit none

  character(len=*), intent(in) :: config_file
  integer(kind=i4b) :: coolant_id
  type(run_config) :: config

  call read_run_config(config_file, config)

  runtime%input_file = config%input
  runtime%output_prefix = config%output
  level = config%level
  theta_crit = config%theta_crit
  runtime%chemistry_iterations = config%chemiterations
  runtime%total_iterations = config%itertot
  rho_min = config%rho_min
  rho_max = config%rho_max
  runtime%cosmic_ray_ionization_rate = config%zeta/1.3d-17
  runtime%input_turbulent_velocity = config%v_turb_inp
  dust_temperature = config%dust_temperature
  end_time = config%end_time
  g2d = config%g2d
  metallicity = config%metallicity
  omega = config%omega
  grain_radius = config%grain_radius
  if (.not.allocated(coolant)) allocate(coolant(1:coolant_count))
  coolant(coolant_c12o)%input_file = config%c12oinput
  coolant(coolant_cii)%input_file = config%ciiinput
  coolant(coolant_ci)%input_file = config%ciinput
  coolant(coolant_oi)%input_file = config%oiinput
  runtime%temperature_guess = config%tguess
  tlow0 = config%tlow0
  thigh0 = config%thigh0
  tmin = config%tmin
  tmax = config%tmax
  fcrit = config%fcrit
  tdiff = config%tdiff
  fieldchoice = config%fieldchoice
  gext = config%gext

  maxpoints = 600
  runtime%av_scale = 6.289e-22*metallicity
  runtime%uv_scale = 3.02

#ifdef REDUCED
  nspec = count_data_records('data/species_reduced.d')
  nreac = count_data_records('data/rates_reduced.d')
#endif

#ifdef FULL
  nspec = count_data_records('data/species_full.d')
  nreac = count_data_records('data/rates_full.d')
#endif

#ifdef MYNETWORK
  nspec = count_data_records('data/species_mynetwork.d')
  nreac = count_data_records('data/rates_mynetwork.d')
#endif

  do coolant_id = 1, coolant_count
    coolant(coolant_id)%nlevels = coolant_default_nlevels(coolant_id)
    coolant(coolant_id)%ntemperatures = coolant_default_ntemps(coolant_id)
  end do

  return
end subroutine readparams
