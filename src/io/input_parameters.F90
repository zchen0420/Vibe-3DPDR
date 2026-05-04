MODULE run_config_module

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
    character(len=128) :: C12Oinput
    character(len=128) :: CIIinput
    character(len=128) :: CIinput
    character(len=128) :: OIinput
    real(kind=dp) :: Tguess
    real(kind=dp) :: Tlow0
    real(kind=dp) :: Thigh0
    real(kind=dp) :: Tmin
    real(kind=dp) :: Tmax
    real(kind=dp) :: Fcrit
    real(kind=dp) :: Tdiff
    character(len=3) :: fieldchoice
    real(kind=dp) :: Gext(1)
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
    read(params_unit,*) config%C12Oinput
    read(params_unit,*) config%CIIinput
    read(params_unit,*) config%CIinput
    read(params_unit,*) config%OIinput

    call skip_lines(params_unit, 3)
    read(params_unit,*) config%Tguess
    read(params_unit,*) config%Tlow0
    read(params_unit,*) config%Thigh0
    read(params_unit,*) config%Tmin
    read(params_unit,*) config%Tmax
    read(params_unit,*) config%Fcrit
    read(params_unit,*) config%Tdiff

    call skip_lines(params_unit, 3)
    read(params_unit,*) config%fieldchoice
    read(params_unit,*) config%Gext

    close(params_unit)
  end subroutine read_run_config

  subroutine skip_lines(unit, nlines)
    integer, intent(in) :: unit
    integer, intent(in) :: nlines
    integer :: line

    do line = 1, nlines
      read(unit,*)
    enddo
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
    enddo
    close(data_unit)
  end function count_data_records

END MODULE run_config_module


SUBROUTINE readparams(config_file)

  !T.Bisbas, T.Bell
  use definitions
  use healpix_types
  use run_config_module
  use maincode_module, only : input, level, Tguess, v_turb_inp, iterstep,&
      & theta_crit, ITERTOT, output, rho_min, rho_max,&
      & fieldchoice, Gext, AV_fac, UV_fac, nspec, nreac, &
      & CII_NLEV, CII_NTEMP, CI_NLEV, CI_NTEMP, &
      & OI_NLEV, OI_NTEMP, C12O_NLEV, C12O_NTEMP, maxpoints, &
      & Tlow0, Thigh0, Tmin, Tmax, Fcrit, Tdiff, dust_temperature,writeiterations,&
      & chemiterations, zeta, end_time, C12Oinput, CIIinput, CIinput, OIinput
  use global_module, only : g2d, metallicity, omega, grain_radius
  implicit none

  character(len=*), intent(in) :: config_file
  type(run_config) :: config

  call read_run_config(config_file, config)

  input = config%input
  output = config%output
  level = config%level
  theta_crit = config%theta_crit
  chemiterations = config%chemiterations
  ITERTOT = config%itertot
  rho_min = config%rho_min
  rho_max = config%rho_max
  zeta = config%zeta/1.3d-17
  v_turb_inp = config%v_turb_inp
  dust_temperature = config%dust_temperature
  end_time = config%end_time
  g2d = config%g2d
  metallicity = config%metallicity
  omega = config%omega
  grain_radius = config%grain_radius
  C12Oinput = config%C12Oinput
  CIIinput = config%CIIinput
  CIinput = config%CIinput
  OIinput = config%OIinput
  Tguess = config%Tguess
  Tlow0 = config%Tlow0
  Thigh0 = config%Thigh0
  Tmin = config%Tmin
  Tmax = config%Tmax
  Fcrit = config%Fcrit
  Tdiff = config%Tdiff
  fieldchoice = config%fieldchoice
  Gext = config%Gext

  maxpoints = 600
  AV_fac = 6.289E-22*metallicity
  UV_fac = 3.02

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

  CII_NLEV = 5
  CII_NTEMP = 18
  CI_NLEV = 5
  CI_NTEMP = 29
  OI_NLEV = 5
  OI_NTEMP = 27
  C12O_NLEV = 41
  C12O_NTEMP = 25

  return
END SUBROUTINE readparams
