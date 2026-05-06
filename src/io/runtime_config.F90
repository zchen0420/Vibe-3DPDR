module runtime_config_module
  use healpix_types
  use maincode_module
  use uclpdr_module, only : start
  use global_module

contains

  subroutine load_runtime_configuration(config_file)
    character(len=*), intent(inout) :: config_file

    config_file = 'configs/default.params'
    if (command_argument_count() >= 1) call get_command_argument(1, config_file)

    call print_banner

    write(6,*) 'Reading runtime%input_file [',trim(config_file),']'
    call readparams(trim(config_file))
    runtime%turbulent_velocity=runtime%input_turbulent_velocity*1.0D5

    call validate_runtime_configuration
    call print_runtime_configuration
    call print_enabled_features
  end subroutine load_runtime_configuration

  subroutine print_banner
    write(6,*) '=============================================================================='
    write(6,*) '*********     *********             *********     *********      *********'
    write(6,*) '        **    **       **           **       **   **       **    **       **'
    write(6,*) '         **   **        **          **        **  **        **   **        **'
    write(6,*) '        **    **         **         **       **   **         **  **       **'
    write(6,*) '  *******     **         **  *****  *********     **         **  *********'
    write(6,*) '        **    **         **         **            **         **  **     **'
    write(6,*) '         **   **        **          **            **        **   **      **'
    write(6,*) '        **    **       **           **            **       **    **       **'
    write(6,*) '*********     *********             **            *********      **        **'
    write(6,*) '=============================================================================='
    write(6,*) '********************   Coders:   T.G.Bisbas, T.A.Bell   **********************'
    write(6,*) '*************   Collaborators:   S.Viti, J.Yates, M.Barlow   *****************'
    write(6,*) '*************************        Version 1.0          ************************'
    write(6,*) '=============================================================================='
    write(6,*) ''
  end subroutine print_banner

  subroutine validate_runtime_configuration
#ifdef PSEUDO_1D
    if (level.gt.0) STOP "HEALPix level must be set to 0 in PSEUDO_1D mode"
#endif
#ifdef PSEUDO_2D
    if (level.gt.0) STOP "HEALPix level must be set to 0 in PSEUDO_1D mode"
#endif
    if (theta_crit.ge.pi/2.0D0) stop 'theta_crit must be less than pi/2'
  end subroutine validate_runtime_configuration

  subroutine print_runtime_configuration
    write(6,*) 'Input file:               ',runtime%input_file
    write(6,*) 'HEALPix level:            ',level
    write(6,*) 'Theta critical:           ',theta_crit
    write(6,*) 'Angle between rays:       ',sqrt(pi/3.0D0/4.0D0**(real(level)))
    write(6,*) 'Maxpoints                 ',maxpoints
    write(6,*) 'Guess Temperature (K):    ',runtime%temperature_guess
    write(6,*) 'Dust  Temperature (K):    ',dust_temperature
    write(6,*) 'Turbulent velocity (cm/s):',runtime%turbulent_velocity
    write(6,*) 'minimum density (cm^-3):  ',rho_min
    write(6,*) 'maximum density (cm^-3):  ',rho_max
#ifdef THERMALBALANCE
    write(6,*) 'Tlow:                     ',Tlow0
    write(6,*) 'Thigh:                    ',Thigh0
#endif
    write(6,*) 'Tmin:                     ',Tmin
    write(6,*) 'Tmax:                     ',Tmax
    write(6,*) 'Fcrit:                    ',Fcrit
    write(6,*) 'Tdiff:                    ',Tdiff
    start = .true.
    write(6,*) 'Form of field:            ',fieldchoice
    write(6,*) 'Gext:                     ',Gext(1)
    write(6,*) 'AV factor:                ',runtime%av_scale
    write(6,*) 'UV factor:                ',runtime%uv_scale
#ifdef REDUCED
    write(6,*) 'Chemical network:           REDUCED'
#endif
#ifdef FULL
    write(6,*) 'Chemical network:           FULL'
#endif
#ifdef MYNETWORK
    write(6,*) 'Chemical network:           MYNETWORK'
#endif
    write(6,*) 'Number of chemistry%network%species:        ',nspec
    write(6,*) 'Number of reactions:      ',nreac
    write(6,*) 'Total iterations:         ',runtime%total_iterations
    write(6,*) 'Output interval / iter.:  ',runtime%output_interval
    write(6,*) 'runtime%chemistry_iterations:           ',runtime%chemistry_iterations
    write(6,*) 'runtime%cosmic_ray_ionization_rate:                     ',runtime%cosmic_ray_ionization_rate*1.3d-17
    write(6,*) 'Gas-to-dust               ',g2d!*100.0
    write(6,*) 'Metallicity               ',metallicity
    write(6,*) 'Omega                     ',omega
    write(6,*) 'Grain radius              ',grain_radius
    close(12)
    write(6,*) '============================================='
    write(6,*) ''
  end subroutine print_runtime_configuration

  subroutine print_enabled_features
    write(6,*) '------FLAGS-----'
#ifdef THERMALBALANCE
    write(6,*) 'THERMALBALANCE'
#endif
#ifdef PSEUDO_1D
    write(6,*) 'PSEUDO_1D'
#elif PSEUDO_2D
    write(6,*) 'PSEUDO_2D'
#else
    write(6,*) 'FULL_3D'
#endif
#ifdef DUST
    write(6,*) 'DUST'
#endif
#ifdef DUST2
    write(6,*) 'DUST2'
#endif
#ifdef CO_FIX
    write(6,*) 'CO_FIX'
#endif
#ifdef H2FORM
    write(6,*) 'H2FORM'
#endif
#ifdef TEMP_FIX
    write(6,*) 'TEMP_FIX'
#endif
#ifdef GUESS_TEMP
    write(6,*) 'GUESS_TEMP'
#endif
    write(6,*) ''
  end subroutine print_enabled_features

end module runtime_config_module
