module escape_probability_module
  use coolants_module, only : coolant_data
  use definitions, only : dp
  use global_module, only : g2d
  use healpix_types, only : c, hp, i4b, kb, mp, pc, pi
  use maincode_module, only : geometry
  implicit none

  private
  public :: lvg_local_conditions
  public :: calculate_lvg_transition_rates

  type :: lvg_local_conditions
    real(kind=dp) :: gas_temperature
    real(kind=dp) :: dust_temperature
    real(kind=dp) :: turbulent_velocity
    real(kind=dp) :: gas_density
    real(kind=dp) :: metallicity
  end type lvg_local_conditions

contains

  subroutine calculate_lvg_transition_rates(coolant_table, conditions, ray_point_count, level_population, &
        &evaluation_points, evaluation_populations, collision_coefficients, transition_rates, cooling_rate, &
        &line_emission, optical_depth, beta)
    type(coolant_data), intent(in) :: coolant_table
    type(lvg_local_conditions), intent(in) :: conditions
    integer(kind=i4b), intent(in) :: ray_point_count(0:)
    real(kind=dp), intent(in) :: level_population(1:)
    real(kind=dp), intent(in) :: evaluation_points(1:,0:,0:)
    real(kind=dp), intent(in) :: evaluation_populations(0:,0:,1:)
    real(kind=dp), intent(in) :: collision_coefficients(1:,1:)
    real(kind=dp), intent(inout) :: transition_rates(1:,1:)
    real(kind=dp), intent(out) :: cooling_rate
    real(kind=dp), intent(out) :: line_emission(1:,1:)
    real(kind=dp), intent(out) :: optical_depth(1:,1:,0:)
    real(kind=dp), intent(out) :: beta(1:,1:,0:)

    integer(kind=i4b) :: eval_index
    integer(kind=i4b) :: lower_level
    integer(kind=i4b) :: nlevels
    integer(kind=i4b) :: ray_index
    integer(kind=i4b) :: ray_lower
    integer(kind=i4b) :: ray_upper
    integer(kind=i4b) :: upper_level
    real(kind=dp) :: background_intensity
    real(kind=dp) :: beta_ij
    real(kind=dp) :: beta_ij_sum
    real(kind=dp) :: dust_background
    real(kind=dp) :: emissivity
    real(kind=dp) :: line_prefactor
    real(kind=dp) :: line_width_factor
    real(kind=dp) :: source_function
    real(kind=dp) :: step_length
    real(kind=dp) :: tau_increment
    real(kind=dp) :: transition_population
    real(kind=dp) :: upper_population
    real(kind=dp), allocatable :: beta_by_ray(:)
    real(kind=dp), allocatable :: field(:,:)
    real(kind=dp), allocatable :: tau_by_ray(:)

    nlevels = size(level_population)
    ray_lower = lbound(ray_point_count,1)
    ray_upper = ubound(ray_point_count,1)

    call assert_lvg_dimensions(coolant_table, nlevels, ray_lower, ray_upper, evaluation_points, &
        &evaluation_populations, collision_coefficients, transition_rates, line_emission, optical_depth, beta)

    line_emission = 0.0d0
    cooling_rate = 0.0d0
    allocate(tau_by_ray(ray_lower:ray_upper))
    allocate(beta_by_ray(ray_lower:ray_upper))
    allocate(field(1:nlevels,1:nlevels))
    field = 0.0d0

    line_width_factor = 1.0d0/sqrt(8.0d0*kb*conditions%gas_temperature/pi/mp + &
        &conditions%turbulent_velocity**2)

    do upper_level=1,nlevels
      do lower_level=1,nlevels
        if (lower_level.ge.upper_level) exit

        tau_by_ray = 0.0d0
        beta_ij = 0.0d0
        beta_by_ray = 0.0d0
        beta_ij_sum = 0.0d0
        line_prefactor = (coolant_table%a_coeffs(upper_level,lower_level)*(c**3)) &
            &/(8.0d0*pi*(coolant_table%frequencies(upper_level,lower_level)**3))
        upper_population = 2.0d0*hp*(coolant_table%frequencies(upper_level,lower_level)**3)/(c**2)

        background_intensity = upper_population*(1.0d0/(exp(hp*coolant_table%frequencies(upper_level,lower_level) &
            &/kb/2.7d0)-1.0d0))
        emissivity = (2.0d0*2.0d-12*conditions%gas_density*conditions%metallicity*100.0d0/g2d) &
            &*(0.01d0*(1.3d0*coolant_table%frequencies(upper_level,lower_level)/3.0d11))
        dust_background = upper_population*(1.0d0/(exp(hp*coolant_table%frequencies(upper_level,lower_level) &
            &/kb/conditions%dust_temperature)-1.0d0)*emissivity)
        background_intensity = background_intensity + dust_background

        if (level_population(upper_level).eq.0.0d0) then
          source_function = 0.0d0
          beta_ij = 1.0d0
          goto 2
        end if

        transition_population = (level_population(lower_level)*coolant_table%weights(upper_level)) &
            &/(level_population(upper_level)*coolant_table%weights(lower_level))-1.0d0
        if (abs(transition_population).lt.1.0d-50) then
          source_function = hp*coolant_table%frequencies(upper_level,lower_level) &
              &*level_population(upper_level)*coolant_table%a_coeffs(upper_level,lower_level)/(4.0d0*pi)
          beta_ij = 1.0d0
          goto 1
        else
          source_function = upper_population/transition_population
        end if

        do ray_index=ray_lower,ray_upper
#ifdef PSEUDO_1D
          if (ray_index.ne.6) then
            tau_by_ray(ray_index) = 1.0d50
          else
#endif
#ifdef PSEUDO_2D
            if (abs(geometry%ray_vectors(3,ray_index)).gt.1.0d-10) tau_by_ray(ray_index) = 1.0d50
#endif
            do eval_index=1,ray_point_count(ray_index)
              tau_increment = ((evaluation_populations(ray_index,eval_index-1,lower_level) &
                  &*coolant_table%weights(upper_level) &
                  &-evaluation_populations(ray_index,eval_index-1,upper_level)*coolant_table%weights(lower_level)) &
                  &+(evaluation_populations(ray_index,eval_index,lower_level)*coolant_table%weights(upper_level) &
                  &-evaluation_populations(ray_index,eval_index,upper_level)*coolant_table%weights(lower_level))) &
                  &/(2.0d0*coolant_table%weights(lower_level))
              step_length = sqrt((evaluation_points(1,ray_index,eval_index-1) &
                  &-evaluation_points(1,ray_index,eval_index))**2 &
                  &+(evaluation_points(2,ray_index,eval_index-1) &
                  &-evaluation_points(2,ray_index,eval_index))**2 &
                  &+(evaluation_points(3,ray_index,eval_index-1) &
                  &-evaluation_points(3,ray_index,eval_index))**2)
              tau_by_ray(ray_index) = tau_by_ray(ray_index) &
                  &+line_prefactor*line_width_factor*tau_increment*step_length*pc
            end do
#ifdef PSEUDO_1D
          end if
#endif
          beta_by_ray(ray_index) = escape_probability_for_tau(tau_by_ray(ray_index))
          optical_depth(upper_level,lower_level,ray_index) = tau_by_ray(ray_index)
          beta(upper_level,lower_level,ray_index) = beta_by_ray(ray_index)
        end do

        beta_ij_sum = sum(beta_by_ray)
#ifdef PSEUDO_1D
        beta_ij = beta_ij_sum
#elif PSEUDO_2D
        beta_ij = beta_ij_sum / 4.0d0
#else
        beta_ij = beta_ij_sum / real(size(ray_point_count),kind=dp)
#endif

        1       continue
        line_emission(upper_level,lower_level) = coolant_table%a_coeffs(upper_level,lower_level)*hp &
            &*coolant_table%frequencies(upper_level,lower_level)*level_population(upper_level)*beta_ij &
            &*(source_function-background_intensity)/source_function
        cooling_rate = cooling_rate + line_emission(upper_level,lower_level)
        2       continue
        field(upper_level,lower_level) = (1.0d0-beta_ij)*source_function + beta_ij*background_intensity
        field(lower_level,upper_level) = field(upper_level,lower_level)
      end do
    end do

    do upper_level=1,nlevels
      do lower_level=1,nlevels
        transition_rates(upper_level,lower_level) = coolant_table%a_coeffs(upper_level,lower_level) &
            &+coolant_table%b_coeffs(upper_level,lower_level)*field(upper_level,lower_level) &
            &+collision_coefficients(upper_level,lower_level)
        if (abs(transition_rates(upper_level,lower_level)).lt.1.0d-50) then
          transition_rates(upper_level,lower_level) = 0.0d0
        end if
      end do
    end do

    deallocate(tau_by_ray)
    deallocate(beta_by_ray)
    deallocate(field)
  end subroutine calculate_lvg_transition_rates

  subroutine assert_lvg_dimensions(coolant_table, nlevels, ray_lower, ray_upper, evaluation_points, &
        &evaluation_populations, collision_coefficients, transition_rates, line_emission, optical_depth, beta)
    type(coolant_data), intent(in) :: coolant_table
    integer(kind=i4b), intent(in) :: nlevels
    integer(kind=i4b), intent(in) :: ray_lower
    integer(kind=i4b), intent(in) :: ray_upper
    real(kind=dp), intent(in) :: evaluation_points(1:,0:,0:)
    real(kind=dp), intent(in) :: evaluation_populations(0:,0:,1:)
    real(kind=dp), intent(in) :: collision_coefficients(1:,1:)
    real(kind=dp), intent(in) :: transition_rates(1:,1:)
    real(kind=dp), intent(in) :: line_emission(1:,1:)
    real(kind=dp), intent(in) :: optical_depth(1:,1:,0:)
    real(kind=dp), intent(in) :: beta(1:,1:,0:)

    if (coolant_table%nlevels.ne.nlevels) stop 'calculate_lvg_transition_rates received inconsistent levels'
    if (size(collision_coefficients,1).ne.nlevels .or. size(collision_coefficients,2).ne.nlevels) then
      stop 'calculate_lvg_transition_rates received inconsistent collision coefficients'
    end if
    if (size(transition_rates,1).ne.nlevels .or. size(transition_rates,2).ne.nlevels) then
      stop 'calculate_lvg_transition_rates received inconsistent transition rates'
    end if
    if (size(line_emission,1).ne.nlevels .or. size(line_emission,2).ne.nlevels) then
      stop 'calculate_lvg_transition_rates received inconsistent line output'
    end if
    if (lbound(evaluation_points,2).ne.ray_lower .or. ubound(evaluation_points,2).ne.ray_upper) then
      stop 'calculate_lvg_transition_rates received inconsistent evaluation point rays'
    end if
    if (lbound(evaluation_populations,1).ne.ray_lower .or. ubound(evaluation_populations,1).ne.ray_upper) then
      stop 'calculate_lvg_transition_rates received inconsistent evaluation population rays'
    end if
    if (lbound(optical_depth,3).ne.ray_lower .or. ubound(optical_depth,3).ne.ray_upper) then
      stop 'calculate_lvg_transition_rates received inconsistent optical-depth rays'
    end if
    if (lbound(beta,3).ne.ray_lower .or. ubound(beta,3).ne.ray_upper) then
      stop 'calculate_lvg_transition_rates received inconsistent beta rays'
    end if
  end subroutine assert_lvg_dimensions

  real(kind=dp) function escape_probability_for_tau(tau)
  real(kind=dp), intent(in) :: tau

  if (tau.lt.-5.0d0) then
    escape_probability_for_tau = (1.0d0-exp(5.0d0))/(-5.0d0)
  else if (abs(tau).lt.1.0d-8) then
    escape_probability_for_tau = 1.0d0
  else
    escape_probability_for_tau = (1.0d0-exp(-tau))/tau
  end if
end function escape_probability_for_tau

end module escape_probability_module
