MODULE maincode_module
  USE ISO_C_BINDING
  use definitions
  use healpix_types

  integer(kind=I4B) :: i            ! counter
  integer :: ii           ! counter
  integer(kind=I4B) :: k            ! counter
  integer(kind=I4B) :: kk           ! counter
  integer(kind=I4B) :: p            ! counter
  integer(kind=I4B) :: ipix         ! pix id
  integer(kind=I4B) :: level        ! current level
  integer(kind=I4B) :: nrays        ! no. of rays on current level
  integer(kind=I4B) :: nside        ! refer to healpix manual
  integer(kind=I4B) :: itot         ! total number of grid points
  integer(kind=I4B) :: ktot         ! total number of grid points
  integer(kind=I4B) :: tot_eval     ! total number of evaluation points
  INTEGER(KIND=I4B) :: NLEV,NTEMP
  integer(kind=i4b) :: iteration, ITERTOT
  integer(kind=i4b) :: iterstep     ! output interval (per how many iterations)
  integer(kind=i4b) :: NSPEC, NREAC
  integer(kind=i4b) :: CII_NLEV, CII_NTEMP   !CII cooling variables
  integer(kind=i4b) :: CI_NLEV, CI_NTEMP     !CI cooling variables
  integer(kind=i4b) :: OI_NLEV, OI_NTEMP     !OI cooling variables
  integer(kind=i4b) :: C12O_NLEV, C12O_NTEMP !12CO cooling variables
  integer(kind=i4b) :: maxpoints
  integer(kind=i4b) :: suma
  integer :: CHEMITERATIONS
  integer :: pdr_tot

  logical::writeiterations
  logical::SPH

  integer(kind=I4B), allocatable :: rb(:)       ! ray with the IDs of the input file
  integer(kind=I4B), allocatable :: rrb(:)       ! ray with the IDs of the input file

  integer(kind=I4B) :: iwork(1:1000)
  real(kind=dp) :: rwork(1:1000000)

  real(kind=dp) :: ZETA
  real(kind=DP) :: angle_los   ! line-of-sight angle
  real(kind=dp) :: radius      ! maximum distance from the origin(1:3)
  real(kind=dp) :: theta_crit  ! critical theta angle to produce evaluation point
  real(kind=DP) :: rvec(1:3)   ! local 3 column array for the x,y,z coordinates of the grid point
  real(kind=DP) :: theta, phi  ! used from healpix to convert cartesian -> spherical
  real(kind=dp) :: adaptive_step !adaptive step for calculations of UVfield and column density.
  real(kind=dp) :: rhs1        ! right-hand-side 1 (to make integration calculation simpler)
  real(kind=DP) :: rhs2        ! -ditto-
  real(kind=DP) :: points      ! dummy variable used to find itot
  real(kind=DP) :: origin(1:3) ! grid point at which we perform calculation
  real(kind=DP) :: healpixvector(1:3) ! Healpix rays
  real(kind=dp) :: ENERGY,WEIGHT,EINSTEINA,FREQUENCY
  real(kind=dp) :: Tguess
  real(kind=dp) :: n_H
  real(kind=dp) :: CII_Z_function
  real(kind=dp) :: CI_Z_function
  real(kind=dp) :: OI_Z_function
  real(kind=dp) :: C12O_Z_function
  real(kind=dp) :: Z_increment
  real(kind=dp) :: frac1, frac2, frac3, tau_increment
  real(kind=dp) :: beta_ij_ray, beta_ij_sum, beta_ij
  real(kind=dp) :: TPOP, TMP2, BB_ij
  real(kind=dp) :: relch, v_turb, v_turb_inp
  real(kind=dp) :: Tlow0
  real(kind=dp) :: Thigh0
  real(kind=dp) :: Tmin
  real(kind=dp) :: Tmax
  real(kind=dp) :: Fcrit
  real(kind=dp) :: Tdiff
  real(kind=dp) :: dust_temperature
  !  real(kind=dp) :: h
  real(kind=dp) :: avmax

  character(len=3) :: fieldchoice
  real(kind=dp) :: Gext(1)
  real(kind=dp) :: AV_fac, UV_fac

  real(kind=DP), allocatable :: vectors(:,:) ! Healpix rays
  real(kind=DP), allocatable :: ep(:,:)            ! evaluation point along each ray (local)
  real(kind=dp), allocatable :: ra(:)              ! distance
  real(kind=dp), allocatable :: density(:)         ! density of each grid point - from input
  real(kind=DP), allocatable :: c_dens(:)          ! column density
  real(kind=dp), allocatable :: rra(:)             ! distance
  real(kind=dp), allocatable :: COEFF(:)
  real(kind=dp), allocatable :: ENERGIES(:), WEIGHTS(:)
  real(kind=dp), allocatable :: A_COEFFS(:,:), B_COEFFS(:,:), C_COEFFS(:,:)
  real(kind=dp), allocatable :: FREQUENCIES(:,:), TEMPERATURES(:)
  real(kind=dp), allocatable :: H_COL(:,:,:)
  real(kind=dp), allocatable :: EL_COL(:,:,:)
  real(kind=dp), allocatable :: HE_COL(:,:,:)
  real(kind=dp), allocatable :: H2_COL(:,:,:)
  real(kind=dp), allocatable :: PH2_COL(:,:,:)
  real(kind=dp), allocatable :: OH2_COL(:,:,:)
  real(kind=dp), allocatable :: tau_ij(:)
  real(kind=dp), allocatable :: field(:,:)
  !CII cooling variables
  real(kind=dp) :: CII_COOLING
  real(kind=dp), allocatable :: CII_ENERGIES(:),CII_WEIGHTS(:)
  real(kind=dp), allocatable :: CII_A_COEFFS(:,:),CII_B_COEFFS(:,:)
  real(kind=dp), allocatable :: CII_FREQUENCIES(:,:)
  real(kind=dp), allocatable :: CII_TEMPERATURES(:,:), CII_HP(:,:,:)
  real(kind=dp), allocatable :: CII_H(:,:,:),CII_EL(:,:,:)
  real(kind=dp), allocatable :: CII_HE(:,:,:),CII_H2(:,:,:)
  real(kind=dp), allocatable :: CII_PH2(:,:,:),CII_OH2(:,:,:)
  !CI cooling variables
  real(kind=dp) :: CI_COOLING
  real(kind=dp), allocatable :: CI_ENERGIES(:),CI_WEIGHTS(:)
  real(kind=dp), allocatable :: CI_A_COEFFS(:,:),CI_B_COEFFS(:,:)
  real(kind=dp), allocatable :: CI_FREQUENCIES(:,:)
  real(kind=dp), allocatable :: CI_TEMPERATURES(:,:),CI_HP(:,:,:)
  real(kind=dp), allocatable :: CI_H(:,:,:),CI_EL(:,:,:)
  real(kind=dp), allocatable :: CI_HE(:,:,:),CI_H2(:,:,:)
  real(kind=dp), allocatable :: CI_PH2(:,:,:),CI_OH2(:,:,:)
  !OI cooling variables
  real(kind=dp) :: OI_COOLING
  real(kind=dp), allocatable :: OI_ENERGIES(:),OI_WEIGHTS(:)
  real(kind=dp), allocatable :: OI_A_COEFFS(:,:),OI_B_COEFFS(:,:)
  real(kind=dp), allocatable :: OI_FREQUENCIES(:,:)
  real(kind=dp), allocatable :: OI_TEMPERATURES(:,:),OI_HP(:,:,:)
  real(kind=dp), allocatable :: OI_H(:,:,:),OI_EL(:,:,:)
  real(kind=dp), allocatable :: OI_HE(:,:,:),OI_H2(:,:,:)
  real(kind=dp), allocatable :: OI_PH2(:,:,:),OI_OH2(:,:,:)
  !12CO cooling variables
  real(kind=dp) :: C12O_COOLING
  real(kind=dp), allocatable :: C12O_ENERGIES(:),C12O_WEIGHTS(:)
  real(kind=dp), allocatable :: C12O_A_COEFFS(:,:),C12O_B_COEFFS(:,:)
  real(kind=dp), allocatable :: C12O_FREQUENCIES(:,:)
  real(kind=dp), allocatable :: C12O_TEMPERATURES(:,:),C12O_HP(:,:,:)
  real(kind=dp), allocatable :: C12O_H(:,:,:),C12O_EL(:,:,:)
  real(kind=dp), allocatable :: C12O_HE(:,:,:),C12O_H2(:,:,:)
  real(kind=dp), allocatable :: C12O_PH2(:,:,:),C12O_OH2(:,:,:)

  character(len=128) :: input
  character(len=128) :: inputchem
  character(len=128) :: C12Oinput, CIIinput, CIinput, OIinput
  character(len=128)  :: output

  type columndens_node
    real(kind=dp), pointer :: columndens_point(:,:)
  end type columndens_node
  type(columndens_node), allocatable :: column(:)
  integer(kind=i4b), allocatable :: DUPLICATE(:)

  real(kind=dp),allocatable :: ALPHA(:),BETA(:),GAMMA(:),RATE(:),RTMIN(:),RTMAX(:)
  CHARACTER(len=10), allocatable :: REACTANT(:,:),PRODUCT(:,:)
  real(kind=dp), allocatable :: MASS(:),dummyabundance(:)!,abundances(:,:)
  CHARACTER(len=10), allocatable :: SPECIES(:)

  real(kind=dp),bind(c,name='maincode_module_mp_start_time_'):: start_time
  real(kind=dp),bind(c,name='maincode_module_mp_end_time_'):: end_time
  integer(kind=i4b) :: status

  real(kind=dp), allocatable :: dusttemperature(:)
  real(kind=dp), allocatable :: gastemperature(:)
  real(kind=dp), allocatable :: previousgastemperature(:)

#ifdef THERMALBALANCE
  real(kind=dp), allocatable :: Fmean(:)
  real(kind=dp), allocatable :: Fratio(:)
  real(kind=dp), allocatable :: Tlow(:)
  real(kind=dp), allocatable :: Thigh(:)
#endif

  type pdr_node
    integer(kind=i4b), pointer :: epray(:)        !population of evaluation points per ray
    integer(kind=i4b), pointer :: projected(:,:)  !ID of projected grid points in the line of sight
    integer(kind=i4b), pointer :: raytype(:)      !raytype
    real(kind=dp), pointer :: epoint(:,:,:)       !co-ordinates of each evaluation point
    real(kind=dp), pointer :: columndensity(:)    !column density
    real(kind=dp), pointer :: AV(:)               !AV
    real(kind=dp), pointer :: rad_surface(:)      !rad_surface
    real(kind=dp), pointer :: abundance(:)        !abundance of species
    real(kind=dp), pointer :: CII_pop(:)          !CII population density
    real(kind=dp), pointer :: CI_pop(:)           !CI population density
    real(kind=dp), pointer :: OI_pop(:)           !OI population density
    real(kind=dp), pointer :: C12O_pop(:)         !CO population density
    real(kind=dp), pointer :: CII_line(:,:)       !CII line cooling
    real(kind=dp), pointer :: CI_line(:,:)        !CI line cooling
    real(kind=dp), pointer :: OI_line(:,:)        !OI line cooling
    real(kind=dp), pointer :: C12O_line(:,:)      !C12O line cooling
    !===============
    real(kind=dp), pointer :: CII_optdepth(:,:,:)       !CII line cooling
    real(kind=dp), pointer :: CI_optdepth(:,:,:)        !CI line cooling
    real(kind=dp), pointer :: OI_optdepth(:,:,:)        !OI line cooling
    real(kind=dp), pointer :: C12O_optdepth(:,:,:)      !C12O line cooling
    real(kind=dp), pointer :: C12O_beta(:,:,:)      !C12O line escape probability
    real(kind=dp), pointer :: CII_beta(:,:,:)      !C12O line escape probability
    real(kind=dp), pointer :: CI_beta(:,:,:)      !C12O line escape probability
    real(kind=dp), pointer :: OI_beta(:,:,:)      !C12O line escape probability
    real(kind=dp), pointer :: CII_CCOEFFS(:,:)
    real(kind=dp), pointer :: check_step(:)
    real(kind=dp), pointer :: check_coef(:,:,:)
    !===============
    real(kind=dp) :: UVfield                      !attenuated UV field of element
    real(kind=dp) :: rho                          !density of element
    real(kind=dp) :: smoo                         !smoothing length (if SPH = .TRUE. in params.dat)
    real(kind=dp) :: x,y,z                        !position of element
    integer(kind=i4b) :: etype                    !element type (1 = PDR, 2 = ION, 3 = DARK)
#ifdef DUST2
    real(kind=dp) :: dust_t
#endif
  end type pdr_node
  type (pdr_node), allocatable :: pdr(:)      !main 3DPDR array for each grid point p
  real(kind=dp), allocatable :: CII_solution(:,:)
  real(kind=dp), allocatable :: CI_solution(:,:)
  real(kind=dp), allocatable :: OI_solution(:,:)
  real(kind=dp), allocatable :: C12O_solution(:,:)

  integer(kind=i4b)::levpop_iteration
  real(kind=dp) :: rho_min
  real(kind=dp) :: rho_max
#ifdef GUESS_TEMP
  real(kind=dp),allocatable :: Tmin_array(:)
  real(kind=dp),allocatable :: Tmax_array(:)
#endif

  !================================
  !================================
  !================================
  real::time_evalpoints,t2,t3,t3b,t4 !times for the cpu_time subroutine
  real(kind=dp)::percentage,thermal_percentage
  real(kind=dp)::OI_percentage,C12O_percentage
  real(kind=dp)::CII_percentage,CI_percentage
  real(kind=dp)::levpop_percentage
  real(kind=dp) :: rad_tot
  real(kind=dp),allocatable :: CII_cool(:)
  real(kind=dp),allocatable :: CI_cool(:)
  real(kind=dp),allocatable :: OI_cool(:)
  real(kind=dp),allocatable :: C12O_cool(:)
  real(kind=dp),allocatable :: total_cooling_rate(:)
  integer(kind=i4b) :: referee, id
#ifdef THERMALBALANCE
  logical,allocatable :: converged(:),doleveltmin(:)
  !logical :: show_balance
#endif
  logical,allocatable :: level_converged(:)
  logical :: level_conv,first_time
  logical,allocatable :: CII_conv(:),CI_conv(:),OI_conv(:),C12O_conv(:)
  real(kind=dp),allocatable :: CII_RELCH(:,:),CI_RELCH(:,:),OI_RELCH(:,:),C12O_RELCH(:,:)
  logical :: relch_conv
  integer(kind=i4b) :: CII_i,CI_i,OI_i,C12O_i
  logical,allocatable::expanded(:)
  integer(kind=i4b) :: grand_ptot
  integer(kind=i4b) :: pdr_ptot,ion_ptot,dark_ptot
  integer(kind=i4b),allocatable :: IDlist_pdr(:),IDlist_ion(:),IDlist_dark(:)
  real(kind=DP), allocatable :: pdrpoint(:,:)      ! coordinates of pdr element
  integer(kind=i4b) :: pp
  real(kind=dp) :: xpos,ypos,zpos,denst
  !================================
  !================================
  !================================

END MODULE maincode_module
