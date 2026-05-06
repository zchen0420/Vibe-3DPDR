module coolants_module
  use definitions, only : dp
  use healpix_types, only : i4b
  implicit none

  integer(kind=i4b), parameter :: COOLANT_COUNT = 4
  integer(kind=i4b), parameter :: COOLANT_CII = 1
  integer(kind=i4b), parameter :: COOLANT_CI = 2
  integer(kind=i4b), parameter :: COOLANT_OI = 3
  integer(kind=i4b), parameter :: COOLANT_C12O = 4
  integer(kind=i4b), parameter :: COLLISION_PARTNER_COUNT = 7
  integer(kind=i4b), parameter :: COLLIDER_H2 = 1
  integer(kind=i4b), parameter :: COLLIDER_PARA_H2 = 2
  integer(kind=i4b), parameter :: COLLIDER_ORTHO_H2 = 3
  integer(kind=i4b), parameter :: COLLIDER_ELECTRON = 4
  integer(kind=i4b), parameter :: COLLIDER_H = 5
  integer(kind=i4b), parameter :: COLLIDER_HE = 6
  integer(kind=i4b), parameter :: COLLIDER_PROTON = 7

  character(len=4), parameter :: COOLANT_LABELS(COOLANT_COUNT) = (/&
      &'CII ', 'CI  ', 'OI  ', 'CO  '/)
  integer(kind=i4b), parameter :: COOLANT_DEFAULT_NLEVELS_TABLE(COOLANT_COUNT) = (/5, 5, 5, 41/)
  integer(kind=i4b), parameter :: COOLANT_DEFAULT_NTEMPS_TABLE(COOLANT_COUNT) = (/18, 29, 27, 25/)

  type :: coolant_data
    character(len=128) :: input_file
    integer(kind=i4b) :: nlevels
    integer(kind=i4b) :: ntemperatures
    real(kind=dp), allocatable :: energies(:)
    real(kind=dp), allocatable :: weights(:)
    real(kind=dp), allocatable :: a_coeffs(:,:)
    real(kind=dp), allocatable :: b_coeffs(:,:)
    real(kind=dp), allocatable :: frequencies(:,:)
    real(kind=dp), allocatable :: collision_temperatures(:,:)
    real(kind=dp), allocatable :: collision_rates(:,:,:,:)
  end type coolant_data

  type :: point_coolant_state
    real(kind=dp), allocatable :: population(:)
    real(kind=dp), allocatable :: line(:,:)
    real(kind=dp), allocatable :: optical_depth(:,:,:)
  end type point_coolant_state

  type :: coolant_iteration_state
    real(kind=dp), allocatable :: solution(:,:)
    real(kind=dp), allocatable :: cooling_rate(:)
    logical, allocatable :: converged(:)
    real(kind=dp), allocatable :: relative_change(:,:)
  end type coolant_iteration_state

contains

  function coolant_label(coolant_id) result(label)
    integer(kind=i4b), intent(in) :: coolant_id
    character(len=4) :: label

    if (coolant_id.lt.1 .or. coolant_id.gt.COOLANT_COUNT) then
      stop 'Invalid coolant id'
    endif
    label = COOLANT_LABELS(coolant_id)
  end function coolant_label

  integer(kind=i4b) function coolant_default_nlevels(coolant_id)
    integer(kind=i4b), intent(in) :: coolant_id

    call assert_valid_coolant(coolant_id)
    coolant_default_nlevels = COOLANT_DEFAULT_NLEVELS_TABLE(coolant_id)
  end function coolant_default_nlevels

  integer(kind=i4b) function coolant_default_ntemps(coolant_id)
    integer(kind=i4b), intent(in) :: coolant_id

    call assert_valid_coolant(coolant_id)
    coolant_default_ntemps = COOLANT_DEFAULT_NTEMPS_TABLE(coolant_id)
  end function coolant_default_ntemps

  subroutine assert_valid_coolant(coolant_id)
    integer(kind=i4b), intent(in) :: coolant_id

    if (coolant_id.lt.1 .or. coolant_id.gt.COOLANT_COUNT) then
      stop 'Invalid coolant id'
    endif
  end subroutine assert_valid_coolant

end module coolants_module
