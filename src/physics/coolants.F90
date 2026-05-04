module coolants_module
  use healpix_types, only : i4b
  implicit none

  integer(kind=i4b), parameter :: COOLANT_COUNT = 4
  integer(kind=i4b), parameter :: COOLANT_CII = 1
  integer(kind=i4b), parameter :: COOLANT_CI = 2
  integer(kind=i4b), parameter :: COOLANT_OI = 3
  integer(kind=i4b), parameter :: COOLANT_C12O = 4

  character(len=4), parameter :: COOLANT_LABELS(COOLANT_COUNT) = (/&
      &'CII ', 'CI  ', 'OI  ', 'CO  '/)
  integer(kind=i4b), parameter :: COOLANT_DEFAULT_NLEVELS_TABLE(COOLANT_COUNT) = (/5, 5, 5, 41/)
  integer(kind=i4b), parameter :: COOLANT_DEFAULT_NTEMPS_TABLE(COOLANT_COUNT) = (/18, 29, 27, 25/)

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
