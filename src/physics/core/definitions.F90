! Definition of precision kind variables
! ============================================================================

!#include "macros.h"

! ============================================================================
!T.Bisbas -- taken from SEREN/HEALPix

module definitions

  integer, parameter :: dp = selected_real_kind(p=15) ! double precision
  integer, parameter :: sp = selected_real_kind(p=6) ! single precision

#ifdef DOUBLE_PRECISION
  integer, parameter :: pr = dp ! particle precision
#else
  integer, parameter :: pr = sp ! default = single
#endif

  integer, parameter :: ilp = 4 ! Integer long precision

end module definitions


! ============================================================================
