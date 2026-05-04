MODULE healpix_module

  !T.Bisbas
  use definitions
  use healpix_types
  integer(kind=i4b):: ns_max                      ! ..
  integer(kind=i4b), dimension(0:1023) :: pix2x   ! ..
  integer(kind=i4b), dimension(0:1023) :: pix2y   ! ..
  real(kind=DP) :: x2pix(0:1023)          ! ..
  real(kind=DP) :: y2pix(0:1023)          ! ..
  real(kind=DP), allocatable :: vertex(:,:)      ! ..
  real(kind=DP), allocatable :: vector(:)      ! ..
END MODULE healpix_module
