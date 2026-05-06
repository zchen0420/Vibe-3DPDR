module healpix_module

  !T.Bisbas
  use definitions
  use healpix_types
  integer(kind=i4b):: ns_max ! ..
  integer(kind=i4b), dimension(0:1023) :: pix2x ! ..
  integer(kind=i4b), dimension(0:1023) :: pix2y ! ..
  real(kind=dp) :: x2pix(0:1023) ! ..
  real(kind=dp) :: y2pix(0:1023) ! ..
  real(kind=dp), allocatable :: vertex(:,:) ! ..
  real(kind=dp), allocatable :: vector(:) ! ..
end module healpix_module
