module global_module

  USE ISO_C_BINDING
  use definitions
  use healpix_types

  INTEGER(kind=i4b) :: NH,ND,NH2,NHD,NC,NCx,NCO,NO,NPROTON,NH2O,NHe, &
      &        NMG,NMGx,NN,NFE,NFEx,NSI,NSIx,NCA,NCAx,NCAxx,NS,NSx,NCS, &
      &        NOSH,NCL,NCLx,NH2x,NHEx,NOx,NNx,NNA,NNAx,NCH,NCH2,NOH,NO2, &
      &        NH3x, NH3Ox, NHCOx
  integer(kind=i4b),bind(c,name='global_module_mp_nelect_')::NELECT

  ! REAL(kind=dp), save :: ZETA=3.85D0,OMEGA=0.42D0,GRAIN_RADIUS=1.0D-5,METALLICITY=1.0D0
  REAL(kind=dp) :: g2d
  REAL(kind=dp) :: metallicity
  REAL(kind=dp) :: omega
  REAL(kind=dp) :: grain_radius
  ! REAL(kind=dp), save :: OMEGA=0.42D0,GRAIN_RADIUS=1.0D-7,METALLICITY=1.0D0
  ! REAL(kind=dp), save :: ZETA=1.0D0,OMEGA=0.42D0,GRAIN_RADIUS=1.0D-5,METALLICITY=1.0D0

  real(kind=dp),allocatable :: all_heating(:,:)

end module global_module
