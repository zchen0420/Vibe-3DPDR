!-----------------------------------------------------------------------------
!
!  Copyright (C) 1997-2005 Krzysztof M. Gorski, Eric Hivon,
!                          Benjamin D. Wandelt, Anthony J. Banday,
!                          Matthias Bartelmann, Hans K. Eriksen,
!                          Frode K. Hansen, Martin Reinecke
!
!
!  This file is part of HEALPix.
!
!  HEALPix is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  HEALPix is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with HEALPix; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
!
!  For more information about HEALPix see http://healpix.jpl.nasa.gov
!
!-----------------------------------------------------------------------------
module healpix_types
  ! This module sets the types used in the Fortran 90 modules
  ! of the HEALPIX distribution and follows the example of Numerical Recipes
  !
  ! Benjamin D. Wandelt October 1997
  ! Eric Hivon June 1998
  ! Eric Hivon Oct  2001, edited to be compatible with 'F' compiler
  ! Eric Hivon July 2002, addition of i8b, i2b, i1b
  !                       addition of max_i8b, max_i2b and max_i1b
  !            Jan 2005, explicit form of max_i1b because of ifc 8.1.021
  !            June 2005, redefine i8b as 16 digit integer because of Nec f90 compiler

  ! Include definitions module for universal definition of DP, PR, SP
  use definitions

  integer, parameter, public :: i8b = selected_int_kind(16)
  integer, parameter, public :: i4b = selected_int_kind(9)
  integer, parameter, public :: i2b = selected_int_kind(4)
  integer, parameter, public :: i1b = selected_int_kind(2)
  integer, parameter, public :: lgt = kind(.true.)
  integer, parameter, public :: spc = kind((1.0_sp, 1.0_sp))
  integer, parameter, public :: dpc = kind((1.0_dp, 1.0_dp))
  !
  integer(i8b),  parameter, public :: max_i8b = huge(1_i8b)
  integer,       parameter, public :: max_i4b = huge(1_i4b)
  integer,       parameter, public :: max_i2b = huge(1_i2b)
  integer,       parameter, public :: max_i1b = 127
  real(kind=sp), parameter, public :: max_sp  = huge(1.0_sp)
  real(kind=dp), parameter, public :: max_dp  = huge(1.0_dp)

  ! Numerical Constant (Double precision)
  real(kind=dp), parameter, public :: quartpi=0.785398163397448309615660845819875721049_dp
  real(kind=dp), parameter, public :: halfpi= 1.570796326794896619231321691639751442099_dp
  real(kind=dp), parameter, public :: pi    = 3.141592653589793238462643383279502884197_dp
  real(kind=dp), parameter, public :: twopi = 6.283185307179586476925286766559005768394_dp
  real(kind=dp), parameter, public :: fourpi=12.56637061435917295385057353311801153679_dp
  real(kind=dp), parameter, public :: sqrt2 = 1.41421356237309504880168872420969807856967_dp
  real(kind=dp), parameter, public :: sq4pi_inv = 0.2820947917738781434740397257803862929220_dp
  real(kind=dp), parameter, public :: twothird = 0.6666666666666666666666666666666666666666_dp

  real(kind=dp), parameter, public :: rad2deg = 180.0_dp / pi
  real(kind=dp), parameter, public :: deg2rad = pi / 180.0_dp
  real(kind=sp), parameter, public :: hpx_sbadval = -1.6375e30_sp
  real(kind=dp), parameter, public :: hpx_dbadval = -1.6375e30_dp

  ! Maximum length of filenames
  integer, parameter :: filenamelen = 1024


  ! ---- Normalisation and convention ----
  ! normalisation of spin weighted functions
  real(kind=dp), parameter, public ::  kvs = 1.0_dp ! 1.0 : CMBFAST (Healpix 1.2)
  ! sign of Q
  real(kind=dp), parameter, public :: sgq = -1.0_dp ! -1 : CMBFAST (Healpix 1.2)
  ! sign of spin weighted function !
  real(kind=dp), parameter, public :: sw1 = -1.0_dp ! -1 : Healpix 1.2, bug correction
  real(kind=dp), parameter, public :: ikvs = 1.0_dp / kvs ! inverse of KvS

  !parameters for 3DPDR
  real(kind=dp), parameter, public ::  kb = 1.38065040d-16 !Boltzmann constant cgs
  real(kind=dp), parameter, public ::  c  = 2.99792458d+10 !speed of light cgs
  real(kind=dp), parameter, public ::  mp = 1.67262164d-24 !proton mass cgs
  real(kind=dp), parameter, public ::  hp = 6.62606896d-27 !Planck's constant cgs
  real(kind=dp), parameter, public ::  hb = 1.05457163d-27 !Planck's constant / 2pi
  real(kind=dp), parameter, public ::  hk = 4.79923734d-11 !Planck's constant / Boltzmann constant
  real(kind=dp), parameter, public ::  na = 6.02214179d+23 !Avogadro's number
  real(kind=dp), parameter, public ::  au = 1.66053878d-24 !atomic mass unit
  real(kind=dp), parameter, public ::  mh = 1.67372346d-24 !hydrogen mass cgs
  real(kind=dp), parameter, public ::  me = 9.10938215d-28 !electron mass cgs
  real(kind=dp), parameter, public ::  ec = 4.80320427d-10 !elementary charge in esu
  real(kind=dp), parameter, public ::  pc = 3.08568025d+18 !pc in cm
  real(kind=dp), parameter, public ::  ev = 1.60217646d-12 !electron volt in erg


end module healpix_types
