module photo_rate_interfaces_module

  use definitions
  use healpix_types
  implicit none

  interface
    function h2pdrate(k0,g0,av,nh2)
      use definitions
      use healpix_types
      real(kind=dp) :: h2pdrate
      real(kind=dp), intent(in) :: k0, g0, av
      real(kind=dp), intent(in) :: nh2
      real(kind=dp) :: lambda, scatter, h2shield2
    end function h2pdrate

    function copdrate(k0,g0,av,nco,nh2)
      use definitions
      use healpix_types
      real(kind=dp) :: copdrate
      real(kind=dp), intent(in) :: k0, g0, av, nco
      real(kind=dp), intent(in) :: nh2
      real(kind=dp) :: lambda, lbar, coshield, scatter
    end function copdrate

    function cipdrate(k0,g0,av,kav,nci,nh2,tgas)
      use definitions
      use healpix_types
      real(kind=dp) :: cipdrate
      real(kind=dp), intent(in) :: k0,g0,av,kav,nci,tgas
      real(kind=dp), intent(in) :: nh2
      real(kind=dp) :: tauc
    end function cipdrate

    function sipdrate(k0,g0,av,kav,nsi)
      use definitions
      use healpix_types
      real(kind=dp) :: sipdrate
      real(kind=dp), intent(in) :: k0,g0,av,kav,nsi
      real(kind=dp) :: taus
    end function sipdrate

    function h2shield1(nh2,dopw,radw)
      use definitions
      use healpix_types
      real(kind=dp) :: h2shield1
      real(kind=dp), intent(in) :: nh2
      real(kind=dp), intent(in) ::dopw,radw
      real(kind=dp) :: fpara, fosc, taud, r, t, u, jd, jr
    end function h2shield1

    function h2shield2(nh2)
      use definitions
      use healpix_types
      use uclpdr_module, only : start, numh2, col_grid, sh2_grid, sh2_deriv
      real(kind=dp) :: h2shield2
      real(kind=dp), intent(in) :: nh2
    end function h2shield2

    function coshield(nco,nh2)
      use definitions
      use healpix_types
      use uclpdr_module, only : start, nco_grid, nh2_grid, sco_grid, sco_deriv
      real(kind=dp) :: coshield
      real(kind=dp) :: lognco, lognh2
      real(kind=dp), intent(in) :: nco, nh2
    end function coshield

    function scatter(av,lambda)
      use definitions
      use healpix_types
      real(kind=dp) :: scatter
      real(kind=dp), intent(in) :: av, lambda
      real(kind=dp), dimension(0:5), save :: a = (/&
          &1.000d0,2.006d0,-1.438d0,0.7364d0,-0.5076d0,-0.0592d0/)
      real(kind=dp), dimension(0:5), save :: k = (/&
          &0.7514d0,0.8490d0,1.013d0,1.282d0,2.005d0,5.832d0/)
      real(kind=dp) :: exponent, xlambda
    end function scatter

    function xlambda(lambda)
      use definitions
      use healpix_types
      use uclpdr_module, only : start, n_grid, l_grid, x_grid, x_deriv
      real(kind=dp) :: xlambda
      real(kind=dp), intent(in) :: lambda
    end function xlambda

    function lbar(nco,nh2)
      use definitions
      use healpix_types
      real(kind=dp) :: lbar
      real(kind=dp) :: u,w
      real(kind=dp) :: nco, nh2
    end function lbar

#ifdef H2FORM

    function h2_formation_rate(gas_temperature,grain_temperature) result(rate)
      use definitions
      use healpix_types
      implicit none
      real(kind=dp) :: rate
      real(kind=dp), intent(in) :: gas_temperature,grain_temperature
    end function h2_formation_rate

#endif

  end interface

end module photo_rate_interfaces_module
