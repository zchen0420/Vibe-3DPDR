!=======================================================================
!
!  Calculate the dust temperature for each particle using the treatment
!  of Hollenbach, Takahashi & Tielens (1991, ApJ, 377, 192, eqns 5 & 6)
!  for the heating due to the incident FUV photons and the treatment of
!  Meijerink & Spaans (2005, A&A, 436, 397, eqn B.6) for heating due to
!  the incident flux of X-ray photons.
!
!  Among other things, the dust temperature can influence:
!
!     1) Cooling budget by emitting FIR photons that
!        interact with the line radiative transfer;
!     2) Gas-grain collisional heating or cooling rate;
!     3) H2 formation by changing the sticking probability;
!     4) Evaporation and condensation of molecules on grains.
!
!  The formula derived by Hollenbach, Takahashi & Tielens (1991) has
!  been modified to include the attenuation of the IR radiation. The
!  incident FUV radiation is absorbed and re-emitted in the infrared
!  by dust at the surface of the cloud (up to Av ~ 1mag). In the HTT
!  derivation, this IR radiation then serves as a second heat source
!  for dust deeper into the cloud. However, in their treatment, this
!  second re-radiated component is not attenuated with distance into
!  the cloud so it is *undiluted* with depth, leading to higher dust
!  temperatures deep within the cloud which in turn heat the gas via
!  collisions to unrealistically high temperatures. Models with high
!  gas densities and high incident FUV fluxes (e.g. n_H = 10^5 cm-3,
!  X_0 = 10^8 Draine) can produce T_gas ~ 100 K at Av ~ 50 mag!
!
!  Attenuation of the FIR radiation has therefore been introduced by
!  using an approximation for the infrared-only dust temperature from
!  Rowan-Robinson (1980, eqn 30b):
!
!  T_dust = T_0*(r/r_0)^(-0.4)
!
!  where r_0 is the cloud depth at which T_dust = T_0, corresponding
!  to an A_V of ~ 1 mag, the assumed size of the outer region of the
!  cloud that processes the incident FUV radiation and then re-emits
!  it in the FIR (see the original HTT 1991 paper for details). This
!  should prevent the dust temperature from dropping off too rapidly
!  with distance and maintain a larger warm dust region (~50-100 K).
!
!-----------------------------------------------------------------------
subroutine calculate_dust_temperatures

  use healpix_types
  use maincode_module

  implicit none

  integer(kind=i4b) :: j, point_id, point_index

  real(kind=dp) :: nu_0,r_0,t_0,tau_100
  real(kind=dp) :: t_cmb

  !  Parameters used in the HHT equations (see their paper for details)
  nu_0=2.65d15
  tau_100=1.0d-3
  r_0=1.0d0/runtime%av_scale
  t_cmb=2.73d0

  do point_index=1,pdr_ptot ! Loop over particles
    point_id=grid%pdr_ids(point_index)
    !     Calculate the contribution to the dust temperature from the local FUV flux and the CMB background
    grid%points(point_id)%dust_temperature=8.9d-11*nu_0*(1.71d0*grid%points(point_id)%uvfield)+t_cmb**5

    do j=0,nrays-1 ! Loop over rays

      !        The minimum dust temperature is related to the incident FUV flux along each ray
      !        Convert the incident FUV flux from Draine to Habing units by multiplying by 1.7
      t_0=12.2*(1.71d0*grid%points(point_id)%rad_surface(j))**0.2

      !        Add the contribution to the dust temperature from the FUV flux incident along this ray
      if(t_0.gt.0) grid%points(point_id)%dust_temperature=grid%points(point_id)%dust_temperature &
          & + (0.42-log(3.45d-2*tau_100*t_0))*(3.45d-2*tau_100*t_0)*t_0**5

    end do ! End of loop over rays

    !     Convert from total dust emission intensity to dust temperature
    grid%points(point_id)%dust_temperature=grid%points(point_id)%dust_temperature**0.2

    !     Impose a lower limit on the dust temperature, since values below 10 K can dramatically
    !     limit the rate of H2 formation on grains (the molecule cannot desorb from the surface)
    if(grid%points(point_id)%dust_temperature.lt.10.0d0) then
      grid%points(point_id)%dust_temperature=10.0d0
    end if

    !     Check that the dust temperature is physical
    if(grid%points(point_id)%dust_temperature.gt.1000) then
      write(6,*) 'ERROR! Calculated dust temperature exceeds 1000 K'
      stop
    end if

  end do ! End of loop over particles

  return
end subroutine calculate_dust_temperatures
!=======================================================================
