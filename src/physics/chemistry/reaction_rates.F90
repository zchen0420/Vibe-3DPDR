module reaction_rates_module
  use chemistry_network_module, only : reaction_network
  use definitions, only : dp
  use global_module
  use healpix_types, only : i4b
  use maincode_module, only : runtime
  use photo_rate_interfaces_module
  use reaction_rate_kernels_module
  implicit none

  private
  public :: reaction_rate_environment
  public :: reaction_rate_indices
  public :: calculate_reaction_rates

  type :: reaction_rate_environment
    real(kind=dp) :: gas_temperature
    real(kind=dp) :: dust_temperature
    real(kind=dp), pointer :: radiation_surface(:) => null()
    real(kind=dp), pointer :: visual_extinction(:) => null()
    real(kind=dp), pointer :: column_density(:,:) => null()
  end type reaction_rate_environment

  type :: reaction_rate_indices
    integer(kind=i4b) :: grain_surface
    integer(kind=i4b) :: h2_photodissociation
    integer(kind=i4b) :: hd_photodissociation
    integer(kind=i4b) :: co_photodissociation
    integer(kind=i4b) :: carbon_photoionization
    integer(kind=i4b) :: silicon_photoionization
  end type reaction_rate_indices

contains

  subroutine calculate_reaction_rates(environment, network, rate, indices)
    type(reaction_rate_environment), intent(in) :: environment
    type(reaction_network), intent(in) :: network
    real(kind=dp), intent(out) :: rate(:)
    type(reaction_rate_indices), intent(out) :: indices

    real(kind=dp) :: temperature
    real(kind=dp) :: dust_temperature
    real(kind=dp) :: phi_pah,flux,yield
    integer(kind=i4b) :: i,j,k
    integer(kind=i4b) :: nreac
    integer(kind=i4b) :: ray_index_lower
    integer(kind=i4b) :: ray_index_upper

    temperature = environment%gas_temperature
    dust_temperature = environment%dust_temperature
    nreac = size(rate)
    ray_index_lower = lbound(environment%radiation_surface,1)
    ray_index_upper = ubound(environment%radiation_surface,1)

    associate(reactant => network%reactant, product => network%product, alpha => network%alpha, &
          &beta => network%beta, gamma => network%gamma, rtmin => network%rtmin, rtmax => network%rtmax, &
          &duplicate => network%duplicate, rad_surface => environment%radiation_surface, &
          &av => environment%visual_extinction, column => environment%column_density)

      !     Initialize the RATE coefficients.
      rate=0.0d0

      !     Initialize the stored reaction numbers. If they are not assigned
      !     subsequently, any attempt to access that reaction will generate an
      !     error and the code will crash. This is a useful bug catch.
      indices%grain_surface=0
      indices%h2_photodissociation=0
      indices%hd_photodissociation=0
      indices%co_photodissociation=0
      indices%carbon_photoionization=0
      indices%silicon_photoionization=0

      do i=1,nreac
        !        Determine the type of reaction
        if(reactant(i,2).eq."PHOTON") goto 1
        if(reactant(i,2).eq."CRP   ") goto 2
        if(reactant(i,2).eq."CRPHOT") goto 3
        if(reactant(i,2).eq."FREEZE") goto 4
        if(reactant(i,2).eq."ELFRZE") goto 5
        if(reactant(i,2).eq."CRH   ") goto 6
        if(reactant(i,2).eq."PHOTD ") goto 7
        if(reactant(i,2).eq."THERM ") goto 8
        if(reactant(i,2)(1:1).eq."#") goto 9
        if(reactant(i,2).eq."XRAY  ") goto 110
        if(reactant(i,2).eq."XRSEC ") goto 120
        if(reactant(i,2).eq."XRLYA ") goto 130
        if(reactant(i,2).eq."XRPHOT") goto 140

        !-----------------------------------------------------------------------

        !     Thermal reactions:

        !     The RATE of H2 formation on grains is calculated separately
        !     by the function H2_FORMATION_RATE (see function for details)
        if((reactant(i,1).eq."H  " .and. reactant(i,2).eq."H  "   .and. &
            & (reactant(i,3).eq."   " .or.  reactant(i,3).eq."#  ")) .and. &
            &  (product(i,1).eq."H2 " .and. &
            &  (product(i,2).eq."   " .or.  product(i,2).eq."#  "))) then
#ifdef H2FORM
        rate(i)=h2_formation_rate(temperature,dust_temperature)
#else
        !            RATE(I)=3.0D-18*SQRT(TEMPERATURE)*EXP(-(TEMPERATURE/1.0D3))
        rate(i)=3.0d-18*sqrt(temperature)
#endif
        indices%grain_surface=i
        goto 10
      end if

      !     Rates for reactions involving PAHs are calculated according to the
      !     treatment of Wolfire et al. (2003, ApJ, 587, 278; 2008, ApJ, 680, 384)
      if(any(reactant(i,:).eq."PAH  ") .or. any(reactant(i,:).eq."PAH0 ") .or. &
          & any(reactant(i,:).eq."PAH+ ") .or. any(reactant(i,:).eq."PAH- ")) then
      phi_pah=0.4d0
      rate(i)=alpha(i)*(temperature/100.0d0)**beta(i)*phi_pah
      goto 10
    end if

    !C     Check for large negative GAMMA values that might cause discrepant
    !C     rates at low temperatures. Set these rates to zero when T < RTMIN.

    !CODE RESPONSINBLE FOR O + H+ --> O+ + H
    if(duplicate(i).eq.0) then
      if(gamma(i).lt.-200.0d0 .and. temperature.lt.rtmin(i)) then
        rate(i)=0.0d0
      else
        !if (i==76) then
        !RATE(i)=ALPHA(i)
        !write(6,*) '1'
        !else
        rate(i)=arrhenius_reaction_rate(alpha(i),beta(i),gamma(i),temperature,rtmin(i))
        !endif
      end if
    else if(duplicate(i).eq.1) then
      j=i
      do
        if(temperature.le.rtmax(j)) then
          if(gamma(j).lt.-200.0d0 .and. temperature.lt.rtmin(j)) then
            rate(j)=0.0d0
          else
            !if (i==76) then
            !RATE(j)=ALPHA(j)
            !write(6,*) '2'
            !else
            rate(j)=arrhenius_reaction_rate(alpha(j),beta(j),gamma(j),temperature,rtmin(j))
            !endif
          end if
          exit
        else if(duplicate(j+1).lt.duplicate(j)) then
          if(gamma(j).lt.-200.0d0 .and. temperature.lt.rtmin(j)) then
            rate(j)=0.0d0
          else
            !if (i==76) then
            !RATE(j)=ALPHA(j)
            !write(6,*) '3'
            !else
            rate(j)=arrhenius_reaction_rate(alpha(j),beta(j),gamma(j),temperature,rtmin(j))
            !endif
          end if
          exit
        else
          rate(j)=0.0d0
          j=j+1
        end if
      end do
    end if
    goto 10

    !C-----------------------------------------------------------------------

    !C     Photoreactions:

    !C     Store the reaction number for H2 photodissociation. The RATE itself
    !C     is calculated separately by the function H2PDRATE (within shield.f)
    1       if(reactant(i,1).eq."H2 " .and. reactant(i,3).eq."   ") then
      !C           Loop over all rays
      do k=ray_index_lower,ray_index_upper
        rate(i)=rate(i) + h2pdrate(alpha(i),rad_surface(k),av(k),column(k,species_idx%nh2))
      end do
      if(product(i,1).eq."H " .and. product(i,2).eq."H ") indices%h2_photodissociation=i
      goto 10
    end if

    !C     Store the reaction number for HD photodissociation. The RATE itself
    !C     is calculated separately by the function H2PDRATE (within shield.f)
    if(reactant(i,1).eq."HD " .and. reactant(i,3).eq."   ") then
      !C           Loop over all rays
      do k=ray_index_lower,ray_index_upper
        rate(i)=rate(i) + h2pdrate(alpha(i),rad_surface(k),av(k),column(k,species_idx%nhd))
      end do
      if(any(product(i,:).eq."H ") .and. any(product(i,:).eq."D ")) indices%hd_photodissociation=i
      goto 10
    end if

    !C     Store the reaction number for !CO photodissociation. The RATE itself
    !C     is calculated separately by the function !COPDRATE (within shield.f)
    if(reactant(i,1).eq."CO " .and. reactant(i,3).eq."   " .and. &
        & any(product(i,:).eq."C ") .and. any(product(i,:).eq."O ")) then
    !C           Loop over all rays
    do k=ray_index_lower,ray_index_upper
      rate(i)=rate(i) + copdrate(alpha(i),rad_surface(k),av(k),column(k,species_idx%nco),column(k,species_idx%nh2))
    end do
    indices%co_photodissociation=i
    goto 10

  end if

  !C     Store the reaction number for !CI photoionization. The RATE itself
  !C     is calculated separately by the function CIPDRATE (within shield.f)
  if(reactant(i,1).eq."C  " .and. reactant(i,3).eq."   " .and.&
      &    ((product(i,1).eq."C+ " .and. product(i,2).eq."e- ") .or.&
      &     (product(i,1).eq."e- " .and. product(i,2).eq."C+ "))) then
  !C           Loop over all rays
  do k=ray_index_lower,ray_index_upper
    rate(i)=rate(i) + cipdrate(alpha(i),rad_surface(k),av(k),gamma(i),column(k,species_idx%nc),column(k,species_idx%nh2),temperature)
  end do
  indices%carbon_photoionization=i
  goto 10

end if

!C     Store the reaction number for SI photoionization. The RATE itself
!C     is calculated separately by the function SIPDRATE (within shield.f)
if(reactant(i,1).eq."S  " .and. reactant(i,3).eq."   " .and.&
    &    ((product(i,1).eq."S+ " .and. product(i,2).eq."e- ") .or.&
    &     (product(i,1).eq."e- " .and. product(i,2).eq."S+ "))) then
!C           Loop over all rays
do k=ray_index_lower,ray_index_upper
  rate(i)=rate(i) + sipdrate(alpha(i),rad_surface(k),av(k),gamma(i),column(k,species_idx%ns))
end do
indices%silicon_photoionization=i
goto 10
end if

if(duplicate(i).eq.0) then
  !C           Loop over all rays
  do k=ray_index_lower,ray_index_upper
    rate(i)=rate(i) + alpha(i)*rad_surface(k)*exp(-(gamma(i)*av(k)))/2.0
  end do
else if(duplicate(i).eq.1) then
  j=i
  do
    if(temperature.le.rtmax(j)) then
      !C                 Loop over all rays
      do k=ray_index_lower,ray_index_upper
        rate(j)=rate(j) + alpha(j)*rad_surface(k)*exp(-(gamma(j)*av(k)))/2.0
      end do
      exit
    else if(duplicate(j+1).lt.duplicate(j)) then
      !C                 Loop over all rays
      do k=ray_index_lower,ray_index_upper
        rate(j)=rate(j) + alpha(j)*rad_surface(k)*exp(-(gamma(j)*av(k)))/2.0
      end do
      exit
    else
      rate(j)=0.0d0
      j=j+1
    end if
  end do
end if
goto 10

!C-----------------------------------------------------------------------

!C     Cosmic ray-induced ionization:

2       if(duplicate(i).eq.0) then
  rate(i)=alpha(i)*runtime%cosmic_ray_ionization_rate
else if(duplicate(i).eq.1) then
  j=i
  do
    if(temperature.le.rtmax(j)) then
      rate(j)=alpha(j)*runtime%cosmic_ray_ionization_rate
      exit
    else if(duplicate(j+1).lt.duplicate(j)) then
      rate(j)=alpha(j)*runtime%cosmic_ray_ionization_rate
      exit
    else
      rate(j)=0.0d0
      j=j+1
    end if
  end do
end if
goto 10

110     rate(i)=0.0d0; goto 10
120     rate(i)=0.0d0; goto 10
130     rate(i)=0.0d0; goto 10
140     rate(i)=0.0d0; goto 10

!C-----------------------------------------------------------------------

!C     Photoreactions due to cosmic ray-induced secondary photons:

3       if(duplicate(i).eq.0) then
  rate(i)=alpha(i)*runtime%cosmic_ray_ionization_rate*(temperature/300.0d0)**beta(i)&
      &           *gamma(i)/(1.0d0-omega)
else if(duplicate(i).eq.1) then
  j=i
  do
    if(temperature.le.rtmax(j)) then
      rate(j)=alpha(j)*runtime%cosmic_ray_ionization_rate*(temperature/300.0d0)**beta(j)&
          &                 *gamma(j)/(1.0d0-omega)
      exit
    else if(duplicate(j+1).lt.duplicate(j)) then
      rate(j)=alpha(j)*runtime%cosmic_ray_ionization_rate*(temperature/300.0d0)**beta(j)&
          &                 *gamma(j)/(1.0d0-omega)
      exit
    else
      rate(j)=0.0d0
      j=j+1
    end if
  end do
end if
goto 10

!C-----------------------------------------------------------------------

!C     Freeze-out of neutral chemistry%network%species:

4       rate(i)=freeze_out_reaction_rate(alpha(i),beta(i),gamma(i),temperature)
goto 10

!C-----------------------------------------------------------------------

!C     Freeze-out of singly charged positive ions:

5       rate(i)=freeze_out_reaction_rate(alpha(i),beta(i),gamma(i),temperature)
goto 10

!C-----------------------------------------------------------------------

!C     Desorption due to cosmic ray heating:

!CC     Treatment of Hasegawa & Herbst (1993, MNRAS, 261, 83, Equation 15)
!C 6       RATE(I)=ALPHA(I)*runtime%cosmic_ray_ionization_rate

!C     Treatment of Roberts et al. (2007, MNRAS, 382, 773, Equation 3)
6       rate(i)=cosmic_ray_desorption_rate(gamma(i),runtime%cosmic_ray_ionization_rate)
goto 10

!C-----------------------------------------------------------------------

!C     Photodesorption:

7       yield=photodesorption_yield(temperature)
!C         FLUX=1.0D8 ! Flux of FUV photons in the unattenuated Habing field (in photons cm^-2 s^-1)
flux=1.7d8 ! Flux of FUV photons in the unattenuated Draine field (in photons cm^-2 s^-1)
!C        Loop over all rays
do k=ray_index_lower,ray_index_upper
  rate(i)=rate(i) + flux*rad_surface(k)*exp(-(1.8d0*av(k)))*2.4d-22*yield
end do
goto 10

!C-----------------------------------------------------------------------

!C     Thermal desorption:

!C     Treatment of Hasegawa, Herbst & Leung (1992, ApJS, 82, 167, Equations 2 & 3)
8       rate(i)=thermal_desorption_rate(alpha(i),gamma(i),dust_temperature)
goto 10

!C-----------------------------------------------------------------------

!C     Grain mantle reactions:

9       rate(i)=alpha(i)
goto 10

!C-----------------------------------------------------------------------

!C     Check that the RATE is physical (0<RATE(I)<1) and produce an error
!C     message if not. Impose a lower cut-off on all RATE coefficients to
!C     prevent the problem becoming too stiff. Rates less than 1E-99 are
!C     set to zero. Grain-surface reactions and desorption mechanisms are
!C     allowed rates greater than 1.
10      call clamp_reaction_rate(rate(i),i,reactant(i,1))
!C     End of loop over rates
end do

return
end associate
end subroutine calculate_reaction_rates

end module reaction_rates_module
