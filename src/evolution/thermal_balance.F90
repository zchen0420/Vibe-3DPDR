module thermal_balance_module
  use definitions
  use maincode_module
  use global_module

  implicit none

contains

  subroutine calculate_heating_and_temperature_updates(dobinarychop, previouschange)
    logical, allocatable, intent(inout) :: dobinarychop(:)
    character(len=1), allocatable, intent(inout) :: previouschange(:)

    integer(kind=i4b) :: NRGR,NRH2,NRHD,NRCO,NRCI,NRSI
    integer(kind=i4b) :: point_index
    integer(kind=i4b) :: point_id
    real(kind=dp) :: point_heating(1:12)
    real(kind=dp) :: current_temperature

    do point_index=1,pdr_ptot
      point_id=IDlist_pdr(point_index)
#ifdef THERMALBALANCE
      if (converged(point_index)) cycle
#endif

      call calculate_reaction_rates(gastemperature(point_index),dusttemperature(point_index),&
          &nrays,pdr(point_id)%rad_surface(0:nrays-1),pdr(point_id)%AV(0:nrays-1),&
          &column(point_index)%columndens_point(0:nrays-1,1:nspec),&
          &nreac, reactant, product, alpha, beta, gamma, rate, rtmin, rtmax, duplicate, nspec,&
          &NRGR,NRH2,NRHD,NRCO,NRCI,NRSI)
      call calc_heating(pdr(point_id)%rho,gastemperature(point_index),dusttemperature(point_index),&
          &pdr(point_id)%UVfield,v_turb,nspec,pdr(point_id)%abundance(:),nreac,rate,point_heating,&
          &NRGR,NRH2,NRHD,NRCO,NRCI,NRSI)

      all_heating(point_index,:)=point_heating

#ifdef THERMALBALANCE
      Fmean(point_index) = all_heating(point_index,12) - total_cooling_rate(point_index)
      Fratio(point_index) = 2.0D0*abs(Fmean(point_index))/ &
          &abs(all_heating(point_index,12) + total_cooling_rate(point_index))

      if (level_conv.and.first_time) then
        current_temperature = gastemperature(point_index)

        if (Fmean(point_index).eq.0) then
          Tlow(point_index) = gastemperature(point_index)
          Thigh(point_index) = gastemperature(point_index)
        else if (Fmean(point_index).gt.0) then
          Tlow(point_index) = gastemperature(point_index)
          gastemperature(point_index) = 1.3D0*gastemperature(point_index)
          previouschange(point_index) = "H"
        else if (Fmean(point_index).lt.0) then
          Thigh(point_index) = gastemperature(point_index)
          gastemperature(point_index) = 0.7D0*gastemperature(point_index)
          if (gastemperature(point_index).lt.Tmin) gastemperature(point_index)=Tmin
          previouschange(point_index) = "C"
        endif

        previousgastemperature(point_index) = current_temperature
        if (previousgastemperature(point_index).lt.Tmin) previousgastemperature(point_index)=Tmin

      else if (level_conv.and..not.first_time) then
        current_temperature = gastemperature(point_index)

        if (Fratio(point_index).le.Fcrit) converged(point_index) = .true.

        if (.not.dobinarychop(point_index)) then
          if (Fmean(point_index).gt.0.and.previouschange(point_index).eq."H") then
            Tlow(point_index) = gastemperature(point_index)
            gastemperature(point_index) = 1.3D0*gastemperature(point_index)
            Thigh(point_index) = gastemperature(point_index)
            previouschange(point_index) = "H"
          endif

          if (Fmean(point_index).lt.0.and.previouschange(point_index).eq."C") then
            Thigh(point_index) = gastemperature(point_index)
            gastemperature(point_index) = 0.7D0*gastemperature(point_index)
            Tlow(point_index) = gastemperature(point_index)
            previouschange(point_index) = "C"
            if (gastemperature(point_index).lt.Tmin) then
              gastemperature(point_index)=Tmin
              Tlow(point_index)=Tmin
              Thigh(point_index)=Tmin
            endif
          endif

          if (Fmean(point_index).gt.0.and.previouschange(point_index).eq."C") then
            gastemperature(point_index) = (Thigh(point_index) + Tlow(point_index))/2.0D0
            dobinarychop(point_index)=.true.
          endif

          if (Fmean(point_index).lt.0.and.previouschange(point_index).eq."H") then
            gastemperature(point_index) = (Thigh(point_index) + Tlow(point_index))/2.0D0
            dobinarychop(point_index)=.true.
          endif

        else
          if (Fmean(point_index).gt.0) then
            Tlow(point_index) = gastemperature(point_index)
            gastemperature(point_index) = (gastemperature(point_index) + Thigh(point_index)) / 2.0D0
          endif
          if (Fmean(point_index).lt.0) then
            Thigh(point_index) = gastemperature(point_index)
            gastemperature(point_index) = (gastemperature(point_index) + Tlow(point_index)) / 2.0D0
          endif
        endif

#ifdef TEMP_FIX
        if ((abs(gastemperature(point_index)-previousgastemperature(point_index)).le.Tdiff).and. &
            &(Fratio(point_index).gt.Fcrit)) then
          converged(point_index)=.true.
        endif
#endif

        previousgastemperature(point_index) = current_temperature

        if ((current_temperature.lt.Tmin).and.(Fmean(point_index).lt.0)) converged(point_index)=.true.
        if ((current_temperature.gt.Tmax).and.(Fmean(point_index).gt.0)) converged(point_index)=.true.

        if (converged(point_index)) then
          if (current_temperature.lt.Tmin) then
            previousgastemperature(point_index) = Tmin
            gastemperature(point_index) = Tmin
            if (doleveltmin(point_index)) then
              converged(point_index)=.true.
            else
              converged(point_index)=.false.
              level_converged(point_index)=.false.
              current_temperature=Tmin
            endif
            doleveltmin(point_index)=.true.
          endif

          if (current_temperature.gt.Tmax) then
            previousgastemperature(point_index) = Tmax
            gastemperature(point_index) = Tmax
          endif
        endif
      endif
#endif
    enddo
  end subroutine calculate_heating_and_temperature_updates

end module thermal_balance_module
