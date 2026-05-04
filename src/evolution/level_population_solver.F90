module level_population_solver_module
  use definitions
  use maincode_module
  use global_module
  use convergence_module, only : set_lte_populations

  implicit none

contains

  subroutine set_initial_lte_populations
    do pp=1,pdr_ptot
      p=IDlist_pdr(pp)
      call set_lte_populations(p, gastemperature(pp), pdr(p)%rho)
#ifndef GUESS_TEMP
      if (pp.eq.1) then
        write(6,*) ''
        write(6,*) 'Z(CII)  = ',CII_Z_FUNCTION
        write(6,*) 'Z(CI)   = ',CI_Z_FUNCTION
        write(6,*) 'Z(OI)   = ',OI_Z_FUNCTION
        write(6,*) 'Z(C12O) = ',C12O_Z_FUNCTION
        write(6,*) ''
      endif
#endif
    enddo
  end subroutine set_initial_lte_populations

  subroutine refresh_lte_populations_for_unconverged_points
    integer(kind=i4b) :: point_index
    integer(kind=i4b) :: point_id

    do point_index=1,pdr_ptot
      point_id=IDlist_pdr(point_index)
#ifdef THERMALBALANCE
      if (level_converged(point_index).or.converged(point_index)) cycle
#else
      if (level_converged(point_index)) cycle
#endif
      call set_lte_populations(point_id, gastemperature(point_index), pdr(point_id)%rho)
    enddo
  end subroutine refresh_lte_populations_for_unconverged_points

  subroutine solve_level_populations_lvg
    integer(kind=i4b) :: point_index
    integer(kind=i4b) :: point_id

    integer(kind=i4b) :: ray_index
    integer(kind=i4b) :: eval_index
    integer(kind=i4b) :: level_index

    real(kind=dp), allocatable :: CII_C_COEFFS(:,:)
    real(kind=dp), allocatable :: CI_C_COEFFS(:,:)
    real(kind=dp), allocatable :: OI_C_COEFFS(:,:)
    real(kind=dp), allocatable :: C12O_C_COEFFS(:,:)

    real(kind=dp), allocatable :: transition_CII(:,:)
    real(kind=dp), allocatable :: transition_CI(:,:)
    real(kind=dp), allocatable :: transition_OI(:,:)
    real(kind=dp), allocatable :: transition_C12O(:,:)

    real(kind=dp), allocatable :: dummyarray_CII(:,:)
    real(kind=dp), allocatable :: dummyarray_CI(:,:)
    real(kind=dp), allocatable :: dummyarray_OI(:,:)
    real(kind=dp), allocatable :: dummyarray_C12O(:,:)

    real(kind=dp), allocatable :: dummyarray_CII_tau(:,:,:)
    real(kind=dp), allocatable :: dummyarray_CI_tau(:,:,:)
    real(kind=dp), allocatable :: dummyarray_OI_tau(:,:,:)
    real(kind=dp), allocatable :: dummyarray_C12O_tau(:,:,:)

    real(kind=dp), allocatable :: dummyarray_CII_beta(:,:,:)
    real(kind=dp), allocatable :: dummyarray_CI_beta(:,:,:)
    real(kind=dp), allocatable :: dummyarray_OI_beta(:,:,:)
    real(kind=dp), allocatable :: dummyarray_C12O_beta(:,:,:)

    real(kind=dp), allocatable :: CIIsolution(:)
    real(kind=dp), allocatable :: CIsolution(:)
    real(kind=dp), allocatable :: OIsolution(:)
    real(kind=dp), allocatable :: C12Osolution(:)

    real(kind=dp), allocatable :: CIIevalpop(:,:,:)
    real(kind=dp), allocatable :: CIevalpop(:,:,:)
    real(kind=dp), allocatable :: OIevalpop(:,:,:)
    real(kind=dp), allocatable :: C12Oevalpop(:,:,:)

    do point_index=1,pdr_ptot
      point_id=IDlist_pdr(point_index)

#ifdef THERMALBALANCE
      if (level_converged(point_index).or.converged(point_index)) cycle
#else
      if (level_converged(point_index)) cycle
#endif
      if (level_converged(point_index)) cycle

      pdr(point_id)%projected(:,0) = point_id

      call allocate_lvg_workspace
      call populate_evaluation_populations(point_id)

      call find_ccoeff(CII_NTEMP,CII_NLEV,gastemperature(point_index),CII_TEMPERATURES,&
          & CII_H,CII_HP,CII_EL,CII_HE,CII_H2,CII_PH2,CII_OH2,&
          & CII_C_COEFFS,pdr(point_id)%abundance(NH)*pdr(point_id)%rho,&
          & pdr(point_id)%abundance(NPROTON)*pdr(point_id)%rho,&
          & pdr(point_id)%abundance(NELECT)*pdr(point_id)%rho,&
          & pdr(point_id)%abundance(NHE)*pdr(point_id)%rho,&
          & pdr(point_id)%abundance(NH2)*pdr(point_id)%rho,1)
      call escape_probability(transition_CII, dusttemperature(point_index), nrays, CII_nlev, &
          &CII_A_COEFFS, CII_B_COEFFS, CII_C_COEFFS, &
          &CII_frequencies, CIIevalpop, maxpoints, &
          &gastemperature(point_index), v_turb, pdr(point_id)%epray, pdr(point_id)%CII_pop, &
          &pdr(point_id)%epoint, CII_weights,CII_cool(point_index),dummyarray_CII,&
          &dummyarray_CII_tau,1,pdr(point_id)%rho,metallicity,dummyarray_CII_beta)
      pdr(point_id)%CII_line=dummyarray_CII
      pdr(point_id)%CII_optdepth=dummyarray_CII_tau
      call solvlevpop(CII_nlev,transition_CII,pdr(point_id)%abundance(NCx)*pdr(point_id)%rho,CIIsolution)
      CII_solution(point_index,:)=CIIsolution
#ifdef CO_FIX
      if (levpop_iteration.ge.120) then
        CII_solution(point_index,:)=pdr(point_id)%CII_pop
      else if (levpop_iteration.ge.75) then
        CII_solution(point_index,:)=0.5*(CII_solution(point_index,:) + pdr(point_id)%CII_pop)
      endif
#endif

      call find_ccoeff(CI_NTEMP,CI_NLEV,gastemperature(point_index),CI_TEMPERATURES,&
          & CI_H,CI_HP,CI_EL,CI_HE,CI_H2,CI_PH2,CI_OH2,&
          & CI_C_COEFFS,pdr(point_id)%abundance(NH)*pdr(point_id)%rho,&
          & pdr(point_id)%abundance(NPROTON)*pdr(point_id)%rho,&
          & pdr(point_id)%abundance(NELECT)*pdr(point_id)%rho,&
          & pdr(point_id)%abundance(NHE)*pdr(point_id)%rho,&
          & pdr(point_id)%abundance(NH2)*pdr(point_id)%rho,2)
      call escape_probability(transition_CI, dusttemperature(point_index), nrays, CI_nlev, &
          &CI_A_COEFFS, CI_B_COEFFS, CI_C_COEFFS, &
          &CI_frequencies, CIevalpop, maxpoints, &
          &gastemperature(point_index), v_turb, pdr(point_id)%epray, pdr(point_id)%CI_pop, &
          &pdr(point_id)%epoint,CI_weights,CI_cool(point_index),dummyarray_CI,&
          &dummyarray_CI_tau,2,pdr(point_id)%rho,metallicity,dummyarray_CI_beta)
      pdr(point_id)%CI_line=dummyarray_CI
      pdr(point_id)%CI_optdepth=dummyarray_CI_tau
      call solvlevpop(CI_nlev,transition_CI,pdr(point_id)%abundance(NC)*pdr(point_id)%rho,CIsolution)
      CI_solution(point_index,:)=CIsolution
#ifdef CO_FIX
      if (levpop_iteration.ge.120) then
        CI_solution(point_index,:)=pdr(point_id)%CI_pop
      else if (levpop_iteration.ge.75) then
        CI_solution(point_index,:)=0.5*(CI_solution(point_index,:) + pdr(point_id)%CI_pop)
      endif
#endif

      call find_ccoeff(OI_NTEMP,OI_NLEV,gastemperature(point_index),OI_TEMPERATURES,&
          & OI_H,OI_HP,OI_EL,OI_HE,OI_H2,OI_PH2,OI_OH2,&
          & OI_C_COEFFS,pdr(point_id)%abundance(NH)*pdr(point_id)%rho,&
          & pdr(point_id)%abundance(NPROTON)*pdr(point_id)%rho,&
          & pdr(point_id)%abundance(NELECT)*pdr(point_id)%rho,&
          & pdr(point_id)%abundance(NHE)*pdr(point_id)%rho,&
          & pdr(point_id)%abundance(NH2)*pdr(point_id)%rho,3)
      call escape_probability(transition_OI, dusttemperature(point_index), nrays, OI_nlev, &
          &OI_A_COEFFS, OI_B_COEFFS, OI_C_COEFFS, &
          &OI_frequencies, OIevalpop, maxpoints, &
          &gastemperature(point_index), v_turb, pdr(point_id)%epray, pdr(point_id)%OI_pop, &
          &pdr(point_id)%epoint,OI_weights,OI_cool(point_index),dummyarray_OI,&
          &dummyarray_OI_tau,3,pdr(point_id)%rho,metallicity,dummyarray_OI_beta)
      pdr(point_id)%OI_line=dummyarray_OI
      pdr(point_id)%OI_optdepth=dummyarray_OI_tau
      call solvlevpop(OI_nlev,transition_OI,pdr(point_id)%abundance(NO)*pdr(point_id)%rho,OIsolution)
      OI_solution(point_index,:)=OIsolution
#ifdef CO_FIX
      if (levpop_iteration.ge.120) then
        OI_solution(point_index,:)=pdr(point_id)%OI_pop
      else if (levpop_iteration.ge.75) then
        OI_solution(point_index,:)=0.5*(OI_solution(point_index,:) + pdr(point_id)%OI_pop)
      endif
#endif

      call find_ccoeff(C12O_NTEMP,C12O_NLEV,gastemperature(point_index),C12O_TEMPERATURES,&
          & C12O_H,C12O_HP,C12O_EL,C12O_HE,C12O_H2,C12O_PH2,C12O_OH2,&
          & C12O_C_COEFFS,pdr(point_id)%abundance(NH)*pdr(point_id)%rho,&
          & pdr(point_id)%abundance(NPROTON)*pdr(point_id)%rho,&
          & pdr(point_id)%abundance(NELECT)*pdr(point_id)%rho,&
          & pdr(point_id)%abundance(NHE)*pdr(point_id)%rho,&
          & pdr(point_id)%abundance(NH2)*pdr(point_id)%rho,4)
      call escape_probability(transition_C12O, dusttemperature(point_index), nrays, C12O_nlev, &
          &C12O_A_COEFFS, C12O_B_COEFFS, C12O_C_COEFFS, &
          &C12O_frequencies, C12Oevalpop, maxpoints, &
          &gastemperature(point_index), v_turb, pdr(point_id)%epray, pdr(point_id)%C12O_pop, &
          &pdr(point_id)%epoint,C12O_weights,C12O_cool(point_index),dummyarray_C12O,&
          &dummyarray_C12O_tau,4,pdr(point_id)%rho,metallicity,dummyarray_C12O_beta)
      pdr(point_id)%C12O_line=dummyarray_C12O
      pdr(point_id)%C12O_optdepth=dummyarray_C12O_tau
      call solvlevpop(C12O_nlev,transition_C12O,pdr(point_id)%abundance(NCO)*pdr(point_id)%rho,C12Osolution)
      C12O_solution(point_index,:)=C12Osolution

#ifdef CO_FIX
      if (CII_percentage.eq.100.and.CI_percentage.eq.100.and.OI_percentage.eq.100) then
        if (levpop_iteration.ge.120) then
          C12O_solution(point_index,:)=pdr(point_id)%C12O_pop
        else if (levpop_iteration.ge.75) then
          C12O_solution(point_index,:)=0.5*(C12O_solution(point_index,:) + pdr(point_id)%C12O_pop)
        endif
      endif
#endif

      call deallocate_lvg_workspace
    enddo

    total_cooling_rate=CII_cool+CI_cool+OI_cool+C12O_cool

  contains

    subroutine allocate_lvg_workspace
      allocate(CII_C_COEFFS(1:CII_NLEV,1:CII_NLEV))
      allocate(CI_C_COEFFS(1:CI_NLEV,1:CI_NLEV))
      allocate(OI_C_COEFFS(1:OI_NLEV,1:OI_NLEV))
      allocate(C12O_C_COEFFS(1:C12O_NLEV,1:C12O_NLEV))

      allocate(transition_CII(1:CII_nlev,1:CII_nlev))
      allocate(transition_CI(1:CI_nlev,1:CI_nlev))
      allocate(transition_OI(1:OI_nlev,1:OI_nlev))
      allocate(transition_C12O(1:C12O_nlev,1:C12O_nlev))

      allocate(dummyarray_CII(1:CII_nlev,1:CII_nlev))
      allocate(dummyarray_CI(1:CI_nlev,1:CI_nlev))
      allocate(dummyarray_OI(1:OI_nlev,1:OI_nlev))
      allocate(dummyarray_C12O(1:C12O_nlev,1:C12O_nlev))

      allocate(dummyarray_CII_tau(1:CII_nlev,1:CII_nlev,0:nrays-1))
      allocate(dummyarray_CI_tau(1:CI_nlev,1:CI_nlev,0:nrays-1))
      allocate(dummyarray_OI_tau(1:OI_nlev,1:OI_nlev,0:nrays-1))
      allocate(dummyarray_C12O_tau(1:C12O_nlev,1:C12O_nlev,0:nrays-1))

      allocate(dummyarray_CII_beta(1:CII_nlev,1:CII_nlev,0:nrays-1))
      allocate(dummyarray_CI_beta(1:CI_nlev,1:CI_nlev,0:nrays-1))
      allocate(dummyarray_OI_beta(1:OI_nlev,1:OI_nlev,0:nrays-1))
      allocate(dummyarray_C12O_beta(1:C12O_nlev,1:C12O_nlev,0:nrays-1))

      allocate(CIIsolution(1:CII_nlev))
      allocate(CIsolution(1:CI_nlev))
      allocate(OIsolution(1:OI_nlev))
      allocate(C12Osolution(1:C12O_nlev))

      allocate(CIIevalpop(0:nrays-1,0:maxpoints,1:CII_nlev))
      allocate(CIevalpop(0:nrays-1,0:maxpoints,1:CI_nlev))
      allocate(OIevalpop(0:nrays-1,0:maxpoints,1:OI_nlev))
      allocate(C12Oevalpop(0:nrays-1,0:maxpoints,1:C12O_nlev))

      CIIevalpop=0.0D0
      CIevalpop=0.0D0
      OIevalpop=0.0D0
      C12Oevalpop=0.0D0
    end subroutine allocate_lvg_workspace

    subroutine populate_evaluation_populations(source_point_id)
      integer(kind=i4b), intent(in) :: source_point_id
      integer(kind=i4b) :: projected_point_id

      do ray_index=0,nrays-1
        do eval_index=0,pdr(source_point_id)%epray(ray_index)
          projected_point_id=int(pdr(source_point_id)%projected(ray_index,eval_index))

          do level_index=1,CII_nlev
            CIIevalpop(ray_index,eval_index,level_index)=pdr(projected_point_id)%CII_pop(level_index)
          enddo
          do level_index=1,CI_nlev
            CIevalpop(ray_index,eval_index,level_index)=pdr(projected_point_id)%CI_pop(level_index)
          enddo
          do level_index=1,OI_nlev
            OIevalpop(ray_index,eval_index,level_index)=pdr(projected_point_id)%OI_pop(level_index)
          enddo
          do level_index=1,C12O_nlev
            C12Oevalpop(ray_index,eval_index,level_index)=pdr(projected_point_id)%C12O_pop(level_index)
          enddo
        enddo
      enddo
    end subroutine populate_evaluation_populations

    subroutine deallocate_lvg_workspace
      deallocate(CII_C_COEFFS)
      deallocate(CI_C_COEFFS)
      deallocate(OI_C_COEFFS)
      deallocate(C12O_C_COEFFS)

      deallocate(transition_CII)
      deallocate(transition_CI)
      deallocate(transition_OI)
      deallocate(transition_C12O)

      deallocate(dummyarray_CII)
      deallocate(dummyarray_CI)
      deallocate(dummyarray_OI)
      deallocate(dummyarray_C12O)

      deallocate(dummyarray_CII_tau)
      deallocate(dummyarray_CI_tau)
      deallocate(dummyarray_OI_tau)
      deallocate(dummyarray_C12O_tau)

      deallocate(dummyarray_CII_beta)
      deallocate(dummyarray_CI_beta)
      deallocate(dummyarray_OI_beta)
      deallocate(dummyarray_C12O_beta)

      deallocate(CIIsolution)
      deallocate(CIsolution)
      deallocate(OIsolution)
      deallocate(C12Osolution)

      deallocate(CIIevalpop)
      deallocate(CIevalpop)
      deallocate(OIevalpop)
      deallocate(C12Oevalpop)
    end subroutine deallocate_lvg_workspace

  end subroutine solve_level_populations_lvg

end module level_population_solver_module
