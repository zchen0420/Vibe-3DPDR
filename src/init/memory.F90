module memory_module
  use definitions, only : dp
  use healpix_types, only : i4b
  use coolants_module, only : COOLANT_CII, COOLANT_CI, COOLANT_OI, COOLANT_C12O

contains

  subroutine allocations

    use maincode_module
    use uclpdr_module, only : SCO_GRID
    use global_module, only : all_heating
    implicit none

    !load SCO_GRID data [UCL_PDR]
    SCO_GRID(1:8,1) = (/0.000D+00,-1.408D-02,-1.099D-01,-4.400D-01,-1.154D+00,-1.888D+00,-2.760D+00,-4.001D+00/)
    SCO_GRID(1:8,2) = (/-8.539D-02,-1.015D-01,-2.104D-01,-5.608D-01,-1.272D+00,-1.973D+00,-2.818D+00,-4.055D+00/)
    SCO_GRID(1:8,3) = (/-1.451D-01,-1.612D-01,-2.708D-01,-6.273D-01,-1.355D+00,-2.057D+00,-2.902D+00,-4.122D+00/)
    SCO_GRID(1:8,4) = (/-4.559D-01,-4.666D-01,-5.432D-01,-8.665D-01,-1.602D+00,-2.303D+00,-3.146D+00,-4.421D+00/)
    SCO_GRID(1:8,5) = (/-1.303D+00,-1.312D+00,-1.367D+00,-1.676D+00,-2.305D+00,-3.034D+00,-3.758D+00,-5.077D+00/)
    SCO_GRID(1:8,6) = (/-3.883D+00,-3.888D+00,-3.936D+00,-4.197D+00,-4.739D+00,-5.165D+00,-5.441D+00,-6.446D+00/)

    call allocate_coolant_storage(COOLANT_CII)
    call allocate_coolant_storage(COOLANT_CI)
    call allocate_coolant_storage(COOLANT_OI)
    call allocate_coolant_storage(COOLANT_C12O)

    allocate(C_COEFFS(1:NLEV,1:NLEV))

    allocate(species(1:nspec))
    allocate(dummyabundance(1:nspec))
    allocate(mass(1:nspec))
    allocate(reactant(1:nreac,1:3))
    allocate(product(1:nreac,1:4))
    allocate(rate(1:nreac))
    allocate(alpha(1:nreac))
    allocate(beta(1:nreac))
    allocate(gamma(1:nreac))
    allocate(rtmin(1:nreac))
    allocate(rtmax(1:nreac))
    allocate(duplicate(1:nreac))

    !allocations start from 0 to cope with the ONE dark molecular element
    allocate(total_cooling_rate(0:pdr_ptot))
    allocate(CII_cool(0:pdr_ptot))
    allocate(CI_cool(0:pdr_ptot))
    allocate(OI_cool(0:pdr_ptot))
    allocate(C12O_cool(0:pdr_ptot))

    allocate(all_heating(0:pdr_ptot,1:12))

    allocate(dusttemperature(0:pdr_ptot))
    allocate(gastemperature(0:pdr_ptot))
    allocate(previousgastemperature(0:pdr_ptot))

#ifdef THERMALBALANCE
    allocate(Fratio(0:pdr_ptot))
    allocate(Fmean(0:pdr_ptot));Fmean=0
    allocate(Tlow(0:pdr_ptot))
    allocate(Thigh(0:pdr_ptot))
#endif

    return
  end subroutine allocations

  subroutine allocate_coolant_storage(coolant_id)
    use maincode_module
    integer(kind=i4b), intent(in) :: coolant_id

    select case (coolant_id)
    case (COOLANT_CII)
      call allocate_coolant_arrays(CII_NLEV,CII_NTEMP,CII_ENERGIES,CII_WEIGHTS,&
          &CII_A_COEFFS,CII_B_COEFFS,CII_FREQUENCIES,CII_TEMPERATURES,&
          &CII_HP,CII_H,CII_EL,CII_HE,CII_H2,CII_PH2,CII_OH2)
    case (COOLANT_CI)
      call allocate_coolant_arrays(CI_NLEV,CI_NTEMP,CI_ENERGIES,CI_WEIGHTS,&
          &CI_A_COEFFS,CI_B_COEFFS,CI_FREQUENCIES,CI_TEMPERATURES,&
          &CI_HP,CI_H,CI_EL,CI_HE,CI_H2,CI_PH2,CI_OH2)
    case (COOLANT_OI)
      call allocate_coolant_arrays(OI_NLEV,OI_NTEMP,OI_ENERGIES,OI_WEIGHTS,&
          &OI_A_COEFFS,OI_B_COEFFS,OI_FREQUENCIES,OI_TEMPERATURES,&
          &OI_HP,OI_H,OI_EL,OI_HE,OI_H2,OI_PH2,OI_OH2)
    case (COOLANT_C12O)
      call allocate_coolant_arrays(C12O_NLEV,C12O_NTEMP,C12O_ENERGIES,C12O_WEIGHTS,&
          &C12O_A_COEFFS,C12O_B_COEFFS,C12O_FREQUENCIES,C12O_TEMPERATURES,&
          &C12O_HP,C12O_H,C12O_EL,C12O_HE,C12O_H2,C12O_PH2,C12O_OH2)
    case default
      stop 'Invalid coolant storage id'
    end select
  end subroutine allocate_coolant_storage

  subroutine allocate_coolant_arrays(nlevels, ntemps, energies, weights, a_coeffs, b_coeffs, &
        &frequencies, temperatures, hp, h, el, he, h2, ph2, oh2)
    integer(kind=i4b), intent(in) :: nlevels
    integer(kind=i4b), intent(in) :: ntemps
    real(kind=dp), allocatable, intent(inout) :: energies(:)
    real(kind=dp), allocatable, intent(inout) :: weights(:)
    real(kind=dp), allocatable, intent(inout) :: a_coeffs(:,:)
    real(kind=dp), allocatable, intent(inout) :: b_coeffs(:,:)
    real(kind=dp), allocatable, intent(inout) :: frequencies(:,:)
    real(kind=dp), allocatable, intent(inout) :: temperatures(:,:)
    real(kind=dp), allocatable, intent(inout) :: hp(:,:,:)
    real(kind=dp), allocatable, intent(inout) :: h(:,:,:)
    real(kind=dp), allocatable, intent(inout) :: el(:,:,:)
    real(kind=dp), allocatable, intent(inout) :: he(:,:,:)
    real(kind=dp), allocatable, intent(inout) :: h2(:,:,:)
    real(kind=dp), allocatable, intent(inout) :: ph2(:,:,:)
    real(kind=dp), allocatable, intent(inout) :: oh2(:,:,:)

    allocate(energies(1:nlevels))
    allocate(weights(1:nlevels))
    allocate(a_coeffs(1:nlevels,1:nlevels))
    allocate(b_coeffs(1:nlevels,1:nlevels))
    allocate(frequencies(1:nlevels,1:nlevels))
    allocate(temperatures(1:7,1:ntemps))
    allocate(hp(1:nlevels,1:nlevels,1:ntemps))
    allocate(h(1:nlevels,1:nlevels,1:ntemps))
    allocate(el(1:nlevels,1:nlevels,1:ntemps))
    allocate(he(1:nlevels,1:nlevels,1:ntemps))
    allocate(h2(1:nlevels,1:nlevels,1:ntemps))
    allocate(ph2(1:nlevels,1:nlevels,1:ntemps))
    allocate(oh2(1:nlevels,1:nlevels,1:ntemps))
  end subroutine allocate_coolant_arrays

end module memory_module
