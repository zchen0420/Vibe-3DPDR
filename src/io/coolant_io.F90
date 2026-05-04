module coolant_io_module
  use maincode_module

contains

  subroutine load_coolant_data
    call readinput(C12Oinput,C12O_NLEV,C12O_NTEMP,C12O_ENERGIES,C12O_WEIGHTS,&
        &     C12O_A_COEFFS,C12O_B_COEFFS,C12O_FREQUENCIES,C12O_TEMPERATURES,&
        &     C12O_H,C12O_HP,C12O_EL,C12O_HE,C12O_H2,C12O_PH2,C12O_OH2)
    call readinput(CIIinput,CII_NLEV,CII_NTEMP,CII_ENERGIES,CII_WEIGHTS,&
        &     CII_A_COEFFS,CII_B_COEFFS,CII_FREQUENCIES,CII_TEMPERATURES,&
        &     CII_H,CII_HP,CII_EL,CII_HE,CII_H2,CII_PH2,CII_OH2)
    call readinput(CIinput,CI_NLEV,CI_NTEMP,CI_ENERGIES,CI_WEIGHTS,&
        &     CI_A_COEFFS,CI_B_COEFFS,CI_FREQUENCIES,CI_TEMPERATURES,&
        &     CI_H,CI_HP,CI_EL,CI_HE,CI_H2,CI_PH2,CI_OH2)
    call readinput(OIinput,OI_NLEV,OI_NTEMP,OI_ENERGIES,OI_WEIGHTS,&
        &     OI_A_COEFFS,OI_B_COEFFS,OI_FREQUENCIES,OI_TEMPERATURES,&
        &     OI_H,OI_HP,OI_EL,OI_HE,OI_H2,OI_PH2,OI_OH2)

    write(6,*) ''
  end subroutine load_coolant_data

end module coolant_io_module
