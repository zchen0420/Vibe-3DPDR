!C***********************************************************************
!C     Read in the chemical reaction rates and the species, masses and
!C     initial abundances (if specified). The  rates and species files
!C     are assumed to have comma separated values (CSV) format. This is
!C     in line with the Rate05 formatting, removing the need for file-
!C     dependent FORMAT statements.
!C***********************************************************************
!
!C-----------------------------------------------------------------------
!C     Read in the species data, including initial fractional abundances
!C     (3rd column) and their masses (4th column). Check that the value
!C     of NSPEC agrees with the number of species in the file and produce
!C     an error message if not.
!C-----------------------------------------------------------------------
SUBROUTINE READ_SPECIES(NSPEC,SPECIES,ABUNDANCE,MASS)

  !T.Bell
  use definitions
  use healpix_types
  use global_module

  IMPLICIT NONE
  INTEGER(kind=i4b),intent(in) :: NSPEC
  real(kind=dp), intent(out) :: ABUNDANCE(1:nspec),MASS(1:nspec)
  CHARACTER(len=10), intent(out) :: SPECIES(1:nspec)

  INTEGER(kind=i4b) :: I,INDEX,SPECIESFILE

  SPECIESFILE = 3

  !C     Initialize the variables and read in the species data. Check that
  !C     the value of NSPEC agrees with the number of species in the file
  !C     and produce an error message if not.
  SPECIES="          "
  ABUNDANCE=0.0D0
  MASS=0.0D0

  !C     Initialize all the species index labels. If they are not assigned
  !C     subsequently, any attempt to access that species will generate an
  !C     error and the code will crash. This is a useful bug catch.
  species_idx%NH=0
  species_idx%ND=0
  species_idx%NH2=0
  species_idx%NHD=0
  species_idx%NH2x=0
  species_idx%NPROTON=0
  species_idx%NC=0
  species_idx%NCx=0
  species_idx%NO=0
  species_idx%NOx=0
  species_idx%NN=0
  species_idx%NNx=0
  species_idx%NS=0
  species_idx%NSx=0
  species_idx%NHe=0
  species_idx%NHEx=0
  species_idx%NNA=0
  species_idx%NNAx=0
  species_idx%NMG=0
  species_idx%NMGx=0
  species_idx%NSI=0
  species_idx%NSIx=0
  species_idx%NFE=0
  species_idx%NFEx=0
  species_idx%NCL=0
  species_idx%NCLx=0
  species_idx%NCA=0
  species_idx%NCAx=0
  species_idx%NCAxx=0
  species_idx%NCO=0
  species_idx%NCH=0
  species_idx%NCH2=0
  species_idx%NOH=0
  species_idx%NO2=0
  species_idx%NCS=0
  species_idx%NH2O=0
  NELECT=0
  species_idx%NH3x=0
  species_idx%NH3Ox=0
  species_idx%NHCOx=0
#ifdef REDUCED
  OPEN(SPECIESFILE,FILE="data/species_reduced.d",STATUS="OLD")
#endif
#ifdef FULL
  OPEN(SPECIESFILE,FILE="data/species_full.d",STATUS="OLD")
#endif
#ifdef MYNETWORK
  OPEN(SPECIESFILE,FILE="data/species_mynetwork.d",STATUS="OLD")
#endif
  REWIND(SPECIESFILE)
  DO I=1,NSPEC
    READ(SPECIESFILE,*,END=1) INDEX,SPECIES(I),ABUNDANCE(I),MASS(I)

    !C        Assign the various index labels to their correct species.
    IF(SPECIES(I).EQ."H         ") species_idx%NH      = I
    IF(SPECIES(I).EQ."D         ") species_idx%ND      = I
    IF(SPECIES(I).EQ."H2        ") species_idx%NH2     = I
    IF(SPECIES(I).EQ."HD        ") species_idx%NHD     = I
    IF(SPECIES(I).EQ."H2+       ") species_idx%NH2x    = I
    IF(SPECIES(I).EQ."H3+       ") species_idx%NH3x    = I
    IF(SPECIES(I).EQ."H+        ") species_idx%NPROTON = I
    IF(SPECIES(I).EQ."C         ") species_idx%NC      = I
    IF(SPECIES(I).EQ."C+        ") species_idx%NCx     = I
    IF(SPECIES(I).EQ."O         ") species_idx%NO      = I
    IF(SPECIES(I).EQ."O+        ") species_idx%NOx     = I
    IF(SPECIES(I).EQ."N         ") species_idx%NN      = I
    IF(SPECIES(I).EQ."N+        ") species_idx%NNx     = I
    IF(SPECIES(I).EQ."S         ") species_idx%NS      = I
    IF(SPECIES(I).EQ."S+        ") species_idx%NSx     = I
    IF(SPECIES(I).EQ."He        ") species_idx%NHe     = I
    IF(SPECIES(I).EQ."HE        ") species_idx%NHe     = I
    IF(SPECIES(I).EQ."He+       ") species_idx%NHEx    = I
    IF(SPECIES(I).EQ."HE+       ") species_idx%NHEx    = I
    IF(SPECIES(I).EQ."Na        ") species_idx%NNA     = I
    IF(SPECIES(I).EQ."NA        ") species_idx%NNA     = I
    IF(SPECIES(I).EQ."Na+       ") species_idx%NNAx    = I
    IF(SPECIES(I).EQ."NA+       ") species_idx%NNAx    = I
    IF(SPECIES(I).EQ."Mg        ") species_idx%NMG     = I
    IF(SPECIES(I).EQ."MG        ") species_idx%NMG     = I
    IF(SPECIES(I).EQ."Mg+       ") species_idx%NMGx    = I
    IF(SPECIES(I).EQ."MG+       ") species_idx%NMGx    = I
    IF(SPECIES(I).EQ."Si        ") species_idx%NSI     = I
    IF(SPECIES(I).EQ."SI        ") species_idx%NSI     = I
    IF(SPECIES(I).EQ."Si+       ") species_idx%NSIx    = I
    IF(SPECIES(I).EQ."SI+       ") species_idx%NSIx    = I
    IF(SPECIES(I).EQ."Fe        ") species_idx%NFE     = I
    IF(SPECIES(I).EQ."FE        ") species_idx%NFE     = I
    IF(SPECIES(I).EQ."Fe+       ") species_idx%NFEx    = I
    IF(SPECIES(I).EQ."FE+       ") species_idx%NFEx    = I
    IF(SPECIES(I).EQ."Cl        ") species_idx%NCL     = I
    IF(SPECIES(I).EQ."CL        ") species_idx%NCL     = I
    IF(SPECIES(I).EQ."Cl+       ") species_idx%NCLx    = I
    IF(SPECIES(I).EQ."CL+       ") species_idx%NCLx    = I
    IF(SPECIES(I).EQ."Ca        ") species_idx%NCA     = I
    IF(SPECIES(I).EQ."CA        ") species_idx%NCA     = I
    IF(SPECIES(I).EQ."Ca+       ") species_idx%NCAx    = I
    IF(SPECIES(I).EQ."CA+       ") species_idx%NCAx    = I
    IF(SPECIES(I).EQ."Ca++      ") species_idx%NCAxx   = I
    IF(SPECIES(I).EQ."CA++      ") species_idx%NCAxx   = I
    IF(SPECIES(I).EQ."CO        ") species_idx%NCO     = I
    IF(SPECIES(I).EQ."CH        ") species_idx%NCH     = I
    IF(SPECIES(I).EQ."CH2       ") species_idx%NCH2    = I
    IF(SPECIES(I).EQ."OH        ") species_idx%NOH     = I
    IF(SPECIES(I).EQ."O2        ") species_idx%NO2     = I
    IF(SPECIES(I).EQ."CS        ") species_idx%NCS     = I
    IF(SPECIES(I).EQ."H2O       ") species_idx%NH2O    = I
    IF(SPECIES(I).EQ."H3O+      ") species_idx%NH3Ox   = I
    IF(SPECIES(I).EQ."HCO+      ") species_idx%NHCOx   = I
    IF(SPECIES(I).EQ."e-        ") NELECT  = I
    IF(SPECIES(I).EQ."ELECTR    ") NELECT  = I
  ENDDO

  I=I-1
  READ(SPECIESFILE,*,END=1)
  I=I+1
  1    IF(I.NE.NSPEC) THEN
  write(6,*) 'ERROR! Number of species (NSPEC) does not match ',&
      &           'the number of entries in the species file'
  STOP
ENDIF

!C     Check that the final species in the file is e-. Print a warning
!C     message to screen and logfile if not.
IF(SPECIES(NSPEC).NE."e-") THEN
  write(6,*) 'WARNING! Last entry in species file is not e-'
  WRITE(10,*)'WARNING! Last entry in species file is not e-'
ENDIF


!C     Check that the total hydrogen nuclei abundance adds up to 1.
!C     If not, modify the abundance of H2 (only consider H, H+ & H2)
IF((ABUNDANCE(species_idx%NH)+ABUNDANCE(species_idx%NPROTON)+2.0D0*ABUNDANCE(species_idx%NH2)).NE.1.0D0) THEN
  ABUNDANCE(species_idx%NH2)=0.5D0*(1.0D0-ABUNDANCE(species_idx%NH)-ABUNDANCE(species_idx%NPROTON))
ENDIF

!C     Calculate the intial electron abundance, if not
!C     specified, as the sum of the metal ion abundances
IF(ABUNDANCE(NELECT).LE.0.0D0) THEN
  ABUNDANCE(NELECT)=0.0D0
  IF(species_idx%NCx.NE.0) ABUNDANCE(NELECT)=ABUNDANCE(NELECT)+ABUNDANCE(species_idx%NCx)
  IF(species_idx%NSx.NE.0) ABUNDANCE(NELECT)=ABUNDANCE(NELECT)+ABUNDANCE(species_idx%NSx)
  IF(species_idx%NNAx.NE.0) ABUNDANCE(NELECT)=ABUNDANCE(NELECT)+ABUNDANCE(species_idx%NNAx)
  IF(species_idx%NMGx.NE.0) ABUNDANCE(NELECT)=ABUNDANCE(NELECT)+ABUNDANCE(species_idx%NMGx)
  IF(species_idx%NSIx.NE.0) ABUNDANCE(NELECT)=ABUNDANCE(NELECT)+ABUNDANCE(species_idx%NSIx)
  IF(species_idx%NFEx.NE.0) ABUNDANCE(NELECT)=ABUNDANCE(NELECT)+ABUNDANCE(species_idx%NFEx)
  IF(species_idx%NCLx.NE.0) ABUNDANCE(NELECT)=ABUNDANCE(NELECT)+ABUNDANCE(species_idx%NCLx)
  IF(species_idx%NCAx.NE.0) ABUNDANCE(NELECT)=ABUNDANCE(NELECT)+ABUNDANCE(species_idx%NCAx)
  IF(species_idx%NCAxx.NE.0) ABUNDANCE(NELECT)=ABUNDANCE(NELECT)+2.0D0*ABUNDANCE(species_idx%NCAxx)
  IF(species_idx%NPROTON.NE.0) ABUNDANCE(NELECT)=ABUNDANCE(NELECT)+ABUNDANCE(species_idx%NPROTON)
ENDIF

CLOSE(SPECIESFILE)
RETURN
END SUBROUTINE
