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
subroutine read_species(nspec,species,abundance,mass)

  !T.Bell
  use definitions
  use healpix_types
  use global_module

  implicit none
  integer(kind=i4b),intent(in) :: nspec
  real(kind=dp), intent(out) :: abundance(1:nspec),mass(1:nspec)
  character(len=10), intent(out) :: species(1:nspec)

  integer(kind=i4b) :: i,index,speciesfile

  speciesfile = 3

  !C     Initialize the variables and read in the species data. Check that
  !C     the value of NSPEC agrees with the number of species in the file
  !C     and produce an error message if not.
  species="          "
  abundance=0.0d0
  mass=0.0d0

  !C     Initialize all the species index labels. If they are not assigned
  !C     subsequently, any attempt to access that species will generate an
  !C     error and the code will crash. This is a useful bug catch.
  species_idx%nh=0
  species_idx%nd=0
  species_idx%nh2=0
  species_idx%nhd=0
  species_idx%nh2x=0
  species_idx%nproton=0
  species_idx%nc=0
  species_idx%ncx=0
  species_idx%no=0
  species_idx%nox=0
  species_idx%nn=0
  species_idx%nnx=0
  species_idx%ns=0
  species_idx%nsx=0
  species_idx%nhe=0
  species_idx%nhex=0
  species_idx%nna=0
  species_idx%nnax=0
  species_idx%nmg=0
  species_idx%nmgx=0
  species_idx%nsi=0
  species_idx%nsix=0
  species_idx%nfe=0
  species_idx%nfex=0
  species_idx%ncl=0
  species_idx%nclx=0
  species_idx%nca=0
  species_idx%ncax=0
  species_idx%ncaxx=0
  species_idx%nco=0
  species_idx%nch=0
  species_idx%nch2=0
  species_idx%noh=0
  species_idx%no2=0
  species_idx%ncs=0
  species_idx%nh2o=0
  nelect=0
  species_idx%nh3x=0
  species_idx%nh3ox=0
  species_idx%nhcox=0
#ifdef REDUCED
  open(speciesfile,file="data/species_reduced.d",status="OLD")
#endif
#ifdef FULL
  open(speciesfile,file="data/species_full.d",status="OLD")
#endif
#ifdef MYNETWORK
  open(speciesfile,file="data/species_mynetwork.d",status="OLD")
#endif
  rewind(speciesfile)
  do i=1,nspec
    read(speciesfile,*,end=1) index,species(i),abundance(i),mass(i)

    !C        Assign the various index labels to their correct species.
    if(species(i).eq."H         ") species_idx%nh      = i
    if(species(i).eq."D         ") species_idx%nd      = i
    if(species(i).eq."H2        ") species_idx%nh2     = i
    if(species(i).eq."HD        ") species_idx%nhd     = i
    if(species(i).eq."H2+       ") species_idx%nh2x    = i
    if(species(i).eq."H3+       ") species_idx%nh3x    = i
    if(species(i).eq."H+        ") species_idx%nproton = i
    if(species(i).eq."C         ") species_idx%nc      = i
    if(species(i).eq."C+        ") species_idx%ncx     = i
    if(species(i).eq."O         ") species_idx%no      = i
    if(species(i).eq."O+        ") species_idx%nox     = i
    if(species(i).eq."N         ") species_idx%nn      = i
    if(species(i).eq."N+        ") species_idx%nnx     = i
    if(species(i).eq."S         ") species_idx%ns      = i
    if(species(i).eq."S+        ") species_idx%nsx     = i
    if(species(i).eq."He        ") species_idx%nhe     = i
    if(species(i).eq."HE        ") species_idx%nhe     = i
    if(species(i).eq."He+       ") species_idx%nhex    = i
    if(species(i).eq."HE+       ") species_idx%nhex    = i
    if(species(i).eq."Na        ") species_idx%nna     = i
    if(species(i).eq."NA        ") species_idx%nna     = i
    if(species(i).eq."Na+       ") species_idx%nnax    = i
    if(species(i).eq."NA+       ") species_idx%nnax    = i
    if(species(i).eq."Mg        ") species_idx%nmg     = i
    if(species(i).eq."MG        ") species_idx%nmg     = i
    if(species(i).eq."Mg+       ") species_idx%nmgx    = i
    if(species(i).eq."MG+       ") species_idx%nmgx    = i
    if(species(i).eq."Si        ") species_idx%nsi     = i
    if(species(i).eq."SI        ") species_idx%nsi     = i
    if(species(i).eq."Si+       ") species_idx%nsix    = i
    if(species(i).eq."SI+       ") species_idx%nsix    = i
    if(species(i).eq."Fe        ") species_idx%nfe     = i
    if(species(i).eq."FE        ") species_idx%nfe     = i
    if(species(i).eq."Fe+       ") species_idx%nfex    = i
    if(species(i).eq."FE+       ") species_idx%nfex    = i
    if(species(i).eq."Cl        ") species_idx%ncl     = i
    if(species(i).eq."CL        ") species_idx%ncl     = i
    if(species(i).eq."Cl+       ") species_idx%nclx    = i
    if(species(i).eq."CL+       ") species_idx%nclx    = i
    if(species(i).eq."Ca        ") species_idx%nca     = i
    if(species(i).eq."CA        ") species_idx%nca     = i
    if(species(i).eq."Ca+       ") species_idx%ncax    = i
    if(species(i).eq."CA+       ") species_idx%ncax    = i
    if(species(i).eq."Ca++      ") species_idx%ncaxx   = i
    if(species(i).eq."CA++      ") species_idx%ncaxx   = i
    if(species(i).eq."CO        ") species_idx%nco     = i
    if(species(i).eq."CH        ") species_idx%nch     = i
    if(species(i).eq."CH2       ") species_idx%nch2    = i
    if(species(i).eq."OH        ") species_idx%noh     = i
    if(species(i).eq."O2        ") species_idx%no2     = i
    if(species(i).eq."CS        ") species_idx%ncs     = i
    if(species(i).eq."H2O       ") species_idx%nh2o    = i
    if(species(i).eq."H3O+      ") species_idx%nh3ox   = i
    if(species(i).eq."HCO+      ") species_idx%nhcox   = i
    if(species(i).eq."e-        ") nelect  = i
    if(species(i).eq."ELECTR    ") nelect  = i
  end do

  i=i-1
  read(speciesfile,*,end=1)
  i=i+1
  1    if(i.ne.nspec) then
    write(6,*) 'ERROR! Number of species (NSPEC) does not match ',&
        &           'the number of entries in the species file'
    stop
  end if

  !C     Check that the final species in the file is e-. Print a warning
  !C     message to screen and logfile if not.
  if(species(nspec).ne."e-") then
    write(6,*) 'WARNING! Last entry in species file is not e-'
    write(10,*)'WARNING! Last entry in species file is not e-'
  end if


  !C     Check that the total hydrogen nuclei abundance adds up to 1.
  !C     If not, modify the abundance of H2 (only consider H, H+ & H2)
  if((abundance(species_idx%nh)+abundance(species_idx%nproton)+2.0d0*abundance(species_idx%nh2)).ne.1.0d0) then
    abundance(species_idx%nh2)=0.5d0*(1.0d0-abundance(species_idx%nh)-abundance(species_idx%nproton))
  end if

  !C     Calculate the intial electron abundance, if not
  !C     specified, as the sum of the metal ion abundances
  if(abundance(nelect).le.0.0d0) then
    abundance(nelect)=0.0d0
    if(species_idx%ncx.ne.0) abundance(nelect)=abundance(nelect)+abundance(species_idx%ncx)
    if(species_idx%nsx.ne.0) abundance(nelect)=abundance(nelect)+abundance(species_idx%nsx)
    if(species_idx%nnax.ne.0) abundance(nelect)=abundance(nelect)+abundance(species_idx%nnax)
    if(species_idx%nmgx.ne.0) abundance(nelect)=abundance(nelect)+abundance(species_idx%nmgx)
    if(species_idx%nsix.ne.0) abundance(nelect)=abundance(nelect)+abundance(species_idx%nsix)
    if(species_idx%nfex.ne.0) abundance(nelect)=abundance(nelect)+abundance(species_idx%nfex)
    if(species_idx%nclx.ne.0) abundance(nelect)=abundance(nelect)+abundance(species_idx%nclx)
    if(species_idx%ncax.ne.0) abundance(nelect)=abundance(nelect)+abundance(species_idx%ncax)
    if(species_idx%ncaxx.ne.0) abundance(nelect)=abundance(nelect)+2.0d0*abundance(species_idx%ncaxx)
    if(species_idx%nproton.ne.0) abundance(nelect)=abundance(nelect)+abundance(species_idx%nproton)
  end if

  close(speciesfile)
  return
end subroutine
