module chemistry_network_module
  use definitions, only : dp
  use healpix_types, only : i4b
  implicit none

  type :: reaction_network
    character(len=10), pointer :: species(:) => null()
    real(kind=dp), pointer :: mass(:) => null()
    real(kind=dp), pointer :: initial_abundance(:) => null()
    character(len=10), pointer :: reactant(:,:) => null()
    character(len=10), pointer :: product(:,:) => null()
    real(kind=dp), pointer :: alpha(:) => null()
    real(kind=dp), pointer :: beta(:) => null()
    real(kind=dp), pointer :: gamma(:) => null()
    real(kind=dp), pointer :: rtmin(:) => null()
    real(kind=dp), pointer :: rtmax(:) => null()
    integer(kind=i4b), pointer :: duplicate(:) => null()
  end type reaction_network

  type :: chemistry_state
    type(reaction_network) :: network
    real(kind=dp), pointer :: rate(:) => null()
  end type chemistry_state

contains

  subroutine allocate_chemistry_state(chemistry, nspec, nreac)
    type(chemistry_state), intent(inout) :: chemistry
    integer(kind=i4b), intent(in) :: nspec
    integer(kind=i4b), intent(in) :: nreac

    allocate(chemistry%network%species(1:nspec))
    allocate(chemistry%network%initial_abundance(1:nspec))
    allocate(chemistry%network%mass(1:nspec))
    allocate(chemistry%network%reactant(1:nreac,1:3))
    allocate(chemistry%network%product(1:nreac,1:4))
    allocate(chemistry%rate(1:nreac))
    allocate(chemistry%network%alpha(1:nreac))
    allocate(chemistry%network%beta(1:nreac))
    allocate(chemistry%network%gamma(1:nreac))
    allocate(chemistry%network%rtmin(1:nreac))
    allocate(chemistry%network%rtmax(1:nreac))
    allocate(chemistry%network%duplicate(1:nreac))
  end subroutine allocate_chemistry_state

end module chemistry_network_module
