module level_population_system_module
  use definitions, only : dp
  use healpix_types, only : i4b
  implicit none

contains

  subroutine solve_statistical_equilibrium(transition, density, solution)
    real(kind=dp), intent(in) :: transition(:,:)
    real(kind=dp), intent(in) :: density
    real(kind=dp), intent(out) :: solution(:)

    integer(kind=i4b) :: i, j, nlev
    real(kind=dp) :: out1
    real(kind=dp), allocatable :: rate_matrix(:,:)
    logical :: call_writes

    nlev = size(solution)
    if (size(transition,1).ne.nlev .or. size(transition,2).ne.nlev) then
      stop 'solve_statistical_equilibrium received inconsistent dimensions'
    endif

    allocate(rate_matrix(1:nlev,1:nlev))
    rate_matrix = 0.0D0

    do i=1,nlev
      out1 = 0.0D0
      do j=1,nlev
        out1 = out1 + transition(i,j)
        rate_matrix(i,j) = transition(j,i)
      enddo
      rate_matrix(i,i) = -out1
    enddo

    solution = 0.0D0
    do i=1,nlev
      rate_matrix(nlev,i) = 1.0D-8
    enddo
    solution(nlev) = density*1.0D-8

    call gauss_jordan(rate_matrix,nlev,nlev,solution,call_writes)

    do i=1,nlev
      if (solution(i).lt.0.0D0) solution(i)=0.0D0
    enddo

    deallocate(rate_matrix)
  end subroutine solve_statistical_equilibrium

end module level_population_system_module

subroutine solve_level_population_system(nlev, transition, density, solution)
  use definitions, only : dp
  use healpix_types, only : i4b
  use level_population_system_module, only : solve_statistical_equilibrium
  implicit none

  integer(kind=i4b), intent(in) :: nlev
  real(kind=dp), intent(in) :: density
  real(kind=dp), intent(in) :: transition(1:nlev,1:nlev)
  real(kind=dp), intent(out) :: solution(1:nlev)

  call solve_statistical_equilibrium(transition, density, solution)
end subroutine solve_level_population_system
