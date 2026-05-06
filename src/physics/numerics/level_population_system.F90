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
    end if

    allocate(rate_matrix(1:nlev,1:nlev))
    rate_matrix = 0.0d0

    do i=1,nlev
      out1 = 0.0d0
      do j=1,nlev
        out1 = out1 + transition(i,j)
        rate_matrix(i,j) = transition(j,i)
      end do
      rate_matrix(i,i) = -out1
    end do

    solution = 0.0d0
    do i=1,nlev
      rate_matrix(nlev,i) = 1.0d-8
    end do
    solution(nlev) = density*1.0d-8

    call gauss_jordan(rate_matrix,nlev,nlev,solution,call_writes)

    do i=1,nlev
      if (solution(i).lt.0.0d0) solution(i)=0.0d0
    end do

    deallocate(rate_matrix)
  end subroutine solve_statistical_equilibrium

end module level_population_system_module
