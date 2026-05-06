program test_level_population_system
  use definitions, only : dp
  use level_population_system_module, only : solve_statistical_equilibrium
  implicit none

  real(kind=dp), parameter :: density = 30.0D0
  real(kind=dp) :: transition(2,2)
  real(kind=dp) :: solution(2)

  transition = 0.0D0
  transition(1,2) = 2.0D0
  transition(2,1) = 1.0D0

  call solve_statistical_equilibrium(transition, density, solution)

  call assert_close(solution(1), 10.0D0, 1.0D-10, 'solution(1)')
  call assert_close(solution(2), 20.0D0, 1.0D-10, 'solution(2)')
  call assert_close(sum(solution), density, 1.0D-10, 'population sum')

  write(6,*) 'test_level_population_system: ok'

contains

  subroutine assert_close(actual, expected, tolerance, label)
    real(kind=dp), intent(in) :: actual
    real(kind=dp), intent(in) :: expected
    real(kind=dp), intent(in) :: tolerance
    character(len=*), intent(in) :: label

    if (abs(actual - expected) > tolerance) then
      write(6,*) 'test_level_population_system failed: ', trim(label), actual, expected
      stop 1
    endif
  end subroutine assert_close

end program test_level_population_system
