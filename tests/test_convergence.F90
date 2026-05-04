program test_convergence
  use definitions, only : dp
  use healpix_types, only : i4b
  use convergence_module, only : population_has_converged
  implicit none

  real(kind=dp) :: solution(1:3)
  real(kind=dp) :: current(1:3)
  real(kind=dp) :: relch(1:3)

  solution = (/1.0D0, 2.0D0, 3.0D0/)
  current = (/1.0001D0, 1.9999D0, 3.0001D0/)
  if (.not.population_has_converged(solution, current, 1.0D0, relch, 3_i4b)) then
    stop 'Expected converged level populations'
  endif

  current = (/1.5D0, 2.0D0, 3.0D0/)
  if (population_has_converged(solution, current, 1.0D0, relch, 3_i4b)) then
    stop 'Expected non-converged level populations'
  endif

  write(6,*) 'test_convergence: ok'
end program test_convergence
