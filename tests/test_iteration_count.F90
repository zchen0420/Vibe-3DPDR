program test_iteration_count
  use healpix_types, only : i4b
  use output_module, only : iteration_converged, reported_iteration_count
  implicit none

  call assert_true(iteration_converged(228_i4b, 500_i4b), 'Expected convergence before max iterations')
  call assert_int(reported_iteration_count(228_i4b, 500_i4b), 227_i4b, &
      &'Expected converged result to report completed previous iteration')

  call assert_false(iteration_converged(500_i4b, 500_i4b), 'Expected max-iteration status at limit')
  call assert_int(reported_iteration_count(500_i4b, 500_i4b), 500_i4b, &
      &'Expected max-iteration result to report the final attempted iteration')

  call assert_int(reported_iteration_count(1_i4b, 500_i4b), 0_i4b, &
      &'Expected first-iteration convergence to report zero completed iterations')

  write(6,*) 'test_iteration_count: ok'

contains

  subroutine assert_true(actual, label)
    logical, intent(in) :: actual
    character(len=*), intent(in) :: label

    if (.not.actual) then
      write(6,*) 'test_iteration_count failed: ', trim(label)
      stop 1
    endif
  end subroutine assert_true

  subroutine assert_false(actual, label)
    logical, intent(in) :: actual
    character(len=*), intent(in) :: label

    if (actual) then
      write(6,*) 'test_iteration_count failed: ', trim(label)
      stop 1
    endif
  end subroutine assert_false

  subroutine assert_int(actual, expected, label)
    integer(kind=i4b), intent(in) :: actual
    integer(kind=i4b), intent(in) :: expected
    character(len=*), intent(in) :: label

    if (actual.ne.expected) then
      write(6,*) 'test_iteration_count failed: ', trim(label), actual, expected
      stop 1
    endif
  end subroutine assert_int

end program test_iteration_count
