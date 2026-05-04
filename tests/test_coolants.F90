program test_coolants
  use coolants_module
  implicit none

  call assert_label(COOLANT_CII, 'CII ')
  call assert_label(COOLANT_CI, 'CI  ')
  call assert_label(COOLANT_OI, 'OI  ')
  call assert_label(COOLANT_C12O, 'CO  ')
  call assert_int(coolant_default_nlevels(COOLANT_CII), 5, 'CII nlevels')
  call assert_int(coolant_default_nlevels(COOLANT_C12O), 41, 'CO nlevels')
  call assert_int(coolant_default_ntemps(COOLANT_CI), 29, 'CI ntemps')
  call assert_int(coolant_default_ntemps(COOLANT_OI), 27, 'OI ntemps')

  write(6,*) 'test_coolants: ok'

contains

  subroutine assert_label(coolant_id, expected)
    integer, intent(in) :: coolant_id
    character(len=*), intent(in) :: expected

    if (coolant_label(coolant_id).ne.expected) then
      write(6,*) 'test_coolants failed:', coolant_id, coolant_label(coolant_id), expected
      stop 1
    endif
  end subroutine assert_label

  subroutine assert_int(actual, expected, label)
    integer, intent(in) :: actual
    integer, intent(in) :: expected
    character(len=*), intent(in) :: label

    if (actual.ne.expected) then
      write(6,*) 'test_coolants failed:', trim(label), actual, expected
      stop 1
    endif
  end subroutine assert_int

end program test_coolants
