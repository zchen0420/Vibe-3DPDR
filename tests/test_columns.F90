program test_columns
  use definitions, only : dp
  use columns_module, only : column_increment
  implicit none

  real(kind=dp) :: value

  value = column_increment(10.0D0, 30.0D0, 0.5D0, 0.25D0, 2.0D0)
  if (abs(value-12.5D0).gt.1.0D-12) then
    stop 'Unexpected column increment'
  endif

  write(6,*) 'test_columns: ok'
end program test_columns
