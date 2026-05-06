module coolant_io_module
  use healpix_types, only : i4b
  use coolants_module, only : coolant_count
  use coolant_input_module, only : read_lamda_coolant_file
  use maincode_module

contains

  subroutine load_coolant_data
    integer(kind=i4b) :: coolant_id

    do coolant_id = 1, coolant_count
      call read_lamda_coolant_file(coolant(coolant_id))
    end do

    write(6,*) ''
  end subroutine load_coolant_data

end module coolant_io_module
