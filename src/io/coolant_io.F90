module coolant_io_module
  use healpix_types, only : i4b
  use coolants_module, only : COOLANT_COUNT
  use maincode_module

contains

  subroutine load_coolant_data
    integer(kind=i4b) :: coolant_id

    do coolant_id = 1, COOLANT_COUNT
      call readinput(coolant(coolant_id)%input_file,coolant(coolant_id)%nlevels,coolant(coolant_id)%ntemperatures,&
          &     coolant(coolant_id)%energies,coolant(coolant_id)%weights,&
          &     coolant(coolant_id)%a_coeffs,coolant(coolant_id)%b_coeffs,coolant(coolant_id)%frequencies,&
          &     coolant(coolant_id)%temperatures,coolant(coolant_id)%h,coolant(coolant_id)%hp,&
          &     coolant(coolant_id)%el,coolant(coolant_id)%he,coolant(coolant_id)%h2,&
          &     coolant(coolant_id)%ph2,coolant(coolant_id)%oh2)
    enddo

    write(6,*) ''
  end subroutine load_coolant_data

end module coolant_io_module
