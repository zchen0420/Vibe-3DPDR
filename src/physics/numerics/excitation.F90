module excitation_module

contains

  ! Calculate the partition function for the given coolant.
  subroutine calculate_partition_function(partition_function,nlev,energies,weights,temperature)

    use healpix_types
    implicit none

    integer(kind=i4b), intent(in) :: nlev
    real(kind=dp), intent(in)     :: energies(1:nlev),weights(1:nlev)
    real(kind=dp), intent(in)     :: temperature
    real(kind=dp), intent(out)    :: partition_function

    integer(kind=i4b) :: ilevel

    partition_function=0.0d0
    do ilevel=1,nlev
      partition_function=partition_function + weights(ilevel)*exp(-energies(ilevel)/kb/temperature)
    end do

    return
  end subroutine calculate_partition_function

  ! Calculate LTE level populations for the given coolant.
  subroutine calculate_lte_populations(nlev,level_pop,energies,weights,partition_function,density,temperature)

    use healpix_types
    implicit none

    integer(kind=i4b), intent(in) :: nlev
    real(kind=dp), intent(in)     :: energies(1:nlev),weights(1:nlev)
    real(kind=dp), intent(in)     :: partition_function
    real(kind=dp), intent(in)     :: density,temperature
    real(kind=dp), intent(out)    :: level_pop(1:nlev)

    integer(kind=i4b) :: ilevel
    real(kind=dp) :: total_pop

    total_pop=0.0d0
    do ilevel=1,nlev
      level_pop(ilevel)=density*weights(ilevel)*exp(-energies(ilevel)/kb/temperature)/partition_function
      total_pop=total_pop + level_pop(ilevel)
    end do

    ! Check that the sum of the level populations adds up to the total density.
    if(abs(total_pop-density)/density.gt.1.0d-3) then
      write(6,*)'ERROR! Sum of LTE level populations differs from the total density by ', &
          & int(1.0d2*abs(total_pop-density)/density),'%'
      stop
    end if

    return
  end subroutine calculate_lte_populations

end module excitation_module
