module output_module
  use definitions
  use coolants_module, only : coolant_cii, coolant_ci, coolant_oi, coolant_c12o
  use maincode_module
  use global_module

  implicit none

contains

  pure logical function iteration_converged(final_iteration, max_iterations) result(converged)
  integer(kind=i4b), intent(in) :: final_iteration
  integer(kind=i4b), intent(in) :: max_iterations

  converged = final_iteration.lt.max_iterations
end function iteration_converged

pure integer(kind=i4b) function reported_iteration_count(final_iteration, max_iterations) result(iteration_count)
integer(kind=i4b), intent(in) :: final_iteration
integer(kind=i4b), intent(in) :: max_iterations

if (iteration_converged(final_iteration, max_iterations)) then
  iteration_count = max(0_i4b, final_iteration-1_i4b)
else
  iteration_count = final_iteration
end if
end function reported_iteration_count

subroutine write_final_outputs(config_file)
  character(len=*), intent(in) :: config_file
  character(len=32) :: out_file
  character(len=32) :: out_file2
  integer(kind=i4b) :: point_id, point_index
  integer(kind=i4b) :: final_iteration_count

  if (iteration.ge.1) then
    final_iteration_count = reported_iteration_count(iteration, runtime%total_iterations)

    if (iteration_converged(iteration, runtime%total_iterations)) then
      write(6,*) '3DPDR converged after ',final_iteration_count,' iterations'
      write(6,'("RESULT status=converged iterations=",I0)') final_iteration_count
    else
      write(6,*) 'Reached maximum number of iterations without convergence.'
      write(6,*) 'To reach convergence, increase the relative number in [',trim(config_file),']'
      write(6,'("RESULT status=max_iterations iterations=",I0)') final_iteration_count
    end if
    write(6,*) 'Writing final outputs'

    out_file = trim(adjustl(runtime%output_prefix))//".pdr.fin"
    out_file2 = trim(adjustl(out_file))//"]"
    write(6,'(" Writing file [",A)') out_file2
    open(unit=21,file=out_file,status='replace')

    do point_index=1,pdr_ptot-2
      point_id=grid%pdr_ids(point_index)
#ifdef PSEUDO_1D
      write(21,'(I7,4ES11.3,I5,300ES11.3)') point_id,grid%points(point_id)%position(x_axis), grid%points(point_id)%av(6), thermal%previous_gas_temperature(point_index),&
          &thermal%dust_temperature(point_index),grid%points(point_id)%etype,grid%points(point_id)%rho,grid%points(point_id)%uvfield,grid%points(point_id)%abundance
#else
      write(21,'(I7,5ES11.3,I5,300ES11.3)') point_id,grid%points(point_id)%position,&
          &thermal%previous_gas_temperature(point_index),thermal%dust_temperature(point_index),grid%points(point_id)%etype,grid%points(point_id)%rho,grid%points(point_id)%uvfield,&
          &grid%points(point_id)%abundance,grid%points(point_id)%av
#endif
    end do

    point_index=pdr_ptot
    point_id=grid%pdr_ids(point_index)
#ifdef PSEUDO_1D
    write(21,'(I7,4ES11.3,I5,300ES11.3)') point_id,grid%points(point_id)%position(x_axis), grid%points(point_id)%av(6), thermal%previous_gas_temperature(point_index),&
        &thermal%dust_temperature(point_index),grid%points(point_id)%etype,grid%points(point_id)%rho,grid%points(point_id)%uvfield,grid%points(point_id)%abundance
#else
    write(21,'(I7,5ES11.3,I5,300ES11.3)') point_id,grid%points(point_id)%position,&
        &thermal%previous_gas_temperature(point_index),thermal%dust_temperature(point_index),grid%points(point_id)%etype,grid%points(point_id)%rho,grid%points(point_id)%uvfield,&
        &grid%points(point_id)%abundance,grid%points(point_id)%av
#endif

    point_index=pdr_ptot-1
    point_id=grid%pdr_ids(point_index)
#ifdef PSEUDO_1D
    write(21,'(I7,4ES11.3,I5,300ES11.3)') point_id,grid%points(point_id)%position(x_axis), grid%points(point_id)%av(6), thermal%previous_gas_temperature(point_index),&
        &thermal%dust_temperature(point_index),grid%points(point_id)%etype,grid%points(point_id)%rho,grid%points(point_id)%uvfield,grid%points(point_id)%abundance
#else
    write(21,'(I7,5ES11.3,I5,300ES11.3)') point_id,grid%points(point_id)%position,&
        &thermal%previous_gas_temperature(point_index),thermal%dust_temperature(point_index),grid%points(point_id)%etype,grid%points(point_id)%rho,grid%points(point_id)%uvfield,&
        &grid%points(point_id)%abundance,grid%points(point_id)%av
#endif

    if (ion_ptot.gt.0) then
      out_file = trim(adjustl(runtime%output_prefix))//".ion.fin"
      out_file2 = trim(adjustl(out_file))//"]"
      write(6,'(" Writing file [",A)') out_file2
      close(21)
      open(unit=21,file=out_file,status='replace')

      do point_index=1,ion_ptot
        point_id=grid%ion_ids(point_index)
#ifdef PSEUDO_1D
        write(21,'(I7,3ES11.3,I5,300ES11.3)') point_id,grid%points(point_id)%position(x_axis), grid%points(point_id)%av(6),&
            &thermal%previous_gas_temperature(point_index),grid%points(point_id)%etype,grid%points(point_id)%rho,grid%points(point_id)%uvfield,grid%points(point_id)%abundance
#else
        write(21,'(I7,4ES11.3,I5,400ES11.3)') point_id,grid%points(point_id)%position,&
            &thermal%previous_gas_temperature(point_index),grid%points(point_id)%etype,grid%points(point_id)%rho,grid%points(point_id)%uvfield,grid%points(point_id)%abundance,grid%points(point_id)%av(:)
#endif
      end do
    end if

    if (dark_ptot.gt.0) then
      out_file = trim(adjustl(runtime%output_prefix))//".mol.fin"
      out_file2 = trim(adjustl(out_file))//"]"
      write(6,'(" Writing file [",A)') out_file2
      close(21)
      open(unit=21,file=out_file,status='replace')

      do point_index=1,dark_ptot
        point_id=grid%dark_ids(point_index)
#ifdef PSEUDO_1D
        write(21,'(I7,3ES11.3,I5,300ES11.3)') point_id,grid%points(point_id)%position(x_axis), grid%points(point_id)%av(6),&
            &thermal%previous_gas_temperature(point_index),grid%points(point_id)%etype,grid%points(point_id)%rho,grid%points(point_id)%uvfield,grid%points(point_id)%abundance
#else
        write(21,'(I7,4ES11.3,I5,400ES11.3)') point_id,grid%points(point_id)%position,&
            &thermal%previous_gas_temperature(point_index),grid%points(point_id)%etype,grid%points(point_id)%rho,grid%points(point_id)%uvfield,grid%points(point_id)%abundance,grid%points(point_id)%av(:)
#endif
      end do
    end if

    close(21)

    out_file = trim(adjustl(runtime%output_prefix))//trim(adjustl(".cool"))//".fin"
    out_file2 = trim(adjustl(out_file))//"]"
    write(6,'(" Writing file [",A)') out_file2
    open(unit=13,file=out_file,status='replace')

    do point_index=1,pdr_ptot
      point_id=grid%pdr_ids(point_index)
#ifdef PSEUDO_1D
      write(13,'(I7,200ES11.3)') point_id, grid%points(point_id)%position(x_axis), grid%points(point_id)%av(6), coolant_iteration(coolant_cii)%cooling_rate(point_index),coolant_iteration(coolant_ci)%cooling_rate(point_index), &
          &coolant_iteration(coolant_oi)%cooling_rate(point_index),coolant_iteration(coolant_c12o)%cooling_rate(point_index), thermal%total_cooling_rate(point_index)
#else
      write(13,'(I7,200ES11.3)') point_id, grid%points(point_id)%position, coolant_iteration(coolant_cii)%cooling_rate(point_index), coolant_iteration(coolant_ci)%cooling_rate(point_index),&
          &coolant_iteration(coolant_oi)%cooling_rate(point_index), coolant_iteration(coolant_c12o)%cooling_rate(point_index), thermal%total_cooling_rate(point_index), grid%points(point_id)%av(:)
#endif
    end do

    out_file = trim(adjustl(runtime%output_prefix))//trim(adjustl(".heat"))//".fin"
    out_file2 = trim(adjustl(out_file))//"]"
    write(6,'(" Writing file [",A)') out_file2
    open(unit=14,file=out_file,status='replace')

    do point_index=1,pdr_ptot
      point_id=grid%pdr_ids(point_index)
#ifdef PSEUDO_1D
      write(14,'(I7,200ES11.3)') point_id, grid%points(point_id)%position(x_axis), grid%points(point_id)%av(6), all_heating(point_index,:)
#else
      write(14,'(I7,200ES11.3)') point_id, grid%points(point_id)%position, all_heating(point_index,:), grid%points(point_id)%av(:)
#endif
    end do
    close(14)

    out_file = trim(adjustl(runtime%output_prefix))//trim(adjustl(".line"))//".fin"
    out_file2 = trim(adjustl(out_file))//"]"
    write(6,'(" Writing file [",A)') out_file2
    open(unit=16,file=out_file,status='replace')

    do point_index=1,pdr_ptot-2
      point_id=grid%pdr_ids(point_index)
      call write_line_cooling_row(16, point_index, point_id)
    end do

    point_index=pdr_ptot
    point_id=grid%pdr_ids(point_index)
    call write_line_cooling_row(16, point_index, point_id)

    point_index=pdr_ptot-1
    point_id=grid%pdr_ids(point_index)
    call write_line_cooling_row(16, point_index, point_id)
    close(16)

    out_file = trim(adjustl(runtime%output_prefix))//trim(adjustl(".spop"))//".fin"
    out_file2 = trim(adjustl(out_file))//"]"
    write(6,'(" Writing file [",A)') out_file2
    open(unit=16,file=out_file,status='replace')

    do point_index=1,pdr_ptot
      point_id=grid%pdr_ids(point_index)
      write(16,'(I9,200ES11.3)') point_index, grid%points(point_id)%av(6), &
          &grid%points(point_id)%coolant_state(coolant_cii)%population, &
          &grid%points(point_id)%coolant_state(coolant_ci)%population, &
          &grid%points(point_id)%coolant_state(coolant_oi)%population, &
          &grid%points(point_id)%coolant_state(coolant_c12o)%population
    end do
    close(16)

#ifdef PSEUDO_1D
    out_file = trim(adjustl(runtime%output_prefix))//trim(adjustl(".opdp"))//".fin"
    out_file2 = trim(adjustl(out_file))//"]"
    write(6,'(" Writing file [",A)') out_file2
    open(unit=16,file=out_file,status='replace')

    do point_index=1,pdr_ptot
      point_id=grid%pdr_ids(point_index)
      call write_optical_depth_row(16, point_index, point_id)
    end do
    close(16)
#endif
  end if
end subroutine write_final_outputs

subroutine write_line_cooling_row(output_unit, point_index, point_id)
  integer(kind=i4b), intent(in) :: output_unit
  integer(kind=i4b), intent(in) :: point_index
  integer(kind=i4b), intent(in) :: point_id

  associate(cii_emission => grid%points(point_id)%coolant_state(coolant_cii)%line, &
        &ci_emission => grid%points(point_id)%coolant_state(coolant_ci)%line, &
        &oi_emission => grid%points(point_id)%coolant_state(coolant_oi)%line, &
        &co_emission => grid%points(point_id)%coolant_state(coolant_c12o)%line)
#ifdef PSEUDO_1D
    write(output_unit,'(I9,200ES11.3)') point_index, grid%points(point_id)%position(x_axis), grid%points(point_id)%av(6), &
        &cii_emission(2,1), &
        &ci_emission(2,1), ci_emission(3,1), ci_emission(3,2), &
        &oi_emission(2,1), oi_emission(3,1), oi_emission(3,2), &
        &co_emission(2,1), co_emission(3,2), co_emission(4,3), &
        &co_emission(5,4), co_emission(6,5), co_emission(7,6), &
        &co_emission(8,7), co_emission(9,8), co_emission(10,9), &
        &co_emission(11,10)
#else
    write(output_unit,'(I9,200ES11.3)') point_index, grid%points(point_id)%position, &
        &cii_emission(2,1), &
        &ci_emission(2,1), ci_emission(3,1), ci_emission(3,2), &
        &oi_emission(2,1), oi_emission(3,1), oi_emission(3,2), &
        &co_emission(2,1), co_emission(3,2), co_emission(4,3), &
        &co_emission(5,4), co_emission(6,5), co_emission(7,6), &
        &co_emission(8,7), co_emission(9,8), co_emission(10,9), &
        &co_emission(11,10), grid%points(point_id)%av(:)
#endif
  end associate
end subroutine write_line_cooling_row

subroutine write_optical_depth_row(output_unit, point_index, point_id)
  integer(kind=i4b), intent(in) :: output_unit
  integer(kind=i4b), intent(in) :: point_index
  integer(kind=i4b), intent(in) :: point_id

  associate(cii_depth => grid%points(point_id)%coolant_state(coolant_cii)%optical_depth, &
        &ci_depth => grid%points(point_id)%coolant_state(coolant_ci)%optical_depth, &
        &oi_depth => grid%points(point_id)%coolant_state(coolant_oi)%optical_depth, &
        &co_depth => grid%points(point_id)%coolant_state(coolant_c12o)%optical_depth)
    write(output_unit,'(I9,200ES11.3)') point_index, grid%points(point_id)%position(x_axis), &
        &cii_depth(2,1,6), &
        &ci_depth(2,1,6), ci_depth(3,1,6), ci_depth(3,2,6), &
        &oi_depth(2,1,6), oi_depth(3,1,6), oi_depth(3,2,6), &
        &co_depth(2,1,6), co_depth(3,2,6), &
        &co_depth(4,3,6), co_depth(5,4,6), &
        &co_depth(6,5,6), co_depth(7,6,6), &
        &co_depth(8,7,6), co_depth(9,8,6), &
        &co_depth(10,9,6), co_depth(11,10,6)
  end associate
end subroutine write_optical_depth_row

subroutine write_simulation_finished(simulation_start_time)
  real, intent(in) :: simulation_start_time
  real :: simulation_end_time

  write(6,*) ''
  call cpu_time(simulation_end_time)
  write(6,*) 'Simulation time = ',simulation_end_time-simulation_start_time,' seconds.'
  write(6,*) 'Finished !'
end subroutine write_simulation_finished

end module output_module
