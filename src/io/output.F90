module output_module
  use definitions
  use maincode_module
  use global_module

  implicit none

contains

  subroutine write_final_outputs(config_file)
    character(len=*), intent(in) :: config_file
    character(len=32) :: out_file
    character(len=32) :: out_file2

    if (iteration.ge.1) then

      if (iteration.lt.itertot) then
        write(6,*) '3DPDR converged after ',iteration-1,' iterations'
        write(6,'("RESULT status=converged iterations=",I0)') iteration-1
      else
        write(6,*) 'Reached maximum number of iterations without convergence.'
        write(6,*) 'To reach convergence, increase the relative number in [',trim(config_file),']'
        write(6,'("RESULT status=max_iterations iterations=",I0)') iteration
      endif
      write(6,*) 'Writing final outputs'

      out_file = trim(adjustl(output))//".pdr.fin"
      out_file2 = trim(adjustl(out_file))//"]"
      write(6,'(" Writing file [",A)') out_file2
      open(unit=21,file=out_file,status='replace')

      do pp=1,pdr_ptot-2
        p=IDlist_pdr(pp)
#ifdef PSEUDO_1D
        write(21,'(I7,4ES11.3,I5,300ES11.3)') p,pdr(p)%x, pdr(p)%AV(6), previousgastemperature(pp),&
            &dusttemperature(pp),pdr(p)%etype,pdr(p)%rho,pdr(p)%UVfield,pdr(p)%abundance
#else
        write(21,'(I7,5ES11.3,I5,300ES11.3)') p,pdr(p)%x, pdr(p)%y, pdr(p)%z,&
            &previousgastemperature(pp),dusttemperature(pp),pdr(p)%etype,pdr(p)%rho,pdr(p)%UVfield,&
            &pdr(p)%abundance,pdr(p)%AV
#endif
      enddo

      pp=pdr_ptot
      p=IDlist_pdr(pp)
#ifdef PSEUDO_1D
      write(21,'(I7,4ES11.3,I5,300ES11.3)') p,pdr(p)%x, pdr(p)%AV(6), previousgastemperature(pp),&
          &dusttemperature(pp),pdr(p)%etype,pdr(p)%rho,pdr(p)%UVfield,pdr(p)%abundance
#else
      write(21,'(I7,5ES11.3,I5,300ES11.3)') p,pdr(p)%x, pdr(p)%y, pdr(p)%z,&
          &previousgastemperature(pp),dusttemperature(pp),pdr(p)%etype,pdr(p)%rho,pdr(p)%UVfield,&
          &pdr(p)%abundance,pdr(p)%AV
#endif

      pp=pdr_ptot-1
      p=IDlist_pdr(pp)
#ifdef PSEUDO_1D
      write(21,'(I7,4ES11.3,I5,300ES11.3)') p,pdr(p)%x, pdr(p)%AV(6), previousgastemperature(pp),&
          &dusttemperature(pp),pdr(p)%etype,pdr(p)%rho,pdr(p)%UVfield,pdr(p)%abundance
#else
      write(21,'(I7,5ES11.3,I5,300ES11.3)') p,pdr(p)%x, pdr(p)%y, pdr(p)%z,&
          &previousgastemperature(pp),dusttemperature(pp),pdr(p)%etype,pdr(p)%rho,pdr(p)%UVfield,&
          &pdr(p)%abundance,pdr(p)%AV
#endif

      if (ion_ptot.gt.0) then
        out_file = trim(adjustl(output))//".ion.fin"
        out_file2 = trim(adjustl(out_file))//"]"
        write(6,'(" Writing file [",A)') out_file2
        close(21)
        open(unit=21,file=out_file,status='replace')

        do pp=1,ion_ptot
          p=IDlist_ion(pp)
#ifdef PSEUDO_1D
          write(21,'(I7,3ES11.3,I5,300ES11.3)') p,pdr(p)%x, pdr(p)%AV(6),&
              &previousgastemperature(pp),pdr(p)%etype,pdr(p)%rho,pdr(p)%UVfield,pdr(p)%abundance
#else
          write(21,'(I7,4ES11.3,I5,400ES11.3)') p,pdr(p)%x, pdr(p)%y, pdr(p)%z,&
              &previousgastemperature(pp),pdr(p)%etype,pdr(p)%rho,pdr(p)%UVfield,pdr(p)%abundance,pdr(p)%AV(:)
#endif
        enddo
      endif

      if (dark_ptot.gt.0) then
        out_file = trim(adjustl(output))//".mol.fin"
        out_file2 = trim(adjustl(out_file))//"]"
        write(6,'(" Writing file [",A)') out_file2
        close(21)
        open(unit=21,file=out_file,status='replace')

        do pp=1,dark_ptot
          p=IDlist_dark(pp)
#ifdef PSEUDO_1D
          write(21,'(I7,3ES11.3,I5,300ES11.3)') p,pdr(p)%x, pdr(p)%AV(6),&
              &previousgastemperature(pp),pdr(p)%etype,pdr(p)%rho,pdr(p)%UVfield,pdr(p)%abundance
#else
          write(21,'(I7,4ES11.3,I5,400ES11.3)') p,pdr(p)%x, pdr(p)%y, pdr(p)%z,&
              &previousgastemperature(pp),pdr(p)%etype,pdr(p)%rho,pdr(p)%UVfield,pdr(p)%abundance,pdr(p)%AV(:)
#endif
        enddo
      endif

      close(21)

      out_file = trim(adjustl(output))//trim(adjustl(".cool"))//".fin"
      out_file2 = trim(adjustl(out_file))//"]"
      write(6,'(" Writing file [",A)') out_file2
      open(unit=13,file=out_file,status='replace')

      do pp=1,pdr_ptot
        p=IDlist_pdr(pp)
#ifdef PSEUDO_1D
        write(13,'(I7,200ES11.3)') p, pdr(p)%x, pdr(p)%AV(6), CII_cool(pp),CI_cool(pp), &
            &OI_cool(pp),C12O_cool(pp), total_cooling_rate(pp)
#else
        write(13,'(I7,200ES11.3)') p, pdr(p)%x, pdr(p)%y, pdr(p)%z, CII_cool(pp), CI_cool(pp),&
            &OI_cool(pp), C12O_cool(pp), total_cooling_rate(pp), pdr(p)%AV(:)
#endif
      enddo

      out_file = trim(adjustl(output))//trim(adjustl(".heat"))//".fin"
      out_file2 = trim(adjustl(out_file))//"]"
      write(6,'(" Writing file [",A)') out_file2
      open(unit=14,file=out_file,status='replace')

      do pp=1,pdr_ptot
        p=IDlist_pdr(pp)
#ifdef PSEUDO_1D
        write(14,'(I7,200ES11.3)') p, pdr(p)%x, pdr(p)%AV(6), all_heating(pp,:)
#else
        write(14,'(I7,200ES11.3)') p, pdr(p)%x, pdr(p)%y, pdr(p)%z, all_heating(pp,:), pdr(p)%AV(:)
#endif
      enddo
      close(14)

      out_file = trim(adjustl(output))//trim(adjustl(".line"))//".fin"
      out_file2 = trim(adjustl(out_file))//"]"
      write(6,'(" Writing file [",A)') out_file2
      open(unit=16,file=out_file,status='replace')

      do pp=1,pdr_ptot-2
        p=IDlist_pdr(pp)
#ifdef PSEUDO_1D
        write(16,'(I9,200ES11.3)') pp, pdr(p)%x, pdr(p)%AV(6), &
            &pdr(p)%CII_line(2,1),&
            &pdr(p)%CI_line(2,1), pdr(p)%CI_line(3,1), pdr(p)%CI_line(3,2),&
            &pdr(p)%OI_line(2,1), pdr(p)%OI_line(3,1), pdr(p)%OI_line(3,2),&
            &pdr(p)%C12O_line(2,1), pdr(p)%C12O_line(3,2), pdr(p)%C12O_line(4,3),&
            &pdr(p)%C12O_line(5,4), pdr(p)%C12O_line(6,5), pdr(p)%C12O_line(7,6),&
            &pdr(p)%C12O_line(8,7), pdr(p)%C12O_line(9,8), pdr(p)%C12O_line(10,9),&
            &pdr(p)%C12O_line(11,10)
#else
        write(16,'(I9,200ES11.3)') pp, pdr(p)%x, pdr(p)%y, pdr(p)%z, &
            &pdr(p)%CII_line(2,1),&
            &pdr(p)%CI_line(2,1), pdr(p)%CI_line(3,1), pdr(p)%CI_line(3,2),&
            &pdr(p)%OI_line(2,1), pdr(p)%OI_line(3,1), pdr(p)%OI_line(3,2),&
            &pdr(p)%C12O_line(2,1), pdr(p)%C12O_line(3,2), pdr(p)%C12O_line(4,3),&
            &pdr(p)%C12O_line(5,4), pdr(p)%C12O_line(6,5), pdr(p)%C12O_line(7,6),&
            &pdr(p)%C12O_line(8,7), pdr(p)%C12O_line(9,8), pdr(p)%C12O_line(10,9),&
            &pdr(p)%C12O_line(11,10), pdr(p)%AV(:)
#endif
      enddo

      pp=pdr_ptot
      p=IDlist_pdr(pp)
#ifdef PSEUDO_1D
      write(16,'(I9,200ES11.3)') pp, pdr(p)%x, pdr(p)%AV(6), &
          &pdr(p)%CII_line(2,1),&
          &pdr(p)%CI_line(2,1), pdr(p)%CI_line(3,1), pdr(p)%CI_line(3,2),&
          &pdr(p)%OI_line(2,1), pdr(p)%OI_line(3,1), pdr(p)%OI_line(3,2),&
          &pdr(p)%C12O_line(2,1), pdr(p)%C12O_line(3,2), pdr(p)%C12O_line(4,3),&
          &pdr(p)%C12O_line(5,4), pdr(p)%C12O_line(6,5), pdr(p)%C12O_line(7,6),&
          &pdr(p)%C12O_line(8,7), pdr(p)%C12O_line(9,8), pdr(p)%C12O_line(10,9),&
          &pdr(p)%C12O_line(11,10)
#else
      write(16,'(I9,200ES11.3)') pp, pdr(p)%x, pdr(p)%y, pdr(p)%z, &
          &pdr(p)%CII_line(2,1),&
          &pdr(p)%CI_line(2,1), pdr(p)%CI_line(3,1), pdr(p)%CI_line(3,2),&
          &pdr(p)%OI_line(2,1), pdr(p)%OI_line(3,1), pdr(p)%OI_line(3,2),&
          &pdr(p)%C12O_line(2,1), pdr(p)%C12O_line(3,2), pdr(p)%C12O_line(4,3),&
          &pdr(p)%C12O_line(5,4), pdr(p)%C12O_line(6,5), pdr(p)%C12O_line(7,6),&
          &pdr(p)%C12O_line(8,7), pdr(p)%C12O_line(9,8), pdr(p)%C12O_line(10,9),&
          &pdr(p)%C12O_line(11,10), pdr(p)%AV(:)
#endif

      pp=pdr_ptot-1
      p=IDlist_pdr(pp)
#ifdef PSEUDO_1D
      write(16,'(I9,200ES11.3)') pp, pdr(p)%x, pdr(p)%AV(6), &
          &pdr(p)%CII_line(2,1),&
          &pdr(p)%CI_line(2,1), pdr(p)%CI_line(3,1), pdr(p)%CI_line(3,2),&
          &pdr(p)%OI_line(2,1), pdr(p)%OI_line(3,1), pdr(p)%OI_line(3,2),&
          &pdr(p)%C12O_line(2,1), pdr(p)%C12O_line(3,2), pdr(p)%C12O_line(4,3),&
          &pdr(p)%C12O_line(5,4), pdr(p)%C12O_line(6,5), pdr(p)%C12O_line(7,6),&
          &pdr(p)%C12O_line(8,7), pdr(p)%C12O_line(9,8), pdr(p)%C12O_line(10,9),&
          &pdr(p)%C12O_line(11,10)
#else
      write(16,'(I9,200ES11.3)') pp, pdr(p)%x, pdr(p)%y, pdr(p)%z, &
          &pdr(p)%CII_line(2,1),&
          &pdr(p)%CI_line(2,1), pdr(p)%CI_line(3,1), pdr(p)%CI_line(3,2),&
          &pdr(p)%OI_line(2,1), pdr(p)%OI_line(3,1), pdr(p)%OI_line(3,2),&
          &pdr(p)%C12O_line(2,1), pdr(p)%C12O_line(3,2), pdr(p)%C12O_line(4,3),&
          &pdr(p)%C12O_line(5,4), pdr(p)%C12O_line(6,5), pdr(p)%C12O_line(7,6),&
          &pdr(p)%C12O_line(8,7), pdr(p)%C12O_line(9,8), pdr(p)%C12O_line(10,9),&
          &pdr(p)%C12O_line(11,10), pdr(p)%AV(:)
#endif
      close(16)

      out_file = trim(adjustl(output))//trim(adjustl(".spop"))//".fin"
      out_file2 = trim(adjustl(out_file))//"]"
      write(6,'(" Writing file [",A)') out_file2
      open(unit=16,file=out_file,status='replace')

      do pp=1,pdr_ptot
        p=IDlist_pdr(pp)
        write(16,'(I9,200ES11.3)') pp, pdr(p)%Av(6), pdr(p)%CII_pop, pdr(p)%CI_pop, &
            &pdr(p)%OI_pop, pdr(p)%C12O_pop
      enddo
      close(16)

#ifdef PSEUDO_1D
      out_file = trim(adjustl(output))//trim(adjustl(".opdp"))//".fin"
      out_file2 = trim(adjustl(out_file))//"]"
      write(6,'(" Writing file [",A)') out_file2
      open(unit=16,file=out_file,status='replace')

      do pp=1,pdr_ptot
        p=IDlist_pdr(pp)
        write(16,'(I9,200ES11.3)') pp, pdr(p)%x,&
            &pdr(p)%CII_optdepth(2,1,6),&
            &pdr(p)%CI_optdepth(2,1,6), pdr(p)%CI_optdepth(3,1,6), pdr(p)%CI_optdepth(3,2,6),&
            &pdr(p)%OI_optdepth(2,1,6), pdr(p)%OI_optdepth(3,1,6), pdr(p)%OI_optdepth(3,2,6),&
            &pdr(p)%C12O_optdepth(2,1,6), pdr(p)%C12O_optdepth(3,2,6),&
            &pdr(p)%C12O_optdepth(4,3,6), pdr(p)%C12O_optdepth(5,4,6),&
            &pdr(p)%C12O_optdepth(6,5,6), pdr(p)%C12O_optdepth(7,6,6),&
            &pdr(p)%C12O_optdepth(8,7,6), pdr(p)%C12O_optdepth(9,8,6),&
            &pdr(p)%C12O_optdepth(10,9,6), pdr(p)%C12O_optdepth(11,10,6)
      enddo
      close(16)
#endif
    endif
  end subroutine write_final_outputs

  subroutine write_simulation_finished
    write(6,*) ''
    call cpu_time(t4)
    write(6,*) 'Simulation time = ',t4-t2,' seconds.'
    write(6,*) 'Finished !'
  end subroutine write_simulation_finished

end module output_module
