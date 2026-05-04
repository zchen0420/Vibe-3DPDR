module columns_module
  use definitions, only : dp
  use healpix_types, only : i4b
  implicit none

contains

  real(kind=dp) function column_increment(left_density, right_density, left_abundance, right_abundance, step_length_pc)
    real(kind=dp), intent(in) :: left_density
    real(kind=dp), intent(in) :: right_density
    real(kind=dp), intent(in) :: left_abundance
    real(kind=dp), intent(in) :: right_abundance
    real(kind=dp), intent(in) :: step_length_pc

    column_increment = step_length_pc*(left_density*left_abundance + right_density*right_abundance)/2.0D0
  end function column_increment

  subroutine calc_columndens
    use maincode_module

    integer(kind=i4b) :: point_id
    integer(kind=i4b) :: pdr_index
    integer(kind=i4b) :: ray_index
    integer(kind=i4b) :: eval_index
    integer(kind=i4b) :: species_index
    real(kind=dp) :: ray_step

    if (dark_ptot.gt.0) then
      point_id=IDlist_dark(1)
      if (referee.eq.0) allocate(column(0)%columndens_point(0:nrays-1,1:nspec))
      column(0)%columndens_point = 0.0D0
      call calculate_point_columns(point_id, 0)
    endif

    do pdr_index=1,pdr_ptot
      point_id=IDlist_pdr(pdr_index)
      pdr(point_id)%projected(:,0)=point_id
#ifdef THERMALBALANCE
      if (converged(pdr_index)) cycle
#endif
      if (referee.eq.0) allocate(column(pdr_index)%columndens_point(0:nrays-1,1:nspec))
      column(pdr_index)%columndens_point = 0.0D0
      call calculate_point_columns(point_id, pdr_index)
    enddo

  contains

    subroutine calculate_point_columns(point_id, column_index)
      integer(kind=i4b), intent(in) :: point_id
      integer(kind=i4b), intent(in) :: column_index
      integer(kind=i4b) :: left_point_id
      integer(kind=i4b) :: right_point_id

      do ray_index=0,nrays-1
        if (pdr(point_id)%epray(ray_index).eq.0) cycle
        do eval_index=1,pdr(point_id)%epray(ray_index)
          ray_step = sqrt((pdr(point_id)%epoint(1,ray_index,eval_index-1)-&
              &pdr(point_id)%epoint(1,ray_index,eval_index))**2+&
              &(pdr(point_id)%epoint(2,ray_index,eval_index-1)-&
              &pdr(point_id)%epoint(2,ray_index,eval_index))**2 + &
              &(pdr(point_id)%epoint(3,ray_index,eval_index-1)-&
              &pdr(point_id)%epoint(3,ray_index,eval_index))**2)
          left_point_id = int(pdr(point_id)%projected(ray_index,eval_index-1))
          right_point_id = int(pdr(point_id)%projected(ray_index,eval_index))

          do species_index=1,nspec
            column(column_index)%columndens_point(ray_index,species_index) = &
                &column(column_index)%columndens_point(ray_index,species_index) + &
                &column_increment(pdr(left_point_id)%rho, pdr(right_point_id)%rho, &
                &pdr(left_point_id)%abundance(species_index), &
                &pdr(right_point_id)%abundance(species_index), ray_step*pc)
          enddo
        enddo
      enddo
    end subroutine calculate_point_columns

  end subroutine calc_columndens

end module columns_module
