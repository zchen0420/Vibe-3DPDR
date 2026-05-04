module spatial_index_module
  use definitions
  use healpix_types
  use maincode_module

contains

  subroutine prepare_spatial_index_inputs
    call swap_last_two_pdr_points
    call build_pdr_point_index
  end subroutine prepare_spatial_index_inputs

  subroutine swap_last_two_pdr_points
    real(kind=dp), allocatable :: x_rev(:), y_rev(:), z_rev(:), n_rev(:)

    if (pdr_ptot.lt.2) return

    allocate(x_rev(1:pdr_ptot))
    allocate(y_rev(1:pdr_ptot))
    allocate(z_rev(1:pdr_ptot))
    allocate(n_rev(1:pdr_ptot))

    do pp=1,pdr_ptot-2
      p=IDlist_pdr(pp)
      x_rev(pp)=pdr(p)%x
      y_rev(pp)=pdr(p)%y
      z_rev(pp)=pdr(p)%z
      n_rev(pp)=pdr(p)%rho
    enddo

    x_rev(pdr_ptot-1)=pdr(IDlist_pdr(pdr_ptot))%x
    y_rev(pdr_ptot-1)=pdr(IDlist_pdr(pdr_ptot))%y
    z_rev(pdr_ptot-1)=pdr(IDlist_pdr(pdr_ptot))%z
    n_rev(pdr_ptot-1)=pdr(IDlist_pdr(pdr_ptot))%rho

    x_rev(pdr_ptot)=pdr(IDlist_pdr(pdr_ptot-1))%x
    y_rev(pdr_ptot)=pdr(IDlist_pdr(pdr_ptot-1))%y
    z_rev(pdr_ptot)=pdr(IDlist_pdr(pdr_ptot-1))%z
    n_rev(pdr_ptot)=pdr(IDlist_pdr(pdr_ptot-1))%rho

    do pp=1,pdr_ptot
      p=IDlist_pdr(pp)
      pdr(p)%x=x_rev(pp)
      pdr(p)%y=y_rev(pp)
      pdr(p)%z=z_rev(pp)
      pdr(p)%rho=n_rev(pp)
    enddo

    deallocate(x_rev)
    deallocate(y_rev)
    deallocate(z_rev)
    deallocate(n_rev)
  end subroutine swap_last_two_pdr_points

  subroutine build_pdr_point_index
    allocate(rra(0:pdr_ptot))
    allocate(rrb(1:pdr_ptot))
    allocate(pdrpoint(1:3,0:pdr_ptot)) !0 is for the ONE molecular element

    do pp=1,pdr_ptot
      p=IDlist_pdr(pp)
      rra(pp) = sqrt(pdr(p)%x**2+pdr(p)%y**2+pdr(p)%z**2)
      rrb(pp) = p
      pdrpoint(1,pp) = pdr(p)%x
      pdrpoint(2,pp) = pdr(p)%y
      pdrpoint(3,pp) = pdr(p)%z
    enddo

    if (dark_ptot.gt.0) then
      pdrpoint(1,0) = pdr(IDlist_dark(1))%x
      pdrpoint(2,0) = pdr(IDlist_dark(1))%y
      pdrpoint(3,0) = pdr(IDlist_dark(1))%z
    endif
  end subroutine build_pdr_point_index

end module spatial_index_module
