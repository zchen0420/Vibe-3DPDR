program test_ray_path
  use definitions, only : dp
  use healpix_types, only : i4b
  use simulation_grid_module, only : pdr_node, simulation_grid
  use ray_path_module, only : adaptive_step_length, initialize_ray_origin, projected_point_id, &
      & ray_step_length, store_projected_point, update_minimum_adaptive_step
  implicit none

  type(pdr_node) :: point
  type(simulation_grid) :: grid
  real(kind=dp) :: adaptive_minimum
  real(kind=dp) :: origin(1:3)
  logical :: kill_ray

  origin = (/0.0D0, 0.0D0, 0.0D0/)
  allocate(point%epray(0:1))
  allocate(point%epoint(1:3,0:1,0:2))
  allocate(point%projected(0:1,0:2))

  call initialize_ray_origin(point, 7_i4b, origin)
  point%epray(0) = 2_i4b
  point%epoint(1:3,0,1) = (/3.0D0, 4.0D0, 0.0D0/)
  point%epoint(1:3,0,2) = (/3.0D0, 4.0D0, 12.0D0/)
  point%projected(0,1) = 8_i4b
  point%projected(0,2) = 9_i4b

  if (projected_point_id(point, 0_i4b, 0_i4b).ne.7_i4b) stop 'origin projection not initialized'
  if (projected_point_id(point, 0_i4b, 2_i4b).ne.9_i4b) stop 'projected point lookup failed'

  if (abs(ray_step_length(point, 0_i4b, 1_i4b)-5.0D0).gt.1.0D-12) stop 'bad first ray step'
  if (abs(ray_step_length(point, 0_i4b, 2_i4b)-12.0D0).gt.1.0D-12) stop 'bad second ray step'
  if (abs(adaptive_step_length(point, 0_i4b, 2_i4b)-8.0D0).gt.1.0D-12) stop 'bad adaptive step'

  adaptive_minimum = 100.0D0
  call update_minimum_adaptive_step(point, adaptive_minimum)
  if (abs(adaptive_minimum-5.0D0).gt.1.0D-12) stop 'bad adaptive minimum'

  allocate(grid%points(1:2))
  allocate(grid%points(1)%epray(0:0))
  allocate(grid%points(1)%epoint(1:3,0:0,0:1))
  allocate(grid%points(1)%projected(0:0,0:1))
  grid%points(2)%etype = 2_i4b
  call initialize_ray_origin(grid%points(1), 1_i4b, origin)

  kill_ray = .false.
  call store_projected_point(grid, 1_i4b, origin, (/1.0D0, 2.0D0, 3.0D0/), 0_i4b, 2_i4b, 1_i4b, kill_ray)

  if (grid%points(1)%epray(0).ne.1_i4b) stop 'store did not increment epray'
  if (projected_point_id(grid%points(1), 0_i4b, 1_i4b).ne.2_i4b) stop 'store projected id failed'
  if (abs(grid%points(1)%epoint(3,0,1)-3.0D0).gt.1.0D-12) stop 'store projected position failed'
  if (.not.kill_ray) stop 'ionized target should kill ray'

  write(6,*) 'test_ray_path: ok'
end program test_ray_path
