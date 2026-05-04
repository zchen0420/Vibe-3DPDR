program test_radiation
  use definitions, only : dp
  use healpix_types, only : i4b
  use radiation_module, only : attenuated_uv, incident_ray_index
  implicit none

  real(kind=dp) :: field_vector(1:3)
  real(kind=dp) :: ray_vectors(1:3,0:2)
  integer(kind=i4b) :: selected_ray

  field_vector = (/1.0D0, 0.0D0, 0.0D0/)
  ray_vectors(:,0) = (/1.0D0, 0.0D0, 0.0D0/)
  ray_vectors(:,1) = (/-1.0D0, 0.0D0, 0.0D0/)
  ray_vectors(:,2) = (/0.0D0, 1.0D0, 0.0D0/)

  selected_ray = incident_ray_index(field_vector, ray_vectors, 3_i4b)
  if (selected_ray.ne.1_i4b) stop 'Unexpected incident ray index'

  if (abs(attenuated_uv(2.0D0, 0.0D0, 3.0D0)-2.0D0).gt.1.0D-12) then
    stop 'Unexpected unattenuated UV value'
  endif

  if (attenuated_uv(1.0D0, 1000.0D0, 1000.0D0).ne.0.0D0) then
    stop 'Expected UV floor to zero'
  endif

  write(6,*) 'test_radiation: ok'
end program test_radiation
