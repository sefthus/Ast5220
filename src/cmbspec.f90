program cmbspec
  use healpix_types
  use params
  use time_mod
  use rec_mod
  use evolution_mod
  use cl_mod
  implicit none

  ! Initialize time grids
  write(*,*) "initializing time module"
  call initialize_time_mod

  write(*,*) "initializing rec module"
  call initialize_rec_mod

  write(*,*) "initializing evolution module"
  call initialize_perturbation_eqns
  call integrate_perturbation_eqns

  write(*,*) " initialize cl module"
  call compute_cls

  write(*,*) "hello world!"
  ! Output to file desired quantities here
 
  !call write_to_file_time_mod
  !call write_to_file_rec_mod
  !call write_to_file_evolution_mod
  call write_to_file_cl_mod
end program cmbspec
