program cmbspec
  use healpix_types
  use params
  use time_mod
  use rec_mod
  use evolution_mod
  implicit none

  ! Initialize time grids
  call initialize_time_mod
  call initialize_rec_mod

  write(*,*) "initializing evolution module"
  call initialize_perturbation_eqns
  call integrate_perturbation_eqns
  call write_to_file_evolution_mod
  write(*,*) "hello world!"
  ! Output to file desired quantities here

end program cmbspec
