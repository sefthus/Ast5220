program cmbspec
  use healpix_types
  use params
  use time_mod
  use rec_mod
  use evolution_mod
  use cl_mod
  implicit none
  !Time keeping variables
  real(dp) :: start_time,end_time

  !Start timer
  call cpu_time(start_time)

  ! Initialize time grids
  write(*,*) "initializing time module"
  call initialize_time_mod

  write(*,*) "initializing rec module"
  call initialize_rec_mod

  write(*,*) "initializing evolution module"
  !call initialize_perturbation_eqns
  !call integrate_perturbation_eqns

  write(*,*) " initialize cl module"
  call compute_cls

  write(*,*) "hello world!"
  ! Output to file desired quantities here
 
  !call write_to_file_time_mod
  !call write_to_file_rec_mod
  !call write_to_file_evolution_mod
  !call write_to_file_cl_mod

  !Print time used
  call cpu_time(end_time)
  print'("Time used = ",f7.2," seconds.")',end_time-start_time

end program cmbspec
