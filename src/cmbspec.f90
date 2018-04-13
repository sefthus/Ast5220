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
  call initialize_evolution_mod
  write(*,*) "hello wordl!"
  ! Output to file desired quantities here

end program cmbspec
