program harmonic_oscillator
  use healpix_types
  use ode_solver
  use spline_1D_mod
  implicit none

  integer(i4b) :: i, n, n_hr
  real(dp)     :: dt, x_0, v_0, eps, dt_hr
  real(dp), dimension(2) :: y
  real(dp), allocatable, dimension(:) :: t, v, x, x2, v2, t_hr, x_hr, v_hr

  n = 50                     ! number of points to output functions
  n_hr = 1000                ! number of point in high-resolution grid
  dt = 10.d0 * pi / (n - 1)  ! time step in ode_solver d 10 exp double precision
  dt_hr = 10.d0 * pi / (n_hr - 1)  ! time step of hr-grid
  x_0 = 0.d0                 ! initial position
  v_0 = 1.d0                 ! initial velocity
  eps = 1.d-8                ! error tolerance in ode_solver function
  
  allocate(t(n))
  allocate(x(n))
  allocate(v(n))
  do i=1,n
     t(i) = (i-1) * dt
  end do
  ! solve ode-system:
  x(1) = x_0
  v(1) = v_0
  y(1) = x(1)
  y(2) = v(1)

  do i=1,n-1
     call odeint(y, t(i), t(i+1), eps, dt / 10.d0, dt / 10000.d0, derivs, bsstep, output)
     x(i+1) = y(1)
     v(i+1) = y(2)
  end do
  ! output to file
  open(54, file='lowres_data.dat')
  do i=1,n
     write(54, '(3(E17.8))') t(i), x(i), v(i)
  end do
  close(54)
  write(*,*) "Done writing lowres data to file"

  ! spline lowres data
  allocate(x2(n))
  allocate(v2(n))
  ! spline takes arguments (x, y, dy/dx|_1, dy/dx|_n, d^2y/dx^2),
  ! where the last argument is the output of the subroutine.
  ! Setting the two first derivatives to 1e30 corresponds to 
  ! choosing the "natural spline". 
  call spline(t, x, 1d30, 1d30, x2)
  call spline(t, v, 1d30, 1d30, v2)

  ! interpolate to high resolution grid: 
  allocate(t_hr(n_hr))
  allocate(x_hr(n_hr))
  allocate(v_hr(n_hr))
  do i=1,n_hr
     t_hr(i) = (i-1) * dt_hr
     x_hr(i) = splint(t, x, x2, t_hr(i)) ! spline interpolation
     v_hr(i) = splint(t, v, v2, t_hr(i)) ! spline interpolation
  end do

  ! to save memory, we deallocate arrays we no longer need
  deallocate(t)
  deallocate(x)
  deallocate(v)
  deallocate(x2)
  deallocate(v2)

  ! output to file
  open(54, file='highres_data.dat')
  do i=1,n_hr
     write(54, '(3(E17.8))') t_hr(i), x_hr(i), v_hr(i)
  end do
  close(54)
  write(*,*) "Done writing highres data to file"

  deallocate(t_hr)
  deallocate(x_hr)
  deallocate(v_hr)

contains

  subroutine derivs(x, y, dydx)
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: y
    real(dp), dimension(:), intent(out) :: dydx
    real(dp)  :: k_over_m
    k_over_m = 1.d0
    dydx(1) = y(2) ! = v(1)
    dydx(2) = -k_over_m * y(1) 
  end subroutine derivs

  subroutine output(x, y)
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: y
  end subroutine output
  
end program harmonic_oscillator
