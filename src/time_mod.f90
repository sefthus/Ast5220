module time_mod
  use healpix_types
  use params
  use ode_solver
  use spline_1D_mod
  implicit none

  integer(i4b)                           :: n_t                ! Number of x-values
  real(dp),    allocatable, dimension(:) :: x_t                ! Grid of relevant x-values
  real(dp),    allocatable, dimension(:) :: a_t                ! Grid of relevant a-values

  integer(i4b)                           :: n_eta              ! Number of eta grid poins
  real(dp),    allocatable, dimension(:) :: x_eta              ! Grid points for eta
  real(dp),    allocatable, dimension(:) :: eta, eta2          ! Eta and eta'' at each grid point

contains

  subroutine initialize_time_mod
    implicit none

    integer(i4b) :: i, n, n1, n2, n0
    real(dp)     :: z_start_rec, z_end_rec, z_0, x_start_rec, x_end_rec, x_0, dx, x_eta1, x_eta2, a_init
    real(dp)     :: x_init, x_int1, x_int2, x_eta_int, eps, step, stepmin, eta_init, rho_c
    real(dP)     :: H_scale, Omega_mx, Omega_bx, Omega_rx, Omega_lambdax, z
    real(dp), dimension(1) :: y 

    ! Define two epochs, 1) during and 2) after recombination.
    n0          = 300                       ! Number of grid points before recombination
    n1          = 200                       ! Number of grid points during recombination
    n2          = 300                       ! Number of grid points after recombination
    n_t         = n1 + n2 + n0                 ! Total number of grid points
    z_start_rec = 1630.4d0                  ! Redshift of start of recombination
    z_end_rec   = 614.2d0                   ! Redshift of end of recombination
    z_0         = 0.d0                      ! Redshift today
    x_start_rec = -log(1.d0 + z_start_rec)  ! x of start of recombination
    x_end_rec   = -log(1.d0 + z_end_rec)    ! x of end of recombination
    x_0         = 0.d0                      ! x today
    
    n_eta       = 1000                      ! Number of eta grid points (for spline)
    !a_init      = 1.d-10                    ! Start value of a for eta evaluation
    a_init      = 1.d-8                    ! Start value of a for eta evaluation

    x_eta1      = log(a_init)               ! Start value of x for eta evaluation
    x_init      = log(a_init)
    x_eta2      = 0.d0                      ! End value of x for eta evaluation

    eps         = 1.d-8                        ! spline error limit 
    eta_init    = c/get_H_p(x_eta1)!*a_init/(H0*sqrt(Omega_r)) ! eta initial value at a=0

    ! Task: Fill in x and a grids
    allocate(x_t(0:n_t))
    allocate(a_t(0:n_t))

    x_int1 = (x_end_rec - x_start_rec) /n1
    x_int2 = (x_0 - x_end_rec) /n2

    ! x grid before recombination !! Milestone 3
    do i=0,n0
       x_t(i) = x_init + (i)*(x_start_rec-x_init)/(n0-1)
    end do

    ! x grid during recombination
    !x_t(1) = x_start_rec
    do i=1,n1
       x_t(n0+i) = x_start_rec + i*x_int1
    end do

    ! x grid after recombination
    do i=1,n2-1
       x_t(n0+n1+i) = x_end_rec + i*x_int2
    end do

    ! a grid values
    a_t = exp(x_t) ! x = ln a


    ! Task: 1) Compute the conformal time at each eta time step
    !       2) Spline the resulting function, using the provided "spline" routine in spline_1D_mod.f90
    allocate(x_eta(n_eta))
    allocate(eta(n_eta))
    allocate(eta2(n_eta))
    
    ! x_eta grid
    x_eta_int = (x_eta2 - x_eta1)/(n_eta - 1)
    x_eta(1) = x_eta1
    do i=1,n_eta-1
       x_eta(i+1) = x_eta(i) + x_eta_int
    end do

    ! integrating to find eta

    step = abs((x_eta(1) - x_eta(2))/ 100.d0)
    stepmin = abs((x_eta(1) - x_eta(2))/ 10000.d0)
    
    y(1) = eta_init
    eta(1) = eta_init
    do i=1,n_eta-1       
       call odeint(y, x_eta(i), x_eta(i+1), eps, step, stepmin, derivs, bsstep, output)
       eta(i+1) = y(1)
    end do

    ! calling spline on eta
    call spline(x_eta, eta, 1d30, 1d30, eta2)    



    ! write stuff to file - x_eta. eta. eta splint, H, Omegas
    !open (unit=1, file = 'xt_eta_t.dat', status='replace')
    !open (unit=2, file = 'omega_mbrl.dat', status='replace')
    !open (unit=3, file = 'xeta_eta.dat', status='replace')
    !open (unit=4, file = 'xeta_z_H.dat', status='replace')

    !do i=1,n_t 
    !   write (1,'(2(E17.8))') x_t(i), get_eta(x_t(i))
    !end do
    
    !do i=1, n_eta

      ! calculate and write Omegas
    !  H_scale       = H_0/get_H(x_eta(i))
    !  Omega_mx      = Omega_m      * H_scale**2 * exp(x_eta(i))**(-3)
    !  Omega_bx      = Omega_b      * H_scale**2 * exp(x_eta(i))**(-3)
    !  Omega_rx      = Omega_r      * H_scale**2 * exp(x_eta(i))**(-4)
    !  Omega_lambdax = Omega_lambda * H_scale**2

    !  z = exp(-x_eta(i))-1
    !  write (2,'(4(E17.8))') Omega_mx, Omega_bx, Omega_rx, Omega_lambdax
    !  write (3,'(2(E17.8))') x_eta(i), eta(i)
    !  write (4,'(3(E17.8))') x_eta(i), z, get_H(x_eta(i))

    !end do

    !do i=1,4 ! close files
    !  close(i)
    !end do

  end subroutine initialize_time_mod


  subroutine derivs(x, eta, detadx)
    ! we define d eta/d x
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: eta
    real(dp), dimension(:), intent(out) :: detadx
    
    detadx = c/get_H_p(x) 
  end subroutine derivs

  subroutine output(x, y)
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: y
  end subroutine output

  ! Task: Write a function that computes H at given x
  function get_H(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_H
    real(dp)             :: a
    a = exp(x)
    get_H = H_0*sqrt((Omega_b+Omega_m)*a**(-3) + (Omega_r + Omega_nu)*a**(-4) + Omega_lambda)
  end function get_H

  ! Task: Write a function that computes H' = a*H  at given x
  function get_H_p(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_H_p
    real(dp)             :: a
    a = exp(x)
    get_H_p = a*get_H(x)

  end function get_H_p

  ! Task: Write a function that computes dH'/dx at given x
  function get_dH_p(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_dH_p
    real(dp)             :: a
    a = exp(x)    
    get_dH_p = -H_0**2*(0.5d0*(Omega_b-Omega_m)*a**(-2) + (Omega_r + Omega_nu)*a**(-3) - Omega_lambda*a)/get_H(x)
  end function get_dH_p

  ! Task: Write a function that computes eta(x), using the previously precomputed splined function
  function get_eta(x_in)
    implicit none

    real(dp), intent(in) :: x_in
    real(dp)             :: get_eta
    get_eta = splint(x_eta, eta, eta2, x_in)

  end function get_eta

end module time_mod
