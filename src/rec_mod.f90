module rec_mod
  use healpix_types
  use params
  use time_mod
  use ode_solver
  use spline_1D_mod
  implicit none

  integer(i4b),                        private :: n                 ! Number of grid points
  real(dp), allocatable, dimension(:), private :: x_rec             ! Grid
  real(dp), allocatable, dimension(:), private :: a_rec             ! Grid

  real(dp), allocatable, dimension(:), private :: tau, tau2, tau22  ! Splined tau and second derivatives
  real(dp), allocatable, dimension(:), private :: n_e, n_e2         ! Splined (log of) electron density, n_e
  real(dp), allocatable, dimension(:), private :: g, g2, g22        ! Splined visibility function

contains

  subroutine initialize_rec_mod
    implicit none
    
    integer(i4b) :: i, j, k
    real(dp)     :: saha_limit, y, T_b, n_b, dydx, xmin, xmax, dx, f, n_e0, X_e0, xstart, xstop
    logical(lgt) :: use_saha
    real(dp), allocatable, dimension(:) :: X_e ! Fractional electron density, n_e / n_H
    real(dp)     :: dx_rec, step, stepmin
    real(dp), dimension(1) :: y

    saha_limit = 0.99d0       ! Switch from Saha to Peebles when X_e < 0.99
    xstart     = log(1.d-10)  ! Start grids at a = 10^-10
    xstop      = 0.d0         ! Stop  grids at a = 1
    n          = 1000         ! Number of grid points between xstart and xstopo

    eps        = 1.d-8        ! spline error limit    

    allocate(x_rec(n))
    allocate(X_e(n))
    allocate(tau(n))
    allocate(tau2(n))
    allocate(tau22(n))
    allocate(n_e(n))
    allocate(n_e2(n))
    allocate(g(n))
    allocate(g2(n))
    allocate(g22(n))

    ! Task: Fill in x (rec) grid
    x_rec(1) = xstart
    dx_rec   = (xstop-xstart)/(n-1)

    do i = 1, n
       x_rec(i+1) = xstart + i*dx_rec

    a_rec = exp(x_rec)


    ! Task: Compute X_e and n_e at all grid times
    step       = abs((x_rec(1) - x_rec(2))/ 100.d0)     !integration step length
    stepmin    = abs((x_rec(1) - x_rec(2))/ 10000.d0)

    use_saha = .true.
    do i = 1, n-1
       T_b = T_0/a_rec(i)
       n_b = Omega_b*rho_c/(m_H*a_rec(i)**3)
       if (use_saha) then
          ! Use the Saha equation
          K = 1.d0/(n_b * hbar**3) * ((m_e*k_b*T_b)/(2*pi))**(3.d0/2.d0) * exp(-epsilon_0/(k_b * T_b)s)
          X_e(i)= (-K + sqrt(K**2 + 4.d0*K))/2.d0

          if (X_e(i) < saha_limit) use_saha = .false.
       else
          ! Use the Peebles equation
          !y(1)=X_e(i)
          X_e(i+1) = X_e(i)
          call odeint(X_e(i+1:i+1), x_rec(i), x_rec(i+1), eps, step, stepmin, derivs, bsstep, output)


       end if
    end do
    
    n_e = X_e*n_b       ! electron density

    ! Task: Compute splined (log of) electron density function
    n_e = log(n_e)
    call spline(x_rec, n_e, 1d30, 1d30, n_e2)

    ! Task: Compute optical depth at all grid points


    ! Task: Compute splined (log of) optical depth
    ! Task: Compute splined second derivative of (log of) optical depth


    ! Task: Compute splined visibility function
    ! Task: Compute splined second derivative of visibility function


  end subroutine initialize_rec_mod

!---------------------- Peebles equation ----
  subroutine dXe_dx(x, X_e, dydx)
    ! we define dy/dx
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: X_e
    real(dp), dimension(:), intent(out) :: dydx
    real(dp) :: beta, beta2, alpha2, n_b, n_1s, lambda_21s, lambda_alpha
    real(dp) :: C_r, T_b, H, phi2 
    
    H    = get_H(x)
    T_b  = T_0/exp(x)
    n_b  = Omega_b*rho_c/(m_H*exp(x)**3)
    n1s  = (1 - X_e)* n_b

    phi2 = 0.448*log(epsilon_0/(k_b * T_b))
    alpha2 = 64*pi/sqrt(27*pi) * alpha**2/m_e**2 * sqrt(epsilon_0/T_b)*phi2

    beta = alpha2*(me*kb*Tb)**(1.5d0) * exp(-epsilon_0/(k_b * T_b))
    beta2 = beta * exp(3*epsilon_0/(4 * k_b*T_b))

    lambda_alpha = H * (3*epsilon_0)**3/((8*pi)**2 * n1s)
    lambda_21s = 8.227

    C_r = (lambda_21s + lambda_alpha)/(lambda_21s + lambda_alpha + beta2)
    
    dydx = C_r/H * (beta * (1-X_e) - n_H * alpha2 * X_e**2)
  end subroutine dXe_dx

  ! Task: Complete routine for computing n_e at arbitrary x, using precomputed information
  ! Hint: Remember to exponentiate...
  function get_n_e(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_n_e

  end function get_n_e

  ! Task: Complete routine for computing tau at arbitrary x, using precomputed information
  function get_tau(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_tau

  end function get_tau

  ! Task: Complete routine for computing the derivative of tau at arbitrary x, using precomputed information
  function get_dtau(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_dtau

  end function get_dtau

  ! Task: Complete routine for computing the second derivative of tau at arbitrary x, 
  ! using precomputed information
  function get_ddtau(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_ddtau

  end function get_ddtau

  ! Task: Complete routine for computing the visibility function, g, at arbitray x
  function get_g(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_g

  end function get_g

  ! Task: Complete routine for computing the derivative of the visibility function, g, at arbitray x
  function get_dg(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_dg

  end function get_dg

  ! Task: Complete routine for computing the second derivative of the visibility function, g, at arbitray x
  function get_ddg(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_ddg

  end function get_ddg



end module rec_mod
