module cl_mod
  use healpix_types
  use evolution_mod
  use sphbess_mod
  implicit none

contains

  ! Driver routine for (finally!) computing the CMB power spectrum
  subroutine compute_cls
    implicit none

    integer(i4b) :: i, j, l, l_num, x_num, n_spline
    real(dp)     :: dx, S_func, j_func, z, eta, eta0, x0, x_min, x_max, d, e
    integer(i4b), allocatable, dimension(:)       :: ls
    real(dp),     allocatable, dimension(:)       :: integrand
    real(dp),     pointer,     dimension(:,:)     :: j_l, j_l2
    real(dp),     pointer,     dimension(:)       :: x_arg, int_arg, cls, cls2, ls_dp
    real(dp),     pointer,     dimension(:)       :: k, x
    real(dp),     pointer,     dimension(:,:,:,:) :: S_coeff
    real(dp),     pointer,     dimension(:,:)     :: S, S2
    real(dp),     allocatable, dimension(:,:)     :: Theta
    real(dp),     allocatable, dimension(:,:)     :: Theta_l,integrand2
    real(dp),     allocatable, dimension(:)       :: z_spline, j_l_spline, j_l_spline2
    real(dp),     allocatable, dimension(:)       :: x_hires, k_hires, x_lores

    real(dp)           :: t1, t2, integral
    logical(lgt)       :: exist
    character(len=128) :: filename
    real(dp), allocatable, dimension(:) :: y, y2

    ! Set up which l's to compute
    l_num = 44
    allocate(ls(l_num))
    ls = (/ 2, 3, 4, 6, 8, 10, 12, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100, &
         & 120, 140, 160, 180, 200, 225, 250, 275, 300, 350, 400, 450, 500, 550, &
         & 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200 /)

    ! Task: Get source function from evolution_mod
    allocate(S(n_x_highres,n_k_highres))
    allocate(x_hires(n_x_highres))
    allocate(k_hires(n_k_highres))
    call get_hires_source_function(k_hires, x_hires, S)

    ! Task: Initialize spherical Bessel functions for each l; use 5400 sampled points between 
    !       z = 0 and 3500. Each function must be properly splined
    ! Hint: It may be useful for speed to store the splined objects on disk in an unformatted
    !       Fortran (= binary) file, so that these only has to be computed once. Then, if your
    !       cache file exists, read from that; if not, generate the j_l's on the fly.
    n_spline = 5400
    allocate(z_spline(n_spline))    ! Note: z is *not* redshift, but simply the dummy argument of j_l(z)
    allocate(j_l(n_spline, l_num))
    allocate(j_l2(n_spline, l_num))

    do i=1,n_spline
      z_spline(i) = (i-1)*3500.d0/(n_spline-1.d0)
    end do

    ! checking for binary file
    filename = 'j_l.bin'
    inquire(file=filename, exist=exist)
    if (exist) then
       open(10, form='unformatted', file=filename)
       read(10) j_l, j_l2
       close(10)
    else
      write(*,*) "compute spherical Bessel functions"
      do i=1, n_spline
        do l=1, l_num
          if (z_spline(i) > 2.d0) then
              call sphbes(ls(l),z_pline(i),j_l(i,l))
          end if
        end do
      end do

        ! spline bessel functions across z for all l
      write(*,*) "spline bessel functions"
      do l=1,l_num
        call spline(z_spline, j_l(:,1), 1.0d30, 1.0d30, j_l2(:,1))
      end do
      ! Write to file
      open(10, form='unformatted', file=filename)
      write(10) j_l, j_l2
      close(10)
    end if

    allocate(Theta(l_num,n_k_highres))
    allocate(int_arg(n_hires))
    allocate(integrand(n_x_highres))
    allocate(integrand2(l_num,n_k_highres))
    allocate(cls(l_num))
    allocate(cls2(l_num))
    allocate(x_lores(n_hires/10)) ! creating low-res x grid for fast Theta-integration

    ! Overall task: Compute the C_l's for each given l
    do l = 1, l_num
       ! Task: Compute the transfer function, Theta_l(k)
       do j = 1, n_hires
          do i =1, n_hires/10
            if (i<=300) then
                m = 1 + (1-a)*(i_rec-1)/(300-1)
            else
                m = i_rec + (i-301)*(n_hires-i_rec)/199
            end if

            x_lores(i) = x_hires(m)
            integrand(i) = S(m,j)*splint(z_spline,j_l(:,l),j_l2(:,l), k_hires(j)*(get_eta(0.d0)-get_eta(x_hires(m))))
            end do
          call trapz(x_lores,integrand(1:500),Theta(l,j))

       ! Task: Integrate P(k) * (Theta_l^2 / k) over k to find un-normalized C_l's


       ! Task: Store C_l in an array. Optionally output to file

    end do


    ! Task: Spline C_l's found above, and output smooth C_l curve for each integer l


  end subroutine compute_cls
  
  subroutine trapz(x,y,integral) ! Trapezoidal method
    implicit none
    real(dp), dimension(:), intent(in)  :: x,y
    real(dp),               intent(out) :: integral
    integer(i4b)                        :: n, i
    real(dp)                            :: h

    if (size(x) .ne. size(y)) then
       write(*,*) 'x and y does not have same shape'
       return
    end if

    integral = 0.d0
    n = size(x)
    h = (x(n) - x(1))/n
    do i=1, n-1
       integral = integral + h*(y(i)+y(i+1))/2.d0
    end do

  end subroutine trapz

end module cl_mod
