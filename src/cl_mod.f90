module cl_mod
  use healpix_types
  use evolution_mod
  use sphbess_mod
  use spline_1D_mod
  implicit none

  real(dp),     pointer,     dimension(:,:)     :: S, S2
  integer(i4b), allocatable, dimension(:)       :: ls
  real(dp), allocatable, dimension(:)       :: ls_hires, cls_hires
  real(dp),     allocatable, dimension(:)       :: x_hires, k_hires
  real(dp),     allocatable, dimension(:,:)     :: Theta_l
contains

  ! Driver routine for (finally!) computing the CMB power spectrum
  subroutine compute_cls
    implicit none

    integer(i4b) :: i, j, l, l_num, x_num, n_spline
    real(dp)     :: dx, S_func, j_func, z, eta, eta0, x0, x_min, x_max, d, e

    real(dp),     allocatable, dimension(:)       :: integrand, integrand2
    real(dp),     pointer,     dimension(:,:)     :: j_l, j_l2
    real(dp),     pointer,     dimension(:)       :: x_arg, int_arg, cls, cls2, ls_dp
    real(dp),     pointer,     dimension(:)       :: k, x
    real(dp),     pointer,     dimension(:,:,:,:) :: S_coeff

    !real(dp),     allocatable, dimension(:,:)     :: integrand2
    real(dp),     allocatable, dimension(:)       :: z_spline, j_l_spline, j_l_spline2
    real(dp),     allocatable, dimension(:)       :: x_lores

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
    write(*,*) "hei 1"
    ! Task: Get source function from evolution_mod
    allocate(S(n_x_hires,n_k_hires))
    allocate(x_hires(n_x_hires))
    allocate(k_hires(n_k_hires))
    call get_hires_source_function(k_hires, x_hires, S)
    write(*,*) 'made the call to evol_mod'
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
            call sphbes(ls(l),z_spline(i),j_l(i,l))
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

    allocate(Theta_l(l_num,n_k_hires))
    allocate(int_arg(n_x_hires))
    allocate(integrand(n_x_hires))
    allocate(integrand2(n_k_hires))
    allocate(cls(l_num))
    allocate(cls2(l_num))
    allocate(x_lores(n_x_hires/10)) ! creating low-res x grid for fast Theta_l-integration

    ! Overall task: Compute the C_l's for each given l
    do l = 1, l_num
       ! Task: Compute the transfer function, Theta_l(k)
      ! start trapezoidal x
      write(*,*) 'start trapezoidal for theta_l'
      do j=1, n_k_hires       
	!do j = 1, n_x_hires
          !do i =1, n_x_hires/10
          !  if (i<=300) then
          !      m = 1 + (1-a)*(i_rec-1)/(300-1)
          !  else
          !      m = i_rec + (i-301)*(n_x_hires-i_rec)/199
          !  end if

          !  x_lores(i) = x_hires(m)
        do i=1, n_x_hires
          integrand(i) = S(i,j)*splint(z_spline,j_l(:,l),j_l2(:,l), k_hires(j)*(get_eta(0.d0)-get_eta(x_hires(i))))
        end do

        call trapz(x_hires,integrand,Theta_l(l,j))
      end do
        ! end trapezoidal x

       ! Task: Integrate P(k) * (Theta_l^2 / k) over k to find un-normalized C_l's
       ! start trapezoidal k
      write(*,*) 'start trapezoidal for cl'
      do j=1,n_k_hires
        integrand2(j) = (c*k_hires(j)/H_0)**(n_s-1.d0)*Theta_l(l,j)**2/k_hires(j)
      end do 
      call trapz(k_hires,integrand2, integral)
       ! Task: Store C_l in an array. Optionally output to file
       cls(l) = integral*ls(l)*(ls(l)+1.d0)/(2.d0*pi)
        ! end trapezoidal method k
    end do


    ! Task: Spline C_l's found above, and output smooth C_l curve for each integer l
    allocate(ls_dp(l_num))
    allocate(ls_hires(int(maxval(ls))))
    allocate(cls_hires(int(maxval(ls))))

    do l=1, l_num   ! spline requires double precision
        ls_dp(l) = ls(l)
    end do

    write(*,*) 'splining cls'
    call spline(ls_dp, cls, 1.d30, 1.d30, cls2)

    ! new unit stepsize l-grid
    write(*,*) 'making new highresolution l-grid'
    do l=1, int(maxval(ls))
      ls_hires(l) = l
    end do

    ! find Cls for all ls_hires
    write(*,*) 'saving splined cls'
    do l=1, int(maxval(ls))
      cls_hires(l) = splint(ls_dp, cls, cls2, ls_hires(l))
    end do

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
    integral = integral*h
  end subroutine trapz


  subroutine write_to_file_cl_mod
     use healpix_types
    implicit none

    integer(i4b) :: i, l
    integer(i4b), dimension(6) :: j
    j(1:6)=(/50, 250, 500, 2000, 3000, 5000 /)
    !j(1:6)=(/1, 2, 3, 4, 5, 10 /)
    write(*,*) "writing to file; cl_mod"    

!---------- write to file ---
    write(*,*) "opening files "
    open (unit=0, file = 'x_hires.dat', status='replace')
    open (unit=1, file = 'source_func.dat', status='replace')
    open (unit=2, file = 'Theta_l.dat', status='replace')
    open (unit=3, file = 'ls.dat', status='replace')
    open (unit=4, file = 'C_l.dat', status='replace')
    !open (unit=4, file = 'C_l_integrand.dat', status='replace')
    
    write(*,*) "writing stuff"
    do i=1, n_x_hires
      write (0,*) x_hires(i)
      write (1,'(*(2X, ES14.6E3))') S(i,j(1)),S(i,j(2)),S(i,j(3)),S(i,j(4)),S(i,j(5)),S(i,j(6))
      !write(11,'(*(2X, ES14.6E3))') Theta_l(l,j(1)),Theta_l(l,j(2)),Theta_l(l,j(3)),Theta_l(l,j(4)),Theta_l(l,j(5)),Theta_l(l,j(6))
    end do

    do l=1, 44
       write (1,'(*(2X, ES14.6E3))') Theta_l(l,j(1)),Theta_l(l,j(2)),Theta_l(l,j(3)),Theta_l(l,j(4)),Theta_l(l,j(5)),Theta_l(l,j(6))
        write (2,*) ls(l)
    end do

    do i=1,1200
      write (3,'(*(2X, ES14.6E3))') ls_hires(i),cls_hires(i)
    end do
    
    write(*,*) "closing files "
    do i=0, 4
      close(i)
    end do

  end subroutine write_to_file_cl_mod
end module cl_mod
