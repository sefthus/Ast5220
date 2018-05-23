module cl_mod
  use healpix_types
  use evolution_mod
  use sphbess_mod
  use spline_1D_mod
  implicit none

  real(dp),     allocatable, dimension(:,:)     :: S, S2
  integer(i4b), allocatable, dimension(:)       :: ls
  real(dp),     allocatable, dimension(:)       :: ls_dp, ls_hires, cls_hires
  real(dp),     allocatable, dimension(:)       :: x_hires, k_hires,integrand1
  real(dp),     allocatable, dimension(:,:)     :: Theta_l
  real(dp),     allocatable, dimension(:,:)     :: integrand2
  !real(dp),     parameter, private :: n_s=0.96

contains

  ! Driver routine for (finally!) computing the CMB power spectrum
  subroutine compute_cls
    implicit none

    integer(i4b) :: i, j, l, l_num, x_num, n_spline, j_loc
    real(dp)     :: dx, S_func, j_func, z, eta, eta0, x0, x_min, x_max, d, e

    real(dp),     allocatable,     dimension(:,:)     :: j_l, j_l2
    real(dp),     allocatable,     dimension(:)       :: x_arg, int_arg, cls, cls2
    real(dp),     allocatable,     dimension(:)       :: k, x
    real(dp),     allocatable,     dimension(:,:,:,:) :: S_coeff

    real(dp),     allocatable, dimension(:)       :: z_spline, j_l_spline, j_l_spline2
    real(dp),     allocatable, dimension(:)       :: x_lores, besseltest
    real(dp),     allocatable, dimension(:)       :: eta_arr

    real(dp)           :: t1, t2, integral, integral1, integral2, h1, h2, int1

    logical(lgt)       :: exist
    character(len=128) :: filename1,filename
    real(dp), allocatable, dimension(:) :: y, y2
    real(dp) :: start_time, end_time

    ! Set up which l's to compute
    l_num = 44
    allocate(ls(l_num))
    ls = (/ 2, 3, 4, 6, 8, 10, 12, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100, &
         & 120, 140, 160, 180, 200, 225, 250, 275, 300, 350, 400, 450, 500, 550, &
         & 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200 /)
    !ls = (/6,100, 200,500,1000,1200/)

    ! Task: Get source function from evolution_mod
    allocate(S(n_k_hires,n_x_hires))
    allocate(x_hires(n_x_hires))
    allocate(k_hires(n_k_hires))

    filename1 = 'source.bin'
    inquire(file=filename1, exist=exist)
    if (exist) then
       write(*,*) 'reading Source function from file'
       open(0, form='unformatted', file=filename1)
       read(0) x_hires, k_hires, S
       close(0)
    else
        write(*,*) "initializing evolution module"
        call initialize_perturbation_eqns
        call integrate_perturbation_eqns

        call get_hires_source_function(k_hires, x_hires, S)
        write(*,*) 'made the call to evol_mod'
        open(0, form='unformatted', file=filename1)
        write(0) x_hires, k_hires, S
        close(0)
    end if

    ! Task: Initialize spherical Bessel functions for each l; use 5400 sampled points between 
    !       z = 0 and 3500. Each function must be properly splined
    ! Hint: It may be useful for speed to store the splined objects on disk in an unformatted
    !       Fortran (= binary) file, so that these only has to be computed once. Then, if your
    !       cache file exists, read from that; if not, generate the j_l's on the fly.
    n_spline = 5400
    allocate(z_spline(n_spline))    ! Note: z is *not* redshift, but simply the dummy argument of j_l(z)
    allocate(j_l(n_spline, l_num))
    allocate(j_l2(n_spline, l_num))

    write(*,*) 'making z_spline'
    do i=1,n_spline
      z_spline(i) = 0.d0 + (i-1)*(3400.d0-0.d0)/(n_spline-1.d0)
    end do
    
    ! checking for binary file
    filename = 'j_l.bin'
    inquire(file=filename, exist=exist)
    if (exist) then
       write(*,*) 'reading Bessel functions from file'
       open(1, form='unformatted', file=filename)
       read(1) j_l, j_l2
       close(1)
    else
      write(*,*) "compute spherical Bessel functions"
      do l=1, l_num
        j_l(1,l) = 0.d0
        do i=2, n_spline
          call sphbes(ls(l),z_spline(i), j_l(i,l))
        end do
      end do

      ! spline bessel functions across z for all l
      write(*,*) "spline bessel functions"
      do l=1,l_num
        call spline(z_spline, j_l(:,l), 1.0d30, 1.0d30, j_l2(:,l))
      end do

      ! Write to file
      open(1, form='unformatted', file=filename)
      write(1) j_l, j_l2
      close(1)
    end if

    j_loc = locate_dp(k_hires,340.d0*H_0/c)

    allocate(besseltest(n_x_hires))
    open (unit=3 ,file="besseltest.dat",action="write",status="replace")
    do i =1,n_x_hires 
      besseltest(i) = splint(z_spline,j_l(:,17),j_l2(:,17),k_hires(j_loc)*(get_eta(0.d0)-get_eta(x_hires(i))))
      write (3 ,*) besseltest(i)
    end do
    close (3)
    !stop

    allocate(Theta_l(l_num,n_k_hires))
    allocate(int_arg(n_x_hires))
    allocate(integrand1(n_x_hires))
    allocate(integrand2(l_num, n_k_hires))
    allocate(cls(l_num))
    allocate(cls2(l_num))

    !allocate(x_lores(n_x_hires/10)) ! creating low-res x grid for fast Theta_l-integration

    ! Overall task: Compute the C_l's for each given l
    ! Precompute eta0-eta(i)
    allocate(eta_arr(n_x_hires))
    eta0 = get_eta(0.d0)      
    do i=1, n_x_hires
      eta_arr(i) = eta0-get_eta(x_hires(i))
    end do

    ! For integration
    h1 = (x_hires(n_x_hires) - x_hires(1))/n_x_hires
    h2 = (k_hires(n_k_hires) - k_hires(1))/n_k_hires

    do l = 1, l_num
	  write(*,*) 'l=', l
      ! Task: Compute the transfer function, Theta_l(k)
      integral2 = 0.d0               ! Reset C_l integration
      write(*,*) 'integration for theta_l'
      !Start timer
      call cpu_time(start_time)

      do j=1, n_k_hires

        integral1= 0.d0              ! Reset Theta integration
        do i=1, n_x_hires
          integrand1(i) = S(j,i)*splint(z_spline,j_l(:,l),j_l2(:,l), k_hires(j)*eta_arr(i))
          !int1 = S(j,i)*splint(z_spline,j_l(:,l),j_l2(:,l), k_hires(j)*k_eta(i))
          integral1 = integral1 + integrand1(i) !int1
        end do
        Theta_l(l,j) = h1*integral1  !-0.5d0*(integrand1(1)+integrand1(n_x_hires)))

        if(l==17 .and. j==j_loc) then   ! Save l=100 and k=340 compare with Callin
          write(*,*)'writing integrand to file for l=17,k=',j_loc
          open (unit=2 ,file="Sj_l.dat",action="write",status="replace")
             do i=1,n_x_hires
               write (2 ,*) integrand1(i)
             end do 
             close (2)
            !stop
         end if

       ! Task: Integrate P(k) * (Theta_l^2 / k) over k to find un-normalized C_l's
        integrand2(l,j) = (c*k_hires(j)/H_0)**(n_s-1.d0)*Theta_l(l,j)**2/k_hires(j)
        integral2 = integral2 + integrand2(l,j)
      end do 
      write(*,*) 'integration for cl'
      integral2 = h2*integral2! - 0.5d0*(integrand2(l,1)+integrand2(l,n_k_hires)))

      ! Task: Store C_l in an array. Optionally output to file
      cls(l) = integral2*ls(l)*(ls(l)+1.d0)/(2.d0*pi)
      !Print time used
      call cpu_time(end_time)

      print'("Time used = ",f7.2," seconds.")',end_time-start_time
      write(*,*) 'a11'

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
  
  !call write_to_file_cl_mod
  end subroutine compute_cls

  subroutine write_to_file_cl_mod
     use healpix_types
    implicit none

    integer(i4b) :: i,j
    integer(i4b), dimension(6) :: l_val,l
    integer(i4b), dimension(44) :: ll

    !l_val(1:6)=(/6, 100, 200, 500, 900, 1200 /) ! k
    !!l=(/4,17,22,30,38,44/)
    ! Finds index of chosen l values
    l_val(1:6)=(/6, 100, 200, 500, 1000, 1200/)
    do j=1, 6 
        forall (i=1:44) ll(i) = abs(ls(i)-l_val(j))
        l(j) = minloc(ll,1)
    end do
    write(*,*) "writing to file; cl_mod"   
 
!---------- write to file ---
    write(*,*) "opening files "
    open (unit=1, file = 'l_val.dat', status='replace')
    open (unit=2, file = 'x_k_hires.dat', status='replace')
    open (unit=3, file = 'Theta_l.dat', status='replace')
    open (unit=4, file = 'C_l_integrand.dat', status='replace')
    open (unit=5, file = 'C_l.dat', status='replace')

    write(*,*) "writing stuff"
    write(*,*) '    writing chosen k values'
    do i=1, 6 ! write the k values used
      write(1, *) l_val(i)      
    end do

    write(*,*) '    writing x, k, source func., Theta_intl100, C_l100_int,Theta_l, C_l_int'
    do i=1, n_x_hires
      write (2,'(*(2X, ES14.6E3))') x_hires(i), k_hires(i)
      write (3,'(*(2X, ES14.6E3))')&
        Theta_l(l(1), i),Theta_l(l(2), i),Theta_l(l(3), i),Theta_l(l(4),i ),Theta_l(l(5), i),Theta_l(l(6), i)
      write (4,'(*(2X, ES14.6E3))')&
        integrand2(l(1), i)/(c*k_hires(i)/H_0)**(n_s-1.d0),&
        integrand2(l(2), i)/(c*k_hires(i)/H_0)**(n_s-1.d0),&
        integrand2(l(3), i)/(c*k_hires(i)/H_0)**(n_s-1.d0),&
        integrand2(l(4), i)/(c*k_hires(i)/H_0)**(n_s-1.d0),&
        integrand2(l(5), i)/(c*k_hires(i)/H_0)**(n_s-1.d0),&
        integrand2(l(6), i)/(c*k_hires(i)/H_0)**(n_s-1.d0)
    end do

    write(*,*) '    writing l_hires and c_l_hires'
    do i=1,1200
      write (5,'(*(2X, ES14.6E3))') ls_hires(i), cls_hires(i)
    end do
    
    write(*,*) 'closing files'
    do i=1, 5
      close(i)
    end do

  end subroutine write_to_file_cl_mod
end module cl_mod
