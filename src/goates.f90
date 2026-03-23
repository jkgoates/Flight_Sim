module goates_m

    implicit none

    ! CONSTANTS
    real,parameter :: PI = 3.141592653589793
    real,parameter :: PI2 = PI*0.5
    real,parameter :: PI4 = PI*0.25
    real,parameter :: g_ssl = 9.80665 ! Standard Gravity at sea level 
    real,parameter :: r_ez = 6356766.0 ! US standard atmosphere earth radius in meters
    real,parameter :: r_e = 6366707.01949371 ! US standard atmosphere earth radius in meters
    real,parameter :: R_gas = 287.0528
    real,parameter :: gamma = 1.4
    real,parameter :: P_ssl = 101325.0 ! Standard Sea Level Pressure in N/m^2
    real, parameter :: one_sixth = 1./6.

    ! CONVERSION FACTORS
    real,parameter :: Pa_to_lbf_ft2 = 1.0/47.880258
    real,parameter :: kg_m3_to_slug_ft3 = 1.0/515.379   

    ! ATMOSPHERE TABLES
    real, parameter :: Z_i(8)= [0., 11000., 20000., 32000., 47000., 52000., 61000., 79000.]
    real, parameter :: Z_i1(8) = [11000., 20000., 32000., 47000., 52000., 61000., 79000., 90000.]
    real, parameter :: T_i(8)= [288.150, 216.650, 216.650, 228.650, 270.650, 270.650, 252.650, 180.650]
    real, parameter :: Tp_i(8) = [-6.5, 0.0, 1.0, 2.8, 0.0, -2.0, -4.0, 0.0]
    real, parameter :: P_i(8)= [1.01325e5, 2.2632031822221168e4, 5.4748735282708267e3, 8.6801476908672271e2, &
                1.1090558898922531e2, 5.9000524278924367e1, 1.8209924905017658e1, 1.0377004548920223]

    ! OTHER STUFF
    logical :: verbose, save_states
    
contains
    
    !!! COORDINATE TRANSFORM !!!

    function quat_mult(q1, q2) result(q3)

        implicit none
        
        real, intent(in) :: q1(4), q2(4)
        real :: q3(4)

        q3(1) = q1(1)*q2(1) - q1(2)*q2(2) - q1(3)*q2(3) - q1(4)*q2(4)
        q3(2) = q1(1)*q2(2) + q1(2)*q2(1) + q1(3)*q2(4) - q1(4)*q2(3)
        q3(3) = q1(1)*q2(3) - q1(2)*q2(4) + q1(3)*q2(1) + q1(4)*q2(2)
        q3(4) = q1(1)*q2(4) + q1(2)*q2(3) - q1(3)*q2(2) + q1(4)*q2(1)

    end function quat_mult


    function quat_base_to_dependent(vec1, q1) result(vec2)

        implicit none

        real, intent(in) :: q1(4), vec1(3)
        real :: vec2(3)

        real :: t(4)

        t(1) = -vec1(1)*q1(2) - vec1(2)*q1(3) - vec1(3)*q1(4)
        t(2) =  vec1(1)*q1(1) + vec1(2)*q1(4) - vec1(3)*q1(3)
        t(3) = -vec1(1)*q1(4) + vec1(2)*q1(1) + vec1(3)*q1(2)
        t(4) =  vec1(1)*q1(3) - vec1(2)*q1(2) + vec1(3)*q1(1)

        vec2(1) = q1(1)*t(2) - q1(2)*t(1) - q1(3)*t(4) + q1(4)*t(3)
        vec2(2) = q1(1)*t(3) + q1(2)*t(4) - q1(3)*t(1) - q1(4)*t(2)
        vec2(3) = q1(1)*t(4) - q1(2)*t(3) + q1(3)*t(2) - q1(4)*t(1)

    end function quat_base_to_dependent


    function quat_dependent_to_base(vec1, q1) result(vec2)

        implicit none
        
        real, intent(in) :: vec1(3), q1(4)
        real :: vec2(3)

        real :: t(4)

        !t(1) = q1(1)
        !t(2) = -q1(2)
        !t(3) = -q1(3)
        !t(4) = -q1(4)

        !vec2 = quat_base_to_dependent(vec1, t)

        t(1) =  vec1(1)*q1(2) + vec1(2)*q1(3) + vec1(3)*q1(4)
        t(2) =  vec1(1)*q1(1) - vec1(2)*q1(4) + vec1(3)*q1(3)
        t(3) =  vec1(1)*q1(4) + vec1(2)*q1(1) - vec1(3)*q1(2)
        t(4) = -vec1(1)*q1(3) + vec1(2)*q1(2) + vec1(3)*q1(1)

        vec2(1) = q1(1)*t(2) + q1(2)*t(1) + q1(3)*t(4) - q1(4)*t(3)
        vec2(2) = q1(1)*t(3) - q1(2)*t(4) + q1(3)*t(1) + q1(4)*t(2)
        vec2(3) = q1(1)*t(4) + q1(2)*t(3) - q1(3)*t(2) + q1(4)*t(1)
        
    end function quat_dependent_to_base


    subroutine quat_norm(q1)

        implicit none
        
        real, intent(inout) :: q1(4)

        q1 = q1/sqrt(q1(1)**2 + q1(2)**2 + q1(3)**2 + q1(4)**2)

    end subroutine quat_norm


    function euler_to_quat(e1) result(q1)
        implicit none
        
        real, intent(in) :: e1(3)
        real :: q1(4)

        !q1(1) = cos(e1(1)*0.5)*cos(e1(2)*0.5)*cos(e1(3)*0.5) + sin(e1(1)*0.5)*sin(e1(2)*0.5)*sin(e1(3)*0.5)
        !q1(2) = sin(e1(1)*0.5)*cos(e1(2)*0.5)*cos(e1(3)*0.5) - cos(e1(1)*0.5)*sin(e1(2)*0.5)*sin(e1(3)*0.5)
        !q1(3) = cos(e1(1)*0.5)*sin(e1(2)*0.5)*cos(e1(3)*0.5) + sin(e1(1)*0.5)*cos(e1(2)*0.5)*sin(e1(3)*0.5)
        !q1(4) = cos(e1(1)*0.5)*cos(e1(2)*0.5)*sin(e1(3)*0.5) - sin(e1(1)*0.5)*sin(e1(2)*0.5)*cos(e1(3)*0.5)

        !real :: s(3), c(3)

        !s(1) = sin(e1(1)*0.5)
        !s(2) = sin(e1(2)*0.5)
        !s(3) = sin(e1(3)*0.5)

        !c(1) = cos(e1(1)*0.5)
        !c(2) = cos(e1(2)*0.5)
        !c(3) = cos(e1(3)*0.5)

        !q1(1) = c(1)*c(2)*c(3) + s(1)*s(2)*s(3)    
        !q1(2) = s(1)*c(2)*c(3) - c(1)*s(2)*s(3)    
        !q1(3) = c(1)*s(2)*c(3) + s(1)*c(2)*s(3)    
        !q1(4) = c(1)*c(2)*s(3) - s(1)*s(2)*c(3)    

        real :: cphi, ctheta, cpsi, sphi, stheta, spsi

        sphi = sin(e1(1)*0.5)
        stheta = sin(e1(2)*0.5)
        spsi = sin(e1(3)*0.5)

        cphi = cos(e1(1)*0.5)
        ctheta = cos(e1(2)*0.5)
        cpsi = cos(e1(3)*0.5)

        q1(1) = cphi*ctheta*cpsi + sphi*stheta*spsi    
        q1(2) = sphi*ctheta*cpsi - cphi*stheta*spsi    
        q1(3) = cphi*stheta*cpsi + sphi*ctheta*spsi    
        q1(4) = cphi*ctheta*spsi - sphi*stheta*cpsi    
        
    end function euler_to_quat


    function quat_to_euler(q1) result(e1)
        implicit none
        
        real, intent(in) :: q1(4)
        real :: e1(3)

        ! CONSIDER WRITING IT INTO A SINGLE IF STATEMENT
        ! condition = e0*ey-ex*ez
        ! if (abs(abs(condition) - 0.5) < TOLERANCE) then
        !   if(condition > 0) then
        !       ...
        !   else
        !       ...
        !   end if
        ! else
        !   ...
        ! end if

        if (abs(q1(1)*q1(3) - q1(2)*q1(4) - 0.5) < 1e-12) then
            
            e1(1) = 2*asin(q1(2)/cos(PI4))
            e1(2) = PI2
            e1(3) = 0.0

        else if (abs(q1(1)*q1(3) - q1(2)*q1(4) + 0.5) < 1e-12) then

            e1(1) = 2*asin(q1(2)/cos(PI4))
            e1(2) = -PI2
            e1(3) = 0.0

        else

            e1(1) = atan2(2*(q1(1)*q1(2) + q1(3)*q1(4)), (q1(1)**2 + q1(4)**2 - q1(2)**2 - q1(3)**2))
            e1(2) = asin(2*(q1(1)*q1(3) - q1(2)*q1(4)))
            e1(3) = atan2(2*(q1(1)*q1(4) + q1(2)*q1(3)), (q1(1)**2 + q1(2)**2 - q1(3)**2 - q1(4)**2))

            if (e1(3) < 0.0) e1(3) = e1(3) + 2.0*pi

        end if
    
        
    end function quat_to_euler

    !!! STD ATMOSPHERE !!!

    function gravity_SI(h) result(g)
        implicit none
        
        real, intent(in) :: h 
        real :: g

        g = g_ssl*(r_ez/(r_ez + h))**2
        
    end function gravity_SI

    function gravity_English(h_ft) result(g)
        implicit none
        
        real, intent(in) :: h_ft
        real :: g

        g = meters_to_feet(gravity_SI(feet_to_meters(h_ft)))
        
    end function gravity_English


    subroutine std_atm_SI(h, Z, T, P, rho, a)

        implicit none

        real, intent(in) :: h
        real, intent(inout) :: Z, T, P, rho, a
        real :: h_p

        integer :: i


        if (h < 0.0) then
            h_p = 0.0
        else
            h_p = h
        end if

        ! Calculate geopotential altitude
        Z = r_ez*h_p/(r_ez + h_p)

        ! Calculate pressure
        do i = 1, 8

            if (Z >= Z_i(i) .and. Z < Z_i1(i)) then

                T = T_i(i) + 0.001*Tp_i(i)*(Z - Z_i(i))

                if (Tp_i(i) == 0) then
                    P = P_i(i)*exp(-g_ssl*(Z - Z_i(i))/(R_gas * T_i(i)))
                else
                    P = P_i(i)*(T/T_i(i))**(-g_ssl/(R_gas * 0.001*Tp_i(i)))
                end if 

                rho = P/(R_gas*T)

                a = sqrt(gamma*R_gas*T)

                exit
            else
                continue
            end if

        end do

    end subroutine std_atm_SI

    ! FIXED VERSION
    !subroutine std_atm_SI(h, Z, T, P, rho, a)

        !implicit none

        !real, intent(in) :: h
        !real, intent(inout) :: Z, T, P, rho, a
        !real :: h_p, P_0

        !integer :: i

        !P_0 = P_ssl

        !if (h < 0.0) then
            !h_p = 0.0
        !else
            !h_p = h
        !end if

        !! Calculate geopotential altitude
        !Z = r_ez*h_p/(r_ez + h_p)

        !! Calculate pressure
        !do i = 1, 8

            !if (Z >= Z_i(i)) then

                !!T = T_i(i) + 0.001*Tp_i(i)*(Z - Z_i(i))

                !if (Tp_i(i) == 0) then
                    !if (Z < Z_i1(i)) then
                        !T = T_i(i)
                        !P = P_0*exp(-g_ssl*(Z - Z_i(i))/(R_gas * T_i(i)))
                    !else
                        !P_0 = P_0*exp(-g_ssl*(Z_i1(i) - Z_i(i))/(R_gas * T_i(i)))
                    !end if
                !else
                    !if (Z < Z_i1(i)) then
                        !T = T_i(i) + 0.001*Tp_i(i)*(Z - Z_i(i))
                        !P = P_0*(T/T_i(i))**(-g_ssl/(R_gas * 0.001*Tp_i(i)))
                    !else 
                        !P_0 = P_0*((T_i(i) + 0.001*Tp_i(i)*(Z_i1(i) - Z_i(i)))/T_i(i))**(-g_ssl/(R_gas * 0.001*Tp_i(i)))
                    !end if
                !end if 

            !end if

        !end do
        !rho = P/(R_gas*T)
        !a = sqrt(gamma*R_gas*T)

    !end subroutine std_atm_SI

    subroutine std_atm_English(h, Z, T, P, rho, a)

        implicit none
        
        real, intent(in) :: h
        real, intent(inout) :: Z, T, P, rho, a

        real :: h_m
        
        h_m = feet_to_meters(h)

        call std_atm_SI(h_m, Z, T, P, rho, a)

        Z = meters_to_feet(Z)
        T = kelvin_to_rankine(T)
        P = P*1.0/47.880258
        rho = rho*kg_m3_to_slug_ft3
        a = meters_to_feet(a)

    end subroutine std_atm_English

    ! ADD PRINT ATMOSPHERE

    function sutherland_visc_SI(T) result(mu)

        implicit none
        
        real, intent(in) :: T
        real :: mu

        mu = 1.716e-05*((273.15+110.4)/(T + 110.4))*(T/273.15)**(3./2.)

    end function sutherland_visc_SI

    function sutherland_visc_English(T) result(mu)

        implicit none
        
        real, intent(in) :: T
        real :: mu, T_K

        T_K = rankine_to_kelvin(T)

        mu = sutherland_visc_SI(T_K)*1.0/47.880258

    end function sutherland_visc_English


    !!! UNIT CONVERSIONS !!!

    function feet_to_meters(l_ft) result(l_m)

        implicit none
        
        real, intent(in) :: l_ft
        real :: l_m

        l_m = 0.3048*l_ft

    end function feet_to_meters

    function meters_to_feet(l_m) result(l_ft)

        implicit none
        
        real, intent(in) :: l_m
        real :: l_ft

        l_ft = 3.280839895013123*l_m


    end function meters_to_feet

    function kelvin_to_rankine(T) result(T_R)

        implicit none
        
        real, intent(in) :: T
        real :: T_R

        T_R = 1.8*T

    end function kelvin_to_rankine

    function rankine_to_kelvin(T) result(T_K)
        implicit none
        
        real, intent(in) :: T
        real :: T_K
    
        T_K = T*1.0/1.8

    end function rankine_to_kelvin

    function y_English_to_SI(y_eng) result(y_SI)

        implicit none

        real, intent(in) :: y_eng(13)
        real :: y_SI(13)

        y_SI(1:3) = y_eng(1:3)*0.3048
        y_SI(4:6) = y_eng(4:6)
        y_SI(7:9) = y_eng(7:9)*0.3048
        y_SI(10:13) = y_eng(10:13)

    end function y_English_to_SI

    function y_SI_to_English(y_SI) result(y_eng)

        implicit none

        real, intent(in) :: y_SI(13)
        real :: y_eng(13)

        y_eng(1:3) = y_SI(1:3)*1.0/0.3048
        y_eng(4:6) = y_SI(4:6)
        y_eng(7:9) = y_SI(7:9)*1.0/0.3048
        y_eng(10:13) = y_SI(10:13)

    end function y_SI_to_English


    ! Turbulence Functions
    !---------------------------------------------------------------------------------
    function rand_normal() result(z)
        !! Return a single normally distributed random number (mean=0, std=1)
        real :: z
        real, save :: z2 = 0.0
        logical, save :: has_saved = .false.
        real :: u1, u2, r, theta
        
        if (has_saved) then
            !! Return the saved value from last time
            z = z2
            has_saved = .false.
        else
            !! Generate two uniforms
            call random_number(u1)
            call random_number(u2)
            
            !! Guard against u1 = 0
            if (u1 <= 0.0) u1 = 1.0e-12
            if (u1 >  1.0) u1 = 1.0
            
            !! Box-Muller transform
            r = sqrt(-2.0*log(u1))
            theta = 2.0*pi*u2
            z  = r * cos(theta)
            z2 = r * sin(theta)
            
            has_saved = .true.
        end if
    end function rand_normal

    !---------------------------------------------------------------------------------
    function   interpolate_1D(xa, ya, x, saturate) result(y)
        real, intent(in) :: xa(:), ya(:), x
        logical, intent(in), optional :: saturate
        real :: y
        integer :: i, N
        logical :: sat
        N = size(xa)
        if (N /= size(ya)) then
            write(*,*) 'Error interpolating 1D array: sizes do not match!'
            stop
        end if
        if (present(saturate)) then
            sat = saturate
        else
            sat = .true.
        end if
        if (sat) then
            if (x <= xa(1)) then
                y = ya(1)
                return
            end if
            do i=2, N
                if (xa(i-1) <= x .and. x <= xa(i)) then
                    y = interpolate(xa(i-1), xa(i), x, ya(i-1), ya(i))
                    return
                end if
            end do
            y = ya(N)
        else
            do i=2, N
                if (xa(i-1) <= x .and. x <= xa(i)) then
                    y = interpolate(xa(i-1), xa(i), x, ya(i-1), ya(i))
                    return
                end if
            end do
            write(*,*) 'Error interpolating 1D array: value outside of the range'
            stop
        end if
    end function interpolate_1D

    !---------------------------------------------------------------------------------
    function interpolate(x1, x2, x, y1, y2) result(y)
        real, intent(in) :: x1, x2, x, y1, y2
        real :: y
        y = (y2-y1) / (x2-x1) * (x-x1) + y1
    end function interpolate

    !---------------------------------------------------------------------------------
    subroutine psd(x, dt, psd_norm, stdev, filename)
        implicit none
        real, intent(in) :: x(:), dt
        real, intent(out), optional :: psd_norm(:,:)
        real, intent(out), optional :: stdev
        character(len=*), intent(in), optional :: filename
        integer :: n, k, i, iunit
        real :: fs, f, mean
        real, allocatable :: Pxx(:), loc_psd_norm(:,:)
        complex :: eye, temp

        n = size(x)
        if (mod(n,2) /= 0) error stop "periodogram: N must be even for N/2 bin."
        allocate(Pxx(n/2+1))
        allocate(loc_psd_norm(n/2+1,2))
        eye = cmplx(0.0, 1.0)

        ! compute mean and stdev
        mean = 0.0
        stdev = 0.0
        do i=1,n
            mean = mean + x(i)/n
        end do
        do i=1,n
            stdev = stdev + (x(i)-mean)**2
        end do
        stdev = sqrt(stdev/(n-1))

        if(present(filename)) then
            open(newunit=iunit, file=filename, status="replace", action="write")
            write(iunit,'(A)') 'fs[Hz],PSD[ft^2/Hz],Normalized_PSD,Mean,STDEV'
            write(*,*) trim(filename),'   Mean = ',mean,'   Std. Dev = ',stdev
        end if

        Pxx(:) = 0.0
        fs = 1.0/dt
        do k=0,n/2
            f = real(k)*fs/real(n)
            temp = cmplx(0.0, 0.0)
            do i=1,n
                temp = temp + (x(i)-mean)*exp(-eye*2.0*PI*real(i-1)*real(k)/real(n))
            end do
            Pxx(k+1) = 1/fs/real(n)*abs(temp)**2
            if (k /= 0 .and. k /= n/2) Pxx(k+1) = 2.0*Pxx(k+1)
            loc_psd_norm(k+1,1) = f
            loc_psd_norm(k+1,2) = Pxx(k+1)/stdev**2
            if(present(filename)) then
                if(k==0) then
                    write(iunit,'(ES20.12,",",ES20.12,",",ES20.12,",",ES20.12,",",ES20.12)') f, Pxx(k+1), loc_psd_norm(k+1,2),mean,stdev
                else
                    write(iunit,'(ES20.12,",",ES20.12,",",ES20.12)') f, Pxx(k+1), loc_psd_norm(k+1,2)
                end if
            end if
        end do

        if(present(filename)) close(iunit)
        if(present(psd_norm)) psd_norm = loc_psd_norm
    end subroutine psd

    !---------------------------------------------------------------------------------
    subroutine test_rand_normal()
        implicit none
        integer, parameter :: n = 100000
        integer, parameter :: nbins = 100
        real, allocatable :: x(:), hist(:)
        real :: mean, sigma, xmin, xmax, dx, binval
        integer :: i, b, iunit
        character(len=:), allocatable :: filename  

        write(*,*) '   Testing function rand_normal():'

        allocate(x(n))
        allocate(hist(nbins))
        do i = 1, n
            ! call random_number(x(i)) ! this is a uniform distribution from 0 to 1 with a mean of 0.5
            x(i) = rand_normal() ! this is a normal distribution with a mean of 0 and a std dev of 1
        end do

        mean  = sum(x) / n
        sigma = sqrt( sum((x-mean)**2) / (n-1) )

        write(*,*) '     Mean  = ', mean
        write(*,*) '     Std   = ', sigma
        write(*,*) '     Expected mean error ~ ', 1.0/sqrt(real(n))
        write(*,*) '     Expected std error  ~ ', 1.0/sqrt(2.0*n)

        ! Build histogram
        xmin = -5.0; xmax = 5.0
        dx = (xmax-xmin)/nbins
        hist = 0.0

        do i = 1, n
            if (x(i) >= xmin .and. x(i) < xmax) then
                b = int((x(i)-xmin)/dx) + 1
                hist(b) = hist(b) + 1
            end if
        end do

        hist = hist / (n*dx)   ! normalize to PDF

        filename = 'rand_normal_test.csv'
        open(newunit=iunit, file=filename, status="replace", action="write")
        write(iunit,'(A)') 'bin_number,bin_midpoint_value,normalized_histogram,analytic_solution'
        do i=1,nbins
            binval = xmin+dx*i-0.5*dx
            write(iunit,*) i,',',binval,',',hist(i),',',exp(-0.5*(binval**2))/sqrt(2.0*PI)
        end do
        write(*,*) '     Histogram written to file: ',trim(filename)
        write(*,*)
        close(iunit)

        ! Compute PSD
        ! call psd(x(:),0.01,filename='rand_normal_psd.csv')

        deallocate(x)
        deallocate(hist)
    end subroutine test_rand_normal


end module goates_m
