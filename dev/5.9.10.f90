module simulation_m

    use goates_m
    use jsonx_m

    implicit none
    
    real, parameter :: one_sixth = 1./6.
contains
    
    !!! INTEGRATION
    function runge_kutta(t_0, y_0, dt) result(y)

        implicit none
        
        real, intent(in) :: t_0, y_0(13), dt
        real :: y(13)

        real :: k1(13), k2(13), k3(13), k4(13)


        k1 = differential_equations(t_0, y_0)
        k2 = differential_equations(t_0 + 0.5*dt, y_0 + k1*0.5*dt)
        k3 = differential_equations(t_0 + 0.5*dt, y_0 + k2*0.5*dt)
        k4 = differential_equations(t_0 + dt, y_0 + k3*dt)

        y = y_0 + (dt*one_sixth)*(k1 + 2*k2 + 2*k3 + k4)

    end function runge_kutta

    function differential_equations(t, y) result(dy_dt)

        implicit none
        
        real, intent(in) :: t, y(13)
        real :: dy_dt(13)

        real :: mass, I(3,3), F(3), M(3), g, I_inv(3,3), dummy(3)

        g = gravity_English(-y(9))

        call mass_inertia(t, y, mass, I)
        call psuedo_aerodynamics(t, y, F, M)

        dy_dt = 0.0

        ! Sim of Flight Eq. 5.4.5
        
        dy_dt(1:3) = (1.0/mass)*F
        !write(*,*) "F: ", F
        !write(*,*) "dy_dt: ", dy_dt
        dy_dt(1) = dy_dt(1) + g*(2*(y(11)*y(13) - y(12)*y(10))) + (y(6)*y(2) - y(5)*y(3))
        dy_dt(2) = dy_dt(2) + g*(2*(y(12)*y(13) + y(11)*y(10))) + (y(4)*y(3) - y(6)*y(1))
        dy_dt(3) = dy_dt(3) + g*(y(13)**2 + y(10)**2 - y(11)**2 - y(12)**2) + (y(5)*y(1) - y(4)*y(2))
        !write(*,*) "dy_dt: ", dy_dt

        ! Calculate I_inv
        I_inv(1,1) = I(2,2)*I(3,3) - I(2,3)*I(3,2)
        I_inv(1,2) = I(1,3)*I(3,2) - I(1,2)*I(3,3)
        I_inv(1,3) = I(1,2)*I(2,3) - I(1,3)*I(2,2)
        I_inv(2,1) = I(2,3)*I(3,1) - I(2,1)*I(3,3)
        I_inv(2,2) = I(1,1)*I(3,3) - I(1,3)*I(3,1)
        I_inv(2,3) = I(1,3)*I(2,1) - I(1,1)*I(2,3)
        I_inv(3,1) = I(2,1)*I(3,2) - I(2,2)*I(3,1)
        I_inv(3,2) = I(1,2)*I(3,1) - I(1,1)*I(3,2)
        I_inv(3,3) = I(1,1)*I(2,2) - I(1,2)*I(2,1)

        I_inv = I_inv/(I(1,1)*(I(2,2)*I(3,3) - I(2,3)*I(3,2)) - I(1,2)*(I(2,1)*I(3,3) - I(2,3)*I(3,1)) &
                            + I(1,3)*(I(2,1)*I(3,2) - I(2,2)*I(3,1)))

        ! Eq. 5.4.6
        dummy(1) = M(1) + (I(2,2) - I(3,3))*y(5)*y(6) + I(2,3)*(y(5)**2 - y(6)**2) + I(1,3)*y(4)*y(5) - I(1,2)*y(4)*y(6)
        dummy(2) = M(2) + (I(3,3) - I(1,1))*y(4)*y(6) + I(1,3)*(y(6)**2 - y(4)**2) + I(1,2)*y(5)*y(6) - I(2,3)*y(4)*y(5)
        dummy(3) = M(3) + (I(1,1) - I(2,2))*y(4)*y(5) + I(1,2)*(y(4)**2 - y(5)**2) + I(2,3)*y(4)*y(6) - I(1,3)*y(5)*y(6)
        dy_dt(4:6) = matmul(I_inv, dummy)

        ! Eq. 5.4.7
        dy_dt(7:9) = quat_dependent_to_base(y(1:3), y(10:13))
        write(*,*) "TEST: ", (2*(y(11)*y(12) + y(13)*y(10))*y(1) &
                    + (y(12)**2 + y(10)**2 - y(11)**2 - y(13)**2)*y(2) + 2*(y(12)*y(13) - y(11)*y(10))*y(3))

        ! Eq. 5.4.8
        dy_dt(10) = 0.5*(- y(11)*y(4) - y(12)*y(5) - y(13)*y(6))
        dy_dt(11) = 0.5*(  y(10)*y(4) - y(13)*y(5) + y(12)*y(6))
        dy_dt(12) = 0.5*(  y(13)*y(4) + y(10)*y(5) - y(11)*y(6))
        dy_dt(13) = 0.5*(- y(12)*y(4) + y(11)*y(5) + y(10)*y(6))
        
        write(*,*) "dy_dt: ", dy_dt
        write(*,*) "----------------------"


    end function differential_equations

    subroutine psuedo_aerodynamics(t, y, F, M)

        implicit none

        real, intent(in) :: t, y(13)
        real, intent(out) :: F(3), M(3)

        real :: Z, Temp, P, rho, a
        real :: C_L, C_D, C_S, C_ell, C_m, C_n
        real :: C_Lalpha, C_D0, C_D2, C_malpha, C_mqbar, C_ell0, C_ellpbar
        real :: alpha, beta, pbar, qbar, rbar, V, S_w, b, c
        real :: S_alpha, C_alpha, S_beta, C_beta

        print*, "Time: ", t
        print*, "State vector incoming: ", y

        ! Get atmosphere
        call std_atm_English(-y(9), Z, temp, P, rho, a)
        
        !! ARROW CONSTANTS
        C_Lalpha = 4.929
        C_D0 = 5.096
        C_D2 = 48.138
        C_malpha = -2.605
        C_mqbar = -9.06
        C_ellpbar = -5.378
        ! Angled Fletchings
        !C_ell0 = 3.223
        ! Straight Fletchings
        C_ell0 = 0.0

        S_w = 0.000218 ![ft^2]
        b = 2.3 ![ft]
        c = 2.3 ![ft]

        V = sqrt(y(1)**2 + y(2)**2 + y(3)**2)

        pbar = (0.5/V)*y(4)*b
        qbar = (0.5/V)*y(5)*c
        rbar = (0.5/V)*y(6)*b

        alpha = atan2(y(3),y(1))
        beta = atan2(y(2),y(1))

        C_L = C_Lalpha*alpha
        C_S = -C_Lalpha*beta
        C_D = C_D0 + C_D2 * C_L**2 + C_D2 * C_S**2
        C_ell = C_ell0 + C_ellpbar*pbar
        C_m = C_malpha*alpha + C_mqbar*qbar
        C_n = -C_malpha*beta + C_mqbar*rbar

        beta = asin(y(2)/V)
        S_alpha = sin(alpha)
        C_alpha = cos(alpha)
        S_beta = sin(beta)
        C_beta = cos(beta)


        F(1) = 0.5*rho*V**2 * S_w * (C_L*S_alpha - C_S*C_alpha*S_beta -C_D*C_alpha*C_beta)
        F(2) = 0.5*rho*V**2 * S_w * (C_S*C_beta - C_D*S_beta)
        F(3) = 0.5*rho*V**2 * S_w * (-C_L*C_alpha - C_S*S_alpha*S_beta - C_D*S_alpha*C_beta)

        M(1) = 0.5*rho*V**2 * S_w * (b*(C_ell))
        M(2) = 0.5*rho*V**2 * S_w * (c*C_m)
        M(3) = 0.5*rho*V**2 * S_w * (b*(C_n))


        write(*,*) "F: ", F
        write(*,*) "M: ", M

    end subroutine psuedo_aerodynamics


    subroutine mass_inertia(t, y, M, I)

        implicit none

        real, intent(in) :: t, y(13)
        real, intent(out) :: M
        real, intent(out) :: I(3,3)

        real :: r1, r2
        integer :: j

        M = 0.0697 ! lbf at sea level
        M = M*0.3048/g_ssl ! slugs

        I = 0.0
        I(1,1) = 0.0000194 ! slugs-ft^2
        I(2,2) = 0.00097   ! slugs-ft^2
        I(3,3) = 0.00097   ! slugs-ft^2

    end subroutine mass_inertia

    subroutine simulation_main(dt, y_0)

        implicit none

        real, intent(in) :: dt, y_0(13)

        real :: t, y(13)
        integer :: io_unit

        t = 0.0
        y = y_0

        open(newunit=io_unit, file='5.9.10_output.csv', status='replace', action='write')
        write(io_unit,*) 'time[s], u[ft/s], v[ft/s], w[ft/s], p[rad/s], q[rad/s], r[rad/s], xf[ft], yf[ft], zf[ft], e0, ex, ey, ez'
        write(io_unit,*) t,',',y(1),',',y(2),',',y(3),',',y(4),',',y(5),',' &
                            ,y(6),',',y(7),',',y(8),',',y(9),',',y(10),',',y(11),',',y(12),',',y(13)

        call quat_norm(y(10:13))

        ! Run until the arrow reaches the ground
        do while (y(9) < 0.0)

            write(*,*) "time: ", t, " s"
            write(*,*) y
            write(*,*) "----------------------"

            y = runge_kutta(t, y, dt)
            t = t + dt
            write(*,*) "----------------------"
            write(*,*) "----------------------"

            call quat_norm(y(10:13))

            write(io_unit,*) t,',',y(1),',',y(2),',',y(3),',',y(4),',',y(5),',' &
                            ,y(6),',',y(7),',',y(8),',',y(9),',',y(10),',',y(11),',',y(12),',',y(13)

        end do

        close(io_unit)


    end subroutine simulation_main

end module simulation_m 

program main
    use simulation_m
    implicit none


    character(len=100) :: input_file
    type(json_value), pointer :: j_main
    real :: dt, V, H, theta, phi


    real :: y(13)

    ! Get input file from command line
    call get_command_argument(1, input_file)

    ! Load JSON file
    call jsonx_load(input_file, j_main)

    call jsonx_get(j_main, "simulation.time_step[s]", dt, default_value=0.01)
    call jsonx_get(j_main, "initial.airspeed[ft/s]", V)
    call jsonx_get(j_main, "initial.altitude[ft]", H)
    call jsonx_get(j_main, "initial.elevation_angle[deg]", theta)
    call jsonx_get(j_main, "initial.bank_angle[deg]", phi)

    theta = theta*PI/180.
    phi = phi*PI/180.

    y = 0.0
    y(1) = V
    
    y(9) = -H

    y(10:13) = euler_to_quat([phi, theta, 0.0])

    call simulation_main(dt, y)

end program main
