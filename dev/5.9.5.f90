module simulation_m

    use goates_m

    implicit none
    
    real, parameter :: one_sixth = 1./6.
contains
    
    !!! INTEGRATION
    function runge_kutta(t_0, y_0, dt) result(y)

        implicit none
        
        real, intent(in) :: t_0, y_0(13), dt
        real :: y(13)

        real :: k1(13), k2(13), k3(13), k4(13)

        call quat_norm(y(10:13))

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

        real :: mass, I(3,3), F(3), M(3), g, I_inv(3,3)

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
        I_inv(1,2) = -(I(1,2)*I(3,3) - I(1,3)*I(3,2))
        I_inv(1,3) = (I(1,2)*I(2,3) - I(1,3)*I(2,2))
        I_inv(2,1) = -(I(2,1)*I(3,3) - I(2,3)*I(3,1))
        I_inv(2,2) = -(I(1,1)*I(3,3) - I(1,3)*I(3,1))
        I_inv(2,3) = -(I(1,1)*I(2,3) - I(1,3)*I(2,1))
        I_inv(3,1) = (I(2,1)*I(3,2) - I(2,2)*I(3,1))
        I_inv(3,2) = -(I(1,1)*I(3,2) - I(1,2)*I(3,1))
        I_inv(3,3) = (I(1,1)*I(2,2) - I(1,2)*I(2,1))

        I_inv = I_inv/(I(1,1)*(I(2,2)*I(3,3) - I(2,3)*I(3,2)) - I(1,2)*(I(2,1)*I(3,3) - I(2,3)*I(3,1)) &
                            + I(1,3)*(I(2,1)*I(3,2) - I(2,2)*I(3,1)))

        ! Eq. 5.4.6
        ! TO BE IMPLEMENTED LATER


        ! Eq. 5.4.7
        dy_dt(7:9) = quat_dependent_to_base(y(1:3), y(10:13))

        ! Eq. 5.4.8
        dy_dt(10) = 0.5*(- y(11)*y(4) - y(12)*y(5) - y(13)*y(6))
        dy_dt(11) = 0.5*(  y(10)*y(4) - y(13)*y(5) + y(12)*y(6))
        dy_dt(12) = 0.5*(  y(13)*y(4) + y(10)*y(5) - y(11)*y(6))
        dy_dt(13) = 0.5*(- y(12)*y(4) + y(11)*y(5) + y(10)*y(6))
        


    end function differential_equations

    subroutine psuedo_aerodynamics(t, y, F, M)

        implicit none

        real, intent(in) :: t, y(13)
        real, intent(out) :: F(3), M(3)

        real :: Z, Temp, P, rho, a, mu
        real :: radius, Re, CD, V_c, u_c(3)

        ! Sphere properties
        ! TENNIS BALL
        radius = 0.21667/2.0 ! feet

        ! PING PONG BALL
        !radius = 0.13084/2.0 ! feet

        ! Get air density
        call std_atm_English(-y(9), Z, Temp, P, rho, a)

        M = 0.0

        ! Calculate drag force
        mu = sutherland_visc_English(Temp)

        Re = 2.0*rho*norm2(y(1:3))*radius/mu

        if (Re < 0.01) then
            CD = 2405.0
        else if (Re >= 0.01 .and. Re <= 450000.) then
            CD = 24.0/Re + 6.0/(1.0+sqrt(Re)) + 0.4
        else if (Re > 450000. .and. Re <= 560000) then
            CD = 1.0e29*Re**(-5.211)
        else if (Re > 560000. .and. Re <= 14000000.) then
            CD = -2.0e-23*Re**3 - 1.e-16*Re**2 + 9.0e-9*Re + 0.069
        else
            CD = 0.12
        end if

        V_c = norm2(y(1:3))
        u_c = y(1:3)/V_c

        F = -0.5 * rho * V_c**2 * PI * radius**2 * CD * u_c
        
    end subroutine psuedo_aerodynamics


    subroutine mass_inertia(t, y, M, I)

        implicit none

        real, intent(in) :: t, y(13)
        real, intent(out) :: M
        real, intent(out) :: I(3,3)

        real :: r1, r2
        integer :: j

        ! TENNIS BALL
        r2 = 0.21667/2.0
        r1 = r2 - 0.00984252

        M = 0.125 ! lbf

        ! PING PONG BALL
        !r2 = 0.13084/2.0
        !r1 = r2 - 0.00131

        !M = 0.006 ! lbf

        ! Convert
        M = M*0.3048/g_ssl ! slugs

        I = 0.0
        forall(j = 1:3) I(j,j) = 1.0

        I = I*M*(2*(r2**5 - r1**5))/(5*(r2**3 - r1**3))

        !M = M*1.0*14.5939029

        !I = I*1.0*14.5939029
        !I = I*0.3048**2

    end subroutine mass_inertia

    subroutine simulation_main(t_0, t_f, dt, y_0)

        implicit none

        real, intent(in) :: t_0, t_f, dt, y_0(13)

        real :: t, y(13)
        logical :: check
        integer :: io_unit

        t = t_0
        y = y_0

        ! TENNIS BALL
        open(newunit=io_unit, file='tennis_serve.csv', status='replace', action='write')
        write(io_unit,*) 'time[s], u[ft/s], v[ft/s], w[ft/s], p[rad/s], q[rad/s], r[rad/s], xf[ft], yf[ft], zf[ft], e0, ex, ey, ez'
        write(io_unit,*) t,',',y(1),',',y(2),',',y(3),',',y(4),',',y(5),',' &
                            ,y(6),',',y(7),',',y(8),',',y(9),',',y(10),',',y(11),',',y(12),',',y(13)

        check = .false.
        do while (y(9) < 0.0)
            write(*,*) "time: ", t, " s"
            write(*,*) y
            write(*,*) "----------------------"

            y = runge_kutta(t,y,dt)
            t = t+dt
            write(io_unit,*) t,',',y(1),',',y(2),',',y(3),',',y(4),',',y(5),',' &
                            ,y(6),',',y(7),',',y(8),',',y(9),',',y(10),',',y(11),',',y(12),',',y(13)

            if (check) then
                if (y(7) > 39.0) exit
            end if

            if (y(9) > -3.0) then
                if (y(7) < 39.0) check = .true.
            end if

        end do
        close(io_unit)

        ! PING PONG BALL
        !do while (t <= t_f)

            !write(*,*) "time: ", t, " s"
            !write(*,*) y
            !write(*,*) "----------------------"
            !y = runge_kutta(t, y, dt)
            !t = t + dt
            !write(*,*) "----------------------"
            !write(*,*) "----------------------"

            !call quat_norm(y(10:13))

        !end do



    end subroutine simulation_main

end module simulation_m 

program main
    use simulation_m
    implicit none

    real :: y(13)
    real :: M, I(3,3), theta
    y = 0.0

    theta = -7.5*PI/180.0

    ! Tennis Ball
    ! First serve
    !y(1) = 178.933
    ! Second serve
    y(1) = 140.8
    y(9) = -9.0
    y(10:13) = euler_to_quat([0.0, theta, 0.0])


    ! Ping Pong Ball
    !y(1) = 50.0
    !y(9) = -200.0
    !y(10) = 1.0
    
    call simulation_main(0.0, 10.0, 0.005, y)

end program main
