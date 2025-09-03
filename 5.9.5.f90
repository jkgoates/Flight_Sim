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

        real :: mass, I(3,3), F(3), M(3), g

        g = gravity_SI(y(9))

        ! Sim of Flight Eq. 5.4.5
        dy_dt(1:3) = (1/mass)*F
        dy_dt(1) = dy_dt(1) + g*(2*(y(11)*y(13) - y(12)*y(10))) + (y(6)*y(2) - y(5)*y(3))
        dy_dt(2) = dy_dt(2) + g*(2*(y(12)*y(13) + y(11)*y(10))) + (y(4)*y(3) - y(6)*y(1))
        dy_dt(2) = dy_dt(2) + g*(y(13)**2 + y(10)**2 - y(11)**2 - y(12)**2) + (y(5)*y(1) - y(4)*y(2))

        ! Eq. 5.4.6

        ! Eq. 5.4.7
        dy_dt(7:9) = quat_dependent_to_base(y(7:9), y(10:13))

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
        radius = 0.1 ! meters

        ! Get air density
        call std_atm_SI(y(9), Z, Temp, P, rho, a)

        M = 0.0

        ! Calculate drag force
        mu = 1.716e-05*((273.15+110.4)/(Temp + 110.4))*(Temp/273.15)**(3./2.)
        Re = 2*rho*norm2(y(1:3))*radius/mu

        if (Re < 0.01) then
            CD = 2405.0
        else if (Re >= 0.01 .and. Re <= 450000.) then
            CD = 24/Re + 6/(1+sqrt(Re)) + 0.4
        else if (Re > 450000. .and. Re <= 560000) then
            CD = 1.0e29*Re**(-5.211)
        else if (Re > 560000. .and. Re <= 14000000.) then
            CD = -2.0e-23*Re**3 - 1.e-16*Re**2 + 9.0e-9*Re + 0.069
        else if (Re > 14000000.) then
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

        r1 = 0.08
        r2 = 0.1
        M = 0.1

        I = 0.0
        forall(j = 1:3) I(j,j) = 1.0

        I = I*M*(2*(r2**5 - r1**5))/(5*(r2**3 - r1**3))
    
        
    end subroutine mass_inertia

    subroutine simulation_main(t_0, t_f, dt, y_0)

        implicit none

        real, intent(in) :: t_0, t_f, dt, y_0(13)

        real :: t, y(13)        

        t = t_0
        y = y_0

        do while (t < t_f)

            y = runge_kutta(t, y, dt)
            t = t + dt

            call quat_norm(y(10:13))

        end do



    end subroutine simulation_main

end module simulation_m 

program main
    use simulation_m
    implicit none

    real :: y(13)
    y = 0.0
    
    call simulation_main(0.0, 0.1, 0.01, y)

end program main
