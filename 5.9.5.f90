module simulation_m
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

        dy_dt(1) = t + y(2)**2 * sin(y(1))
        dy_dt(2) = t + y(1) * cos(y(2))

    end function differential_equations

    subroutine psuedo_aerodynamics(t, y, F, M)

        implicit none

        real, intent(in) :: t, y(13)
        real, intent(out) :: F(3), M(3)
        
    end subroutine psuedo_aerodynamics


    subroutine mass_inertia(t, y, M, I)

        implicit none

        real, intent(in) :: t, y(13)
        real, intent(out) :: M
        real, intent(out) :: I(3,3)
        
    end subroutine mass_inertia

    subroutine simulation_main(t_0, t_f, dt, y_0)

        implicit none

        real, intent(in) :: t_0, t_f, dt, y_0(13)

        real :: t, y(13)        
        integer :: i

        t = t_0
        y = y_0

        do while (t < t_f)

            y = runge_kutta(t, y, dt)
            t = t + dt

        end do



    end subroutine simulation_main

end module simulation_m 
