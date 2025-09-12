module runge_kutta_m
    implicit none
    
    real, parameter :: one_sixth = 1./6.
contains
    
    !!! INTEGRATION
    function runge_kutta(t_0, y_0, dt) result(y)

        implicit none
        
        real, intent(in) :: t_0, y_0, dt
        real :: y

        real :: k1, k2, k3, k4

        k1 = differential_equations(t_0, y_0)
        k2 = differential_equations(t_0 + 0.5*dt, y_0 + k1*0.5*dt)
        k3 = differential_equations(t_0 + 0.5*dt, y_0 + k2*0.5*dt)
        k4 = differential_equations(t_0 + dt, y_0 + k3*dt)

        y = y_0 + one_sixth*dt*(k1 + 2*k2 + 2*k3 + k4)

        write(*,*) t_0, y_0, k1, k2, k3, k4, y

    end function runge_kutta


    function differential_equations(t, y) result(dy_dt)

        implicit none
        
        real, intent(in) :: t, y
        real :: dy_dt

        dy_dt = 1 + tan(y)

    end function differential_equations

    subroutine test_main()

        implicit none
        
        real :: t, y, dt
        integer :: i

        t = 0.0
        y = 0.0
        dt = 0.025

        do i = 1, 11

            y = runge_kutta(t, y, dt)
            t = t + dt

        end do

    end subroutine test_main

end module runge_kutta_m

program main

    use runge_kutta_m

    implicit none

    call test_main()
    
end program main