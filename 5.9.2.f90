module runge_kutta_m
    implicit none
    
    real, parameter :: one_sixth = 1./6.
contains
    
    !!! INTEGRATION
    function runge_kutta(t_0, y_0, dt) result(y)

        implicit none
        
        real, intent(in) :: t_0, y_0(2), dt
        real :: y(2)

        real :: k1(2), k2(2), k3(2), k4(2)

        k1 = differential_equations(t_0, y_0)
        k2 = differential_equations(t_0 + 0.5*dt, y_0 + k1*0.5*dt)
        k3 = differential_equations(t_0 + 0.5*dt, y_0 + k2*0.5*dt)
        k4 = differential_equations(t_0 + dt, y_0 + k3*dt)

        y = y_0 + (dt*one_sixth)*(k1 + 2*k2 + 2*k3 + k4)

        !write(*,*) t_0, y_0(1), y_0(2)

    end function runge_kutta


    function differential_equations(t, y) result(dy_dt)

        implicit none
        
        real, intent(in) :: t, y(2)
        real :: dy_dt(2)

        dy_dt(1) = t + y(2)**2 * sin(y(1))
        dy_dt(2) = t + y(1) * cos(y(2))

    end function differential_equations

    subroutine test_main()

        implicit none
        
        real :: t, y(2), dt
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

    integer :: i
    real :: start_time, end_time

    call cpu_time(start_time)
    do i = 1, 1000000
        call test_main()
    end do
    call cpu_time(end_time)
    print *, 'RK time total [sec]:          ', end_time - start_time
    
end program main