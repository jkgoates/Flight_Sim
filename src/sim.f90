module simulation_m

    use aircraft_m
    use jsonx_m
    use linalg_mod

    implicit none
    
    type(aircraft) :: vehicle
    ! Private control parameters
    real, private :: da, de, dr, throttle
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

        real :: mass, I(3,3), F(3), M(3), g, I_inv(3,3), dummy(3), h(3)
        real :: controls(4)

        if (verbose) then
            write(*,*) "t: ", t
            write(*,*) "y: ", y
        end if

        g = gravity_English(-y(9))

        controls = get_controls(t, y)

        call vehicle%mass_inertia(t, y, mass, I)
        call vehicle%aerodynamics(t, y, F, M, controls)
        call vehicle%gyroscopic(t, y, h)

        dy_dt = 0.0

        ! Sim of Flight Eq. 5.4.5
        
        dy_dt(1:3) = (1.0/mass)*F
        dy_dt(1) = dy_dt(1) + g*(2*(y(11)*y(13) - y(12)*y(10))) + (y(6)*y(2) - y(5)*y(3))
        dy_dt(2) = dy_dt(2) + g*(2*(y(12)*y(13) + y(11)*y(10))) + (y(4)*y(3) - y(6)*y(1))
        dy_dt(3) = dy_dt(3) + g*(y(13)**2 + y(10)**2 - y(11)**2 - y(12)**2) + (y(5)*y(1) - y(4)*y(2))

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
        dummy(1) = M(1) + (-h(3)*y(5) + h(2)*y(6)) 
        dummy(2) = M(2) + ( h(3)*y(4) - h(1)*y(6)) 
        dummy(3) = M(3) + (-h(2)*y(4) + h(1)*y(5)) 
        dummy(1) = dummy(1) + (I(2,2) - I(3,3))*y(5)*y(6) - I(2,3)*(y(5)**2 - y(6)**2) - I(1,3)*y(4)*y(5) + I(1,2)*y(4)*y(6)
        dummy(2) = dummy(2) + (I(3,3) - I(1,1))*y(4)*y(6) - I(1,3)*(y(6)**2 - y(4)**2) - I(1,2)*y(5)*y(6) + I(2,3)*y(4)*y(5)
        dummy(3) = dummy(3) + (I(1,1) - I(2,2))*y(4)*y(5) - I(1,2)*(y(4)**2 - y(5)**2) - I(2,3)*y(4)*y(6) + I(1,3)*y(5)*y(6)
        dy_dt(4:6) = matmul(I_inv, dummy)


        ! Eq. 5.4.7
        dy_dt(7:9) = quat_dependent_to_base(y(1:3), y(10:13))


        ! Eq. 5.4.8
        dy_dt(10) = 0.5*(- y(11)*y(4) - y(12)*y(5) - y(13)*y(6))
        dy_dt(11) = 0.5*(  y(10)*y(4) - y(13)*y(5) + y(12)*y(6))
        dy_dt(12) = 0.5*(  y(13)*y(4) + y(10)*y(5) - y(11)*y(6))
        dy_dt(13) = 0.5*(- y(12)*y(4) + y(11)*y(5) + y(10)*y(6))
        
        if (verbose) then
            write(*,*) "dy_dt: ", dy_dt
            write(*,*) "----------------------"
        end if


    end function differential_equations

    function get_controls(t, y) result(controls)

        implicit none
        real, intent(in) :: t, y(13)
        real :: controls(4)

        ! Get control inputs
        controls(1) = da*PI/180.
        controls(2) = de*PI/180.
        controls(3) = dr*PI/180.
        controls(4) = throttle

    end function get_controls

    ! Trim Functions
    function calc_R(t, y, G) result(R)

        implicit none
        real, intent(in) :: y(13), t, G(6)
        real :: R(6)




    end function calc_R
        

    subroutine run(j_main)

        implicit none

        type(json_value), pointer, intent(in) :: j_main

        type(json_value), pointer :: j_aircraft

        real :: t, y(13), y_temp(13)
        integer :: io_unit
        real :: dt, tf, V, H, theta, phi, psi
        real :: alpha, beta, p, q, r
        logical :: t_R
        real :: t_p, t_c, start_time, end_time


        call jsonx_get(j_main, "simulation.timestep[s]", dt, default_value=0.01)
        call jsonx_get(j_main, "simulation.total_time[s]", tf)
        call jsonx_get(j_main, "simulation.verbose", verbose, default_value=.false.)
        call jsonx_get(j_main, "initial.airspeed[ft/s]", V)
        call jsonx_get(j_main, "initial.altitude[ft]", H)
        call jsonx_get(j_main, "initial.elevation_angle[deg]", theta, default_value=0.0)
        call jsonx_get(j_main, "initial.bank_angle[deg]", phi, default_value=0.0)
        call jsonx_get(j_main, "initial.heading_angle[deg]", psi, default_value=0.0)
        call jsonx_get(j_main, "initial.alpha[deg]", alpha, default_value=0.0)
        call jsonx_get(j_main, "initial.beta[deg]", beta, default_value=0.0)
        call jsonx_get(j_main, "initial.p[deg/s]", p, default_value=0.0)
        call jsonx_get(j_main, "initial.q[deg/s]", q, default_value=0.0)
        call jsonx_get(j_main, "initial.r[deg/s]", r, default_value =0.0)
        call jsonx_get(j_main, "initial.aileron[deg]", da, default_value=0.0)
        call jsonx_get(j_main, "initial.elevator[deg]", de, default_value=0.0)
        call jsonx_get(j_main, "initial.rudder[deg]", dr, default_value=0.0)
        call jsonx_get(j_main, "initial.throttle", throttle, default_value=0.0)

        ! Initialize the aircraft
        call jsonx_get(j_main, "vehicle", j_aircraft)
        call vehicle%init(j_aircraft)

        ! Set initial conditions

        phi = phi*PI/180.
        theta = theta*PI/180.
        psi = psi*PI/180.
        
        alpha = alpha*PI/180.
        beta = beta*PI/180.

        y = 0.0
        t = 0.0

        y(1) = V*cos(alpha)*cos(beta)
        y(2) = V*sin(beta)
        y(3) = V*sin(alpha)*cos(beta)


        y(4) = p*PI/180.
        y(5) = q*PI/180.
        y(6) = r*PI/180.

        y(9) = -H

        y(10:13) = euler_to_quat([phi, theta, psi])



        open(newunit=io_unit, file='sim_output.csv', status='replace', action='write')
        write(io_unit,*) 'time[s], u[ft/s], v[ft/s], w[ft/s], p[rad/s], q[rad/s], r[rad/s], xf[ft], yf[ft], zf[ft], e0, ex, ey, ez'
        write(io_unit,*) t,',',y(1),',',y(2),',',y(3),',',y(4),',',y(5),',' &
                            ,y(6),',',y(7),',',y(8),',',y(9),',',y(10),',',y(11),',',y(12),',',y(13)

        call quat_norm(y(10:13))

        ! Check for real time sim
        if (dt == 0.0) then
            t_R = .true.
        else
            t_R = .false.
        end if

        if (t_R) then
            y_temp = y
            call cpu_time(t_p)
            y = runge_kutta(t, y, 0.01)
            call cpu_time(t_c)
            dt = t_c - t_p
            y = y_temp
            t_p = t_c
        end if


        call cpu_time(start_time)
        ! Run simulation
        do while (t < tf)

            if (verbose) then
                write(*,*) "time: ", t, " s"
                write(*,*) y
                write(*,*) "----------------------"
            end if

            y = runge_kutta(t, y, dt)
            t = t + dt
            if (t_R) then
                call cpu_time(t_c)
                dt = t_c - t_p
                t_p = t_c
            end if
            if (verbose) then
                write(*,*) "----------------------"
                write(*,*) "----------------------"
            end if

            call quat_norm(y(10:13))

            write(io_unit,*) t,',',y(1),',',y(2),',',y(3),',',y(4),',',y(5),',' &
                            ,y(6),',',y(7),',',y(8),',',y(9),',',y(10),',',y(11),',',y(12),',',y(13)

        end do
        call cpu_time(end_time)

        close(io_unit)

        write(*,*) "Simulation Finished in ", end_time - start_time, " seconds."
        
        if (t_R) then
            write(*,*) "Real time simulation error: ", tf - (end_time - start_time), " seconds."
        end if

    end subroutine run 

end module simulation_m 