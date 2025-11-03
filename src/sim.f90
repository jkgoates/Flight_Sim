module simulation_m

    use aircraft_m
    use jsonx_m
    use linalg_mod
    use micro_time_m

    implicit none
    
    type(aircraft) :: vehicle

    real, parameter :: one_sixth = 1./6.
    real :: controls(4)
    real :: y_init(13)
    type(json_value), pointer :: j_main

    ! NOTES FROM DR. HUNSAKERS CODE
    

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

        if (verbose) then
            write(*,*) "t: ", t
            write(*,'(A,13ES20.12)') "y: ", y
        end if

        g = gravity_English(-y(9))

        call vehicle%mass_inertia(t, y, mass, I)
        call vehicle%aerodynamics(t, y, F, M, controls)
        call vehicle%gyroscopic(t, y, h)

        dy_dt = 0.0

        ! Sim of Flight Eq. 5.4.5
        
        dy_dt(1:3) = (1.0/mass)*F
        dy_dt(1) = dy_dt(1) + g*(2*(y(11)*y(13) - y(12)*y(10))) + (y(6)*y(2) - y(5)*y(3))
        dy_dt(2) = dy_dt(2) + g*(2*(y(12)*y(13) + y(11)*y(10))) + (y(4)*y(3) - y(6)*y(1))
        !dy_dt(3) = dy_dt(3) + g*(y(13)**2 + y(10)**2 - y(11)**2 - y(12)**2) + (y(5)*y(1) - y(4)*y(2))
        dy_dt(3) = dy_dt(3) + g*(y(13)*y(13) + y(10)*y(10) - y(11)*y(11) - y(12)*y(12)) + (y(5)*y(1) - y(4)*y(2))

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
            write(*,'(A,13ES20.12)') "dy_dt: ", dy_dt
            write(*,*) "----------------------"
        end if

    end function differential_equations

    function calc_R(V, H, rot_rates, G, var, theta, psi, solve_bank) result(R)

        implicit none
        real, intent(in) :: V, H, G(6), rot_rates(3), var, theta, psi
        logical, intent(in) :: solve_bank
        real :: alpha, beta, phi
        real :: R(6)

        real :: y_temp(13), dy_dt(13)


        ! Parse G
        alpha = G(1)
        if (solve_bank) then
            beta = var
            phi = G(2)
        else
            beta = G(2)
            phi = var
        end if
        controls = G(3:6)

        ! Set state
        y_temp = 0.0

        y_temp(1) = V*cos(alpha)*cos(beta)
        y_temp(2) = V*sin(beta)
        y_temp(3) = V*sin(alpha)*cos(beta)
        y_temp(4) = rot_rates(1)
        y_temp(5) = rot_rates(2)
        y_temp(6) = rot_rates(3)

        y_temp(9) = -H

        y_temp(10:13) = euler_to_quat((/phi, theta, psi/))

        ! Run diff_eq
        dy_dt = differential_equations(0.0, y_temp)
        
        R = dy_dt(1:6)

    end function calc_R

    subroutine init_trim(V_mag, H)

        implicit none
        
        real, intent(in) :: V_mag, H
        real :: y(13)

        real :: fd_step, relaxation, tol
        integer :: i, k, l, max_iter
        character(len=:), allocatable :: trim_type
        real :: u, v, w, rot_rates(3), G(6), R(6), R_1(6), R_2(6), J(6,6)
        real, dimension(:), allocatable :: dG
        real :: alpha, beta, da, de, dr, throttle
        real :: phi, theta, psi, var
        real :: error, gravity
        logical :: found
        logical :: solve_bank, solve_elev
        real :: gamma, gamma_1, gamma_2
        real :: theta_1, theta_2

        ! DR. HUNSAKER CODE
        ! PULL OUT SOLVER
        !   INCLUDE VERBOSE FLAG
        ! SOLVE BANK AND SOLVE ELEVATION LOGICAL
        ! SCT
        !   GET BANK
        ! SHSS
        !   GET BANK OR SIDESLIP
        ! ONLY ONE FUNCTION (CALC_R IS THE ONLY OTHER FUNCTION BECAUSE IT IS CALLED MULTIPLE TIMES)
        ! CALC R HAS SOLVE BANK LOGICAL

        call jsonx_get(j_main, "initial.trim.type", trim_type)
        call jsonx_get(j_main, "initial.trim.solver.finite_differenc_step_size", fd_step, default_value=0.01)
        call jsonx_get(j_main, "initial.trim.solver.relaxation_factor", relaxation, default_value=0.9)
        call jsonx_get(j_main, "initial.trim.solver.tolerance", tol)
        call jsonx_get(j_main, "initial.trim.solver.max_iterations", max_iter, default_value=100)

        write(*,*) "Trimming aircraft for ", trim_type
        !write(*,'(A,ES20.12)') "  --> Elevation angle set to theta [deg] = ", euler(2)*180./PI
        !write(*,'(A,ES20.12)') "  --> Bank angle set to phi [deg] = ", euler(1)*180./PI


        write(*,*) "Newton Solver Settings:"
        write(*,*) "  --> Finite Difference Step Size = ", fd_step
        write(*,*) "  --> Relaxation Factor = ", relaxation
        write(*,*) "  --> Tolerance = ", tol

        gravity = gravity_English(H)

        ! Initialize
        G = 0.0
        alpha = 0.0
        beta = 0.0
        gamma = 0.0
        phi = 0.0
        theta = 0.0
        psi = 0.0
        u = 0.0
        v = 0.0
        w = 0.0
        rot_rates = 0.0

        ! Check for sideslip angle
        call json_get(j_main, "initial.trim.bank_angle[deg]", phi, found)
        if (.not. found) then
            if (trim_type == "shss") then
                call json_get(j_main, "initial.trim.sideslip[deg]", beta, found)
                if (.not. found) then
                    write(*,*) "User must specify a bank or sideslip angle. Quitting..."
                    stop
                else
                    beta = beta*PI/180.
                    solve_bank = .true.
                end if
            else
                write(*,*) "User must specify a bank angle. Quitting..."
                stop
            end if
        else
            phi = phi*PI/180.
            solve_bank = .false.
        end if

        ! Check for elevation angle
        call json_get(j_main, "initial.trim.elevation_angle[deg]", theta, found)
        if (.not. found) then
            theta = 0.0
            call json_get(j_main, "initial.trim.climb_angle[deg]", gamma, found)
            if (.not. found) then
                solve_elev = .false.
                write(*,*) "User must specify a elevation or climb angle. Quitting..."
                stop
            else
                gamma = gamma*PI/180.
                solve_elev = .true.
            end if
        else
            theta = theta*PI/180.
            solve_elev = .false.
        end if

        write(*,*) "Initial theta [deg] = ", theta*180./PI
        write(*,*) "Initial gamma [deg] = ", gamma*180./PI
        write(*,*) "Initial phi[deg]    = ", phi*180./PI
        write(*,*) "Initial beta[deg]   = ", beta*180./PI

        do i = 1, max_iter

            ! Calculate velocities
            u = V_mag*cos(alpha)*cos(beta)
            v = V_mag*sin(beta)
            w = V_mag*sin(alpha)*cos(beta)

            if (solve_elev) then
                ! Calculate elevation angle
                write(*,*) "Calculating Elevation angle: "
                write(*,*) "phi: ", phi
                write(*,*) "gamma: ", gamma
                write(*,*) "v_mag: ", V_mag
                write(*,*) "u: ", u
                write(*,*) "v: ", v
                write(*,*) "w: ", w
                theta_1 = asin((u*V_mag*sin(gamma) + (v*sin(phi) + w*cos(phi))*sqrt(u**2 + (v*sin(phi) + w*cos(phi))**2 - &
                                                            V_mag**2 * sin(gamma)**2))/(u**2 + (v*sin(phi) + w*cos(phi))**2))
                theta_2 = asin((u*V_mag*sin(gamma) - (v*sin(phi) + w*cos(phi))*sqrt(u**2 + (v*sin(phi) + w*cos(phi))**2 - &
                                                            V_mag**2 * sin(gamma)**2))/(u**2 + (v*sin(phi) + w*cos(phi))**2))
                gamma_1 = asin((u*sin(theta_1) - (v*sin(phi) + w*cos(phi))*cos(theta_1))/V_mag)
                gamma_2 = asin((u*sin(theta_2) - (v*sin(phi) + w*cos(phi))*cos(theta_2))/V_mag)

                write(*,*) "    Theta 1: ", theta_1*180./PI
                write(*,*) "    Gamma 1: ", gamma_1*180./PI
                write(*,*) "    Theta 2: ", theta_2*180./PI
                write(*,*) "    Gamma 2: ", gamma_2*180./PI

                if (abs(gamma_1 - gamma) < 1.e-12) then
                    theta = theta_1
                else if (abs(gamma_2 - gamma) < 1.e-12) then
                    theta = theta_2
                else
                    write(*,*) "Trim solver could not find correct elevation angle. Quitting..."
                    stop
                end if
                write(*,*) "    Correct theta: ", theta
            end if

            ! Calculate rotation rates
            if (trim_type == "sct") then
                rot_rates(1) = -sin(theta)
                rot_rates(2) = sin(phi)*cos(theta)
                rot_rates(3) = cos(phi)*cos(theta)

                rot_rates = rot_rates*gravity*sin(phi)*cos(theta)/(u*cos(theta)*cos(phi) + w*sin(theta))

                write(*,*) "Updating rotation rates for steady coordinated turn:"
                write(*,'(A,ES20.12)') "  --> p [deg/s] = ", rot_rates(1)*180./PI
                write(*,'(A,ES20.12)') "  --> q [deg/s] = ", rot_rates(2)*180./PI
                write(*,'(A,ES20.12)') "  --> r [deg/s] = ", rot_rates(3)*180./PI
            end if


            write(*,*) "G defined as G = [alpha, beta, da, de, dr, throttle]"
            write(*,'(A,6ES20.12)') " G = ", G
            R = calc_R(V_mag, H, rot_rates, G, var, theta, psi, solve_bank)
            write(*,'(A,6ES20.12)') " R = ", R

            ! Solve for trim state
            !if (found_beta) then
                !call newtons_solver(V_mag, H, euler, rot_rates, G, fd_step, relaxation, error, beta)
            !else
                !call newtons_solver(V_mag, H, euler, rot_rates, G, fd_step, relaxation, error)
            !end if

            ! Set condition
            if (solve_bank) then
                G(2) = phi
                var = beta
            else
                G(2) = beta
                var = phi
            end if

            ! Assemble Jacobian
            do k = 1,6
                write(*,*) "Calculating gradient relative to G(", k, ")"
                G(k) = G(k) + fd_step
                write(*,*) "    Positive Finite Difference Step"
                write(*,'(A,6ES20.12)') "        G = ", G
                R_1 = calc_R(V_mag, H, rot_rates, G, var, theta, psi, solve_bank)
                write(*,'(A,6ES20.12)') "        R = ", R_1

                G(k) = G(k) - 2*fd_step
                write(*,*) "    Negative Finite Difference Step"
                write(*,'(A,6ES20.12)') "        G = ", G
                R_2 = calc_R(V_mag, H, rot_rates, G, var, theta, psi, solve_bank)
                write(*,'(A,6ES20.12)') "        R = ", R_2

                do l = 1,6
                    J(l,k) = (R_1(l) - R_2(l))/(2*fd_step)
                end do
                G(k) = G(k) + fd_step
            end do

            write(*,*) "Jacobian J = "
            do k = 1,6
                write(*,'(6ES20.12)') J(k,:)
            end do

            ! Calculate R
            R = calc_R(V_mag, H, rot_rates, G, var, theta, psi, solve_bank)

            call lu_solve(6, J, -R, dG)

            ! Update G
            G = G + relaxation*dG
            if (G(6) < 0.0) G(6) = 0.0

            write(*,*) 
            write(*,'(A,6ES20.12)') "Delta G: ", dG
            write(*,'(A,6ES20.12)') "New G:   ", G
        
            ! Calculate error
            error = maxval(abs(calc_R(V_mag, H, rot_rates, G, var, theta, psi, solve_bank)))


            write(*,'(A,I4,A,ES20.12)') "Iteration: ", i, " Error: ", error

            ! Update alpha, beta, and controls
            alpha = G(1)
            if (solve_bank) then
                phi = G(2)
            else
                beta = G(2)
            end if
            da = G(3)
            de = G(4)
            dr = G(5)
            throttle = G(6)
        
            if (error < tol) exit
        end do 

        write(*,'(A,ES20.12)') "Alpha (deg): ", alpha*180./PI
        write(*,'(A,ES20.12)') "Beta  (deg): ", beta*180./PI
        write(*,'(A,ES20.12)') "p     (deg/s): ", rot_rates(1)*180./PI
        write(*,'(A,ES20.12)') "q     (deg/s): ", rot_rates(2)*180./PI
        write(*,'(A,ES20.12)') "r     (deg/s): ", rot_rates(3)*180./PI
        write(*,'(A,ES20.12)') "Phi   (deg): ", phi*180./PI
        write(*,'(A,ES20.12)') "Theta (deg): ", theta*180./PI
        write(*,'(A,ES20.12)') "da    (deg): ", da*180./PI
        write(*,'(A,ES20.12)') "de    (deg): ", de*180./PI
        write(*,'(A,ES20.12)') "dr    (deg): ", dr*180./PI
        write(*,'(A,ES20.12)') "throttle    : ", throttle
        write(*,'(11A20)') "alpha[deg]", "beta[deg]", "p[deg/s]", "q[deg/s]", "r[deg/s]", "phi[deg]", &
                             "theta[deg]", "da[deg]", "de[deg]", "dr[deg]", "throttle"
        write(*,'(11ES20.12)') alpha*180./PI, beta*180./PI, rot_rates(1)*180./PI, rot_rates(2)*180./PI, &
                    rot_rates(3)*180./PI, phi*180./PI, theta*180./PI, da*180./PI, de*180./PI, dr*180./PI, throttle


        ! Set initial conditions
        controls(1) = da
        controls(2) = de
        controls(3) = dr
        controls(4) = throttle

        y_init = 0.0

        y_init(1) = V_mag*cos(alpha)*cos(beta)
        y_init(2) = V_mag*sin(beta)
        y_init(3) = V_mag*sin(alpha)*cos(beta)


        y_init(4) = rot_rates(1)
        y_init(5) = rot_rates(2)
        y_init(6) = rot_rates(3)

        y_init(9) = -H

        y_init(10:13) = euler_to_quat((/phi, theta, psi/))

    end subroutine init_trim

    subroutine init_state(V_mag, H)
        implicit none
        
        real, intent(inout) :: V_mag, H

        real :: alpha, beta, p, q, r, da, de, dr, throttle, phi, theta, psi
        
        call jsonx_get(j_main, "initial.state.alpha[deg]", alpha, default_value=0.0)
        call jsonx_get(j_main, "initial.state.beta[deg]", beta, default_value=0.0)
        call jsonx_get(j_main, "initial.state.p[deg/s]", p, default_value=0.0)
        call jsonx_get(j_main, "initial.state.q[deg/s]", q, default_value=0.0)
        call jsonx_get(j_main, "initial.state.r[deg/s]", r, default_value =0.0)
        call jsonx_get(j_main, "initial.state.aileron[deg]", da, default_value=0.0)
        call jsonx_get(j_main, "initial.state.elevator[deg]", de, default_value=0.0)
        call jsonx_get(j_main, "initial.state.rudder[deg]", dr, default_value=0.0)
        call jsonx_get(j_main, "initial.state.throttle", throttle, default_value=0.0)
        call jsonx_get(j_main, "initial.state.phi[deg]", phi, default_value=0.0)
        call jsonx_get(j_main, "initial.state.theta[deg]", theta, default_value=0.0)
        call jsonx_get(j_main, "initial.state.psi[deg]", psi, default_value=0.0)

        controls(1) = da*PI/180.
        controls(2) = de*PI/180.
        controls(3) = dr*PI/180.
        controls(4) = throttle

        phi = phi*PI/180.
        theta = theta*PI/180.
        psi = psi*PI/180.

        ! Set initial conditions
        alpha = alpha *PI/180.
        beta  = beta  *PI/180.

        y_init = 0.0

        y_init(1) = V_mag*cos(alpha)*cos(beta)
        y_init(2) = V_mag*sin(beta)
        y_init(3) = V_mag*sin(alpha)*cos(beta)


        y_init(4) = p*PI/180.
        y_init(5) = q*PI/180.
        y_init(6) = r*PI/180.

        y_init(9) = -H

        y_init(10:13) = euler_to_quat((/phi, theta, psi/))

    end subroutine init_state


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! FOR USE WITH CHAPTER 6 EXAMPLES !!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine init_6(filename)

        implicit none
        
        character(len=100), intent(in) :: filename
        type(json_value), pointer :: j_aircraft
        real :: V, H
        real, dimension(:), allocatable :: euler
        real :: alpha, beta, p, q, r
        real :: da, de, dr, throttle
        character(len=:), allocatable :: init_type

        ! Load JSON file
        call jsonx_load(filename, j_main)

        ! Initialize vehicle
        call jsonx_get(j_main, "vehicle", j_aircraft)
        call vehicle%init(j_aircraft)

        call jsonx_get(j_main, "simulation.verbose", verbose, default_value=.false.)
        call jsonx_get(j_main, "initial.airspeed[ft/s]", V)
        call jsonx_get(j_main, "initial.altitude[ft]", H)

        ! Get type of initialization
        call jsonx_get(j_main, "initial.type", init_type)

        select case(init_type)
        case("state")
            call init_state(V,H)
        case("trim")
            call init_trim(V,H)
        case default
            write(*,*) "!!! Type "//init_type//" is not recognized as a valid init type. Quitting..."
            stop
        end select
        
    end subroutine init_6

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! FOR USE WITH CHAPTER 5 EXAMPLES !!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine init_5(filename)

        implicit none
        
        character(len=100), intent(in) :: filename
        type(json_value), pointer :: j_aircraft
        real :: V, H, theta, phi, psi
        real :: alpha, beta, p, q, r
        real :: da, de, dr, throttle

        ! Load JSON file
        call jsonx_load(filename, j_main)

        ! Initialize vehicle
        call jsonx_get(j_main, "vehicle", j_aircraft)
        call vehicle%init(j_aircraft)

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

        controls(1) = da*PI/180.
        controls(2) = de*PI/180.
        controls(3) = dr*PI/180.
        controls(4) = throttle

        ! Set initial conditions
        phi   = phi   *PI/180.
        theta = theta *PI/180.
        psi   = psi   *PI/180.
        
        alpha = alpha *PI/180.
        beta  = beta  *PI/180.

        y_init = 0.0

        y_init = V*cos(alpha)*cos(beta)
        y_init(2) = V*sin(beta)
        y_init(3) = V*sin(alpha)*cos(beta)


        y_init(4) = p*PI/180.
        y_init(5) = q*PI/180.
        y_init(6) = r*PI/180.

        y_init(9) = -H

        y_init(10:13) = euler_to_quat([phi, theta, psi])

    end subroutine init_5

    subroutine run()

        implicit none

        real :: t, y(13), y_temp(13)
        integer :: io_unit
        real :: dt, tf
        logical :: t_R
        real :: t_p, t_c, start_time, end_time


        call jsonx_get(j_main, "simulation.timestep[s]", dt, default_value=0.01)
        call jsonx_get(j_main, "simulation.total_time[s]", tf)

        t = 0.0
        y = y_init

        open(newunit=io_unit, file='sim_output.csv', status='replace', action='write')
        write(io_unit,*) 'time[s], u[ft/s], v[ft/s], w[ft/s], p[rad/s], q[rad/s], r[rad/s], xf[ft], yf[ft], zf[ft], e0, ex, ey, ez'
        write(io_unit,'(ES20.12,A,ES20.12,A,ES20.12,A,ES20.12,A,ES20.12,A,ES20.12,A,ES20.12,A,ES20.12,A,ES20.12, &
                        A,ES20.12,A,ES20.12,A,ES20.12,A,ES20.12,A,ES20.12)') &
                            t,',',y(1),',',y(2),',',y(3),',',y(4),',',y(5),',' &
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
            t_p = get_time()
            y = runge_kutta(t, y, 0.01)
            t_c = get_time()
            dt = t_c - t_p
            y = y_temp
            t_p = t_c
        end if



        start_time = get_time()
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
                t_c = get_time()
                dt = t_c - t_p
                t_p = t_c
            end if
            if (verbose) then
                write(*,*) "----------------------"
                write(*,*) "----------------------"
            end if

            call quat_norm(y(10:13))

            write(io_unit,'(ES20.12,A,ES20.12,A,ES20.12,A,ES20.12,A,ES20.12,A,ES20.12,A,ES20.12,A,&
                            ES20.12,A,ES20.12,A,ES20.12,A,ES20.12,A,ES20.12,A,ES20.12,A,ES20.12)') &
                            t,',',y(1),',',y(2),',',y(3),',',y(4),',',y(5),',' &
                            ,y(6),',',y(7),',',y(8),',',y(9),',',y(10),',',y(11),',',y(12),',',y(13)

        end do
        end_time = get_time()

        close(io_unit)

        write(*,*) "Simulation Finished in ", end_time - start_time, " seconds."
        
        if (t_R) then
            write(*,*) "Real time simulation error: ", tf - (end_time - start_time), " seconds."
        end if

    end subroutine run 

end module simulation_m 