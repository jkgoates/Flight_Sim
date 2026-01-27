module vehicle_m 
    use goates_m
    use linalg_mod
    use jsonx_m

    implicit none

    type stall_settings_t
        real :: alpha_0, alpha_s, lambda_b, minval        
    end type stall_settings_t

    type trim_settings_t
        real :: delta, gamma, tol
        character(:), allocatable :: type
    end type trim_settings_t
    
    type :: aircraft
        
        character(:), allocatable :: type
        logical :: run_physics
        real :: M, Ixx, Iyy, Izz, Ixy, Iyz, Ixz ! Mass and Inertia matrix
        real, dimension(:), allocatable :: CG_shift ! CG shift from reference point (ft)
        real, allocatable :: h(:) ! Gyroscopic components
        real :: T0, a ! Thrust parameters
        real, dimension(:), allocatable :: t_location, t_orientation
        real :: S_w, b, c ! Reference area, span, chord
        real :: CL_0, CL_alpha, CL_alphahat, CL_qbar, CL_de ! Lift coefficients
        real :: CS_beta, CS_pbar, CS_alpha_pbar, CS_rbar, CS_da, CS_dr ! Side force coefficients
        real :: CD_L0, CD_CL1, CD_CL1_CL1, CD_CS_CS, CD_qbar, CD_alpha_qbar, CD_de, CD_alpha_de, CD_de_de ! Drag coefficients
        real :: Cell_0, Cell_1, Cell_beta, Cell_pbar, Cell_alpha_rbar, Cell_rbar, Cell_da, Cell_dr ! Rolling moment coefficients
        real :: Cm_0, Cm_alpha, Cm_qbar, Cm_alphahat, Cm_de ! Pitching moment coefficients
        real :: Cn_beta, Cn_pbar, Cn_alpha_pbar, Cn_rbar, Cn_da, Cn_alpha_da, Cn_dr ! Yawing moment coefficients
        ! real :: rho0

        real, allocatable :: ks(:), kd(:) ! Landing gear spring and damping coefficients
        real, allocatable :: g_b(:,:) ! Landing gear body fixed coordinates
        real, allocatable :: th_b(:) ! tail hook location
        real :: ag_ks, ag_kd, ag_length, ag_distance ! arresting gear coefficients
        logical :: ag_engaged
        real, allocatable :: collision_points(:,:) ! Aircraft collision points
    
        real :: init_V, init_alt, init_states
        real, allocatable :: init_eul(:)

        type(stall_settings_t) :: CD_stall, CL_stall, Cm_stall
        type(trim_settings_t) :: trim
        logical :: include_stall

        real :: states(13), controls(4)

    contains

        procedure :: mass_inertia   => aircraft_mass_inertia
        procedure :: pseudo_aero    => aircraft_pseudo_aero
        procedure :: gyroscopic     => aircraft_gyroscopic
        procedure :: thrust         => aircraft_thrust
        procedure :: init           => aircraft_init
        procedure :: init_to_state  => aircraft_init_to_state
        procedure :: init_to_trim   => aircraft_init_to_trim
        !procedure :: landing_gear   => aircraft_landing_gear
        !procedure :: check_collision => aircraft_check_collision
        !procedure :: arresting_gear => aircraft_arresting_gear
        procedure :: print_aero_table => aircraft_print_aero_table
        procedure :: tick_states => aircraft_tick_states
        procedure :: calc_R        => aircraft_calc_R
        procedure :: runge_kutta => aircraft_runge_kutta
        procedure :: diff_eq       => aircraft_diff_eq
        procedure :: newtons_method => aircraft_newtons_method
    
    end type aircraft
    
contains

    subroutine aircraft_init(this, settings)
    
        class(aircraft), intent(inout) :: this
        type(json_value), pointer, intent(in) :: settings

        type(json_value), pointer :: reference, coefficients, aerodynamics, p1, p2, j_initial
        real, allocatable :: dummy_loc(:)
        integer :: N, cnt, i
        logical :: found

        real :: V, H, lat, long
        real, allocatable :: Euler(:)
        character(len=:), allocatable :: init_type

        ! Parse JSON settings
        call jsonx_get(settings, "initial", j_initial)
        call jsonx_get(settings, "aerodynamics", aerodynamics)
        call jsonx_get(aerodynamics, "reference", reference)

        call jsonx_get(settings, "run_physics", this%run_physics)
        call jsonx_get(settings, "type", this%type)

        if (this%type == "aircraft" .or. this%type == "arrow") then
            call jsonx_get(aerodynamics, "coefficients", coefficients)
        end if
        ! CG Shift
        !call jsonx_get(settings, "CG_shift[ft]", this%CG_shift, 0.0, 3)

        ! Mass properties
        call jsonx_get(settings, "mass.weight[lbf]", this%M)
        call jsonx_get(settings, "mass.Ixx[slug-ft^2]", this%Ixx, 0.0)
        call jsonx_get(settings, "mass.Iyy[slug-ft^2]", this%Iyy, 0.0)
        call jsonx_get(settings, "mass.Izz[slug-ft^2]", this%Izz, 0.0)
        call jsonx_get(settings, "mass.Ixy[slug-ft^2]", this%Ixy, 0.0)
        call jsonx_get(settings, "mass.Iyz[slug-ft^2]", this%Iyz, 0.0)
        call jsonx_get(settings, "mass.Ixz[slug-ft^2]", this%Ixz, 0.0)
        call jsonx_get(settings, "mass.h[slug-ft^2/s]", this%h, 0.0, 3)

        ! Thrust Properties
        call jsonx_get(settings, "thrust.T0[lbf]", this%T0, 0.0)
        call jsonx_get(settings, "thrust.Ta", this%a, 0.0)
        call jsonx_get(settings, "thrust.location[ft]", this%t_location, 0.0, 3)
        call jsonx_get(settings, "thrust.orientation[deg]", this%t_orientation, 0.0, 3)
        this%t_orientation = this%t_orientation*pi/180.0

        ! Reference properties
        call jsonx_get(reference, "area[ft^2]", this%S_w)
        call jsonx_get(reference, "longitudinal_length[ft]", this%c)
        call jsonx_get(reference, "lateral_length[ft]", this%b)
        call jsonx_get(reference, "location[ft]", this%CG_shift, 0.0, 3)

        if (this%type == 'aircraft') then
            ! Lift coefficients
            call jsonx_get(coefficients, "CL.0", this%CL_0)
            call jsonx_get(coefficients, "CL.alpha", this%CL_alpha)
            call jsonx_get(coefficients, "CL.alphahat", this%CL_alphahat)
            call jsonx_get(coefficients, "CL.qbar", this%CL_qbar)
            call jsonx_get(coefficients, "CL.elevator", this%CL_de)

            ! Side force coefficients
            call jsonx_get(coefficients, "CS.beta", this%CS_beta)
            call jsonx_get(coefficients, "CS.pbar", this%CS_pbar)
            call jsonx_get(coefficients, "CS.alpha_pbar", this%CS_alpha_pbar)
            call jsonx_get(coefficients, "CS.rbar", this%CS_rbar)
            call jsonx_get(coefficients, "CS.aileron", this%CS_da)
            call jsonx_get(coefficients, "CS.rudder", this%CS_dr)

            ! Drag Coefficients
            call jsonx_get(coefficients, "CD.L0", this%CD_L0)
            call jsonx_get(coefficients, "CD.CL1", this%CD_CL1)
            call jsonx_get(coefficients, "CD.CL1_CL1", this%CD_CL1_CL1)
            call jsonx_get(coefficients, "CD.CS_CS", this%CD_CS_CS)
            call jsonx_get(coefficients, "CD.qbar", this%CD_qbar)
            call jsonx_get(coefficients, "CD.alpha_qbar", this%CD_alpha_qbar)
            call jsonx_get(coefficients, "CD.elevator", this%CD_de)
            call jsonx_get(coefficients, "CD.alpha_elevator", this%CD_alpha_de)
            call jsonx_get(coefficients, "CD.elevator_elevator", this%CD_de_de)

            ! Rolling moment coefficients
            call jsonx_get(coefficients, "Cl.beta", this%Cell_beta)
            call jsonx_get(coefficients, "Cl.pbar", this%Cell_pbar)
            call jsonx_get(coefficients, "Cl.alpha_rbar", this%Cell_alpha_rbar)
            call jsonx_get(coefficients, "Cl.rbar", this%Cell_rbar)
            call jsonx_get(coefficients, "Cl.aileron", this%Cell_da)
            call jsonx_get(coefficients, "Cl.rudder", this%Cell_dr)

            ! Pitching moment coefficients
            call jsonx_get(coefficients, "Cm.0", this%Cm_0)
            call jsonx_get(coefficients, "Cm.alpha", this%Cm_alpha)
            call jsonx_get(coefficients, "Cm.qbar", this%Cm_qbar)
            call jsonx_get(coefficients, "Cm.alphahat", this%Cm_alphahat)
            call jsonx_get(coefficients, "Cm.elevator", this%Cm_de)

            ! Yawing moment coefficients
            call jsonx_get(coefficients, "Cn.beta", this%Cn_beta)
            call jsonx_get(coefficients, "Cn.pbar", this%Cn_pbar)
            call jsonx_get(coefficients, "Cn.alpha_pbar", this%Cn_alpha_pbar)
            call jsonx_get(coefficients, "Cn.rbar", this%Cn_rbar)
            call jsonx_get(coefficients, "Cn.aileron", this%Cn_da)
            call jsonx_get(coefficients, "Cn.alpha_aileron", this%Cn_alpha_da)
            call jsonx_get(coefficients, "Cn.rudder", this%Cn_dr)
        
        else if (this%type == 'arrow') then

            call jsonx_get(coefficients, "CL.alpha", this%CL_alpha)
            call jsonx_get(coefficients, "CD.L0", this%CD_L0)
            call jsonx_get(coefficients, "CD.CL1_CL1", this%CD_CL1_CL1)
            call jsonx_get(coefficients, "Cl.0", this%Cell_0)
            call jsonx_get(coefficients, "Cl.pbar", this%Cell_pbar)
            call jsonx_get(coefficients, "Cm.alpha", this%Cm_alpha)
            call jsonx_get(coefficients, "Cm.qbar", this%Cm_qbar)
        end if


        ! Stall
        call jsonx_get(aerodynamics, "stall.include_stall", this%include_stall, .false.)
        if (this%include_stall) then
            call jsonx_get(aerodynamics, "stall.CL.alpha_0[deg]", this%CL_stall%alpha_0)
            call jsonx_get(aerodynamics, "stall.CL.alpha_s[deg]", this%CL_stall%alpha_s)
            call jsonx_get(aerodynamics, "stall.CL.lambda_b", this%CL_stall%lambda_b)
            this%CL_stall%alpha_0 = this%CL_stall%alpha_0*pi/180.
            this%CL_stall%alpha_s = this%CL_stall%alpha_s*pi/180.

            call jsonx_get(aerodynamics, "stall.CD.alpha_0[deg]", this%CD_stall%alpha_0)
            call jsonx_get(aerodynamics, "stall.CD.alpha_s[deg]", this%CD_stall%alpha_s)
            call jsonx_get(aerodynamics, "stall.CD.lambda_b", this%CD_stall%lambda_b)
            this%CD_stall%alpha_0 = this%CD_stall%alpha_0*pi/180.
            this%CD_stall%alpha_s = this%CD_stall%alpha_s*pi/180.

            call jsonx_get(aerodynamics, "stall.Cm.alpha_0[deg]", this%Cm_stall%alpha_0)
            call jsonx_get(aerodynamics, "stall.Cm.alpha_s[deg]", this%Cm_stall%alpha_s)
            call jsonx_get(aerodynamics, "stall.Cm.lambda_b", this%Cm_stall%lambda_b)
            call jsonx_get(aerodynamics, "stall.Cm.min", this%Cm_stall%minval)
            this%Cm_stall%alpha_0 = this%Cm_stall%alpha_0*pi/180.
            this%Cm_stall%alpha_s = this%Cm_stall%alpha_s*pi/180.
        end if

        !! Landing Gear
        !call jsonx_get(settings, "landing_gear", p1)
        !N = json_value_count(p1)
        !cnt= 0
        !do i = 1, N
            !call json_value_get(p1, i, p2)
            !if (p2%name(1:1) == 'x') cycle
            !cnt = cnt+1
        !end do
        !allocate(this%ks(cnt))
        !allocate(this%kd(cnt))
        !allocate(this%g_b(cnt, 3))
        !cnt = 1
        !do i = 1, N
            !call json_value_get(p1, i, p2)
            !if (p2%name(1:1) == 'x') cycle
            !call jsonx_get(p2, "location[ft]", dummy_loc, 0.0, 3)
            !this%g_b(cnt, :) = dummy_loc
            !call jsonx_get(p2, "spring_constant[lb/ft]", this%ks(cnt))
            !call jsonx_get(p2, "damping_constant[lb-s/ft]", this%kd(cnt))
            !cnt = cnt+1
        !end do

        !! Collision points
        !call jsonx_get(settings, "collision_points", p1)
        !N = json_value_count(p1)
        !cnt = 0
        !do i = 1, N
            !call json_value_get(p1, i, p2)
            !if (p2%name(1:1) == 'x') cycle
            !cnt = cnt+1
        !end do
        !allocate(this%collision_points(cnt,3))
        !cnt = 1
        !do i = 1, N
            !call json_value_get(p1, i, p2)
            !if (p2%name(1:1) == 'x') cycle
            !call jsonx_get(p1, p2%name, dummy_loc, 0.0, 3)
            !this%collision_points(cnt, :) = dummy_loc
            !cnt = cnt+1
        !end do

        !! Arrestor Gear
        !call jsonx_get(settings, "arresting_gear", p1)
        !call jsonx_get(p1, "spring_constant[lb/ft]", this%ag_ks)
        !call jsonx_get(p1, "damping_constant[lb-s/ft]", this%ag_kd)
        !call jsonx_get(p1, "cable_length[ft]", this%ag_length)
        !call jsonx_get(p1, "cable_distance[ft]", this%ag_distance)
        !call jsonx_get(p1, "tail_hook[ft]", this%th_b)
        !this%ag_engaged = .false.

        ! State
        call jsonx_get(j_initial, "airspeed[ft/s]", this%init_V)
        call jsonx_get(j_initial, "altitude[ft]", this%init_alt)
        call jsonx_get(j_initial, "latitude[deg]", lat)
        call jsonx_get(j_initial, "longitude[deg]", long)
        call jsonx_get(j_initial, "Euler_angles[deg]", this%init_eul, 3)

        ! Get type of initialization
        call jsonx_get(j_initial, "type", init_type)


        select case(init_type)
        case("state")
            call this%init_to_state(j_initial)
        case("trim")
            call this%init_to_trim(j_initial, V,H, Euler)
        case default
            write(*,*) "!!! Type "//init_type//" is not recognized as a valid init type. Quitting..."
            stop
        end select

    end subroutine aircraft_init

    subroutine aircraft_init_to_state(this, j_initial)
        implicit none
        
        class(aircraft), intent(inout) :: this
        type(json_value), pointer, intent(in) :: j_initial

        real :: alpha, beta, p, q, r, da, de, dr, throttle, phi, theta, psi, xf, yf, zf
        
        call jsonx_get(j_initial, "state.angle_of_attack[deg]", alpha, default_value=0.0)
        call jsonx_get(j_initial, "state.sideslip_angle[deg]", beta, default_value=0.0)
        call jsonx_get(j_initial, "state.p[deg/s]", p, default_value=0.0)
        call jsonx_get(j_initial, "state.q[deg/s]", q, default_value=0.0)
        call jsonx_get(j_initial, "state.r[deg/s]", r, default_value =0.0)
        call jsonx_get(j_initial, "state.aileron[deg]", da, default_value=0.0)
        call jsonx_get(j_initial, "state.elevator[deg]", de, default_value=0.0)
        call jsonx_get(j_initial, "state.rudder[deg]", dr, default_value=0.0)
        call jsonx_get(j_initial, "state.throttle", throttle, default_value=0.0)

        this%controls(1) = da*PI/180.
        this%controls(2) = de*PI/180.
        this%controls(3) = dr*PI/180.
        this%controls(4) = throttle

        this%init_eul = this%init_eul*PI/180.
        !phi = phi*PI/180.
        !theta = theta*PI/180.
        !psi = psi*PI/180.

        ! Set initial conditions
        alpha = alpha *PI/180.
        beta  = beta  *PI/180.

        this%states = 0.0

        this%states(1) = this%init_V*cos(alpha)*cos(beta)
        this%states(2) = this%init_V*sin(beta)
        this%states(3) = this%init_V*sin(alpha)*cos(beta)


        this%states(4) = p*PI/180.
        this%states(5) = q*PI/180.
        this%states(6) = r*PI/180.

        this%states(7) = xf
        this%states(8) = yf
        this%states(9) = -this%init_alt

        this%states(10:13) = euler_to_quat(this%init_eul)

        write(*,*) "INITIAL STATES", this%states

    end subroutine aircraft_init_to_state

    subroutine aircraft_init_to_trim(this, j_initial)

        implicit none

        class(aircraft), intent(inout) :: this
        type(json_value), pointer, intent(in) :: j_initial


        call jsonx_get(j_initial, "trim.type", this%trim%type)


    end subroutine aircraft_init_to_trim


    !subroutine aircraft_init_to_trim_old(this, j_initial, V_mag, H)

        !implicit none
        
        !class(aircraft), intent(inout) :: this
        !type(json_value), pointer, intent(in) :: j_initial
        !real, intent(in) :: V_mag, H
        !real :: y(13)

        !real :: fd_step, relaxation, tol
        !integer :: i, k, l, max_iter
        !character(len=:), allocatable :: trim_type
        !real :: u, v, w, rot_rates(3), G(6), R(6), R_1(6), R_2(6), J(6,6)
        !real, dimension(:), allocatable :: dG
        !real :: alpha, beta, da, de, dr, throttle
        !real :: phi, theta, psi, var
        !real :: error, gravity
        !logical :: found
        !logical :: solve_bank, solve_elev
        !real :: gamma, gamma_1, gamma_2
        !real :: theta_1, theta_2

        !call jsonx_get(j_initial, "trim.type", trim_type)
        !call jsonx_get(j_initial, "trim.solver.finite_differenc_step_size", fd_step, default_value=0.01)
        !call jsonx_get(j_initial, "trim.solver.relaxation_factor", relaxation, default_value=0.9)
        !call jsonx_get(j_initial, "trim.solver.tolerance", tol)
        !call jsonx_get(j_initial, "trim.solver.max_iterations", max_iter, default_value=100)

        !write(*,*) "Trimming aircraft for ", trim_type

        !write(*,*) "Newton Solver Settings:"
        !write(*,*) "  --> Finite Difference Step Size = ", fd_step
        !write(*,*) "  --> Relaxation Factor = ", relaxation
        !write(*,*) "  --> Tolerance = ", tol

        !gravity = gravity_English(H)

        !! Initialize
        !G = 0.0
        !alpha = 0.0
        !beta = 0.0
        !gamma = 0.0
        !phi = 0.0
        !theta = 0.0
        !psi = 0.0
        !u = 0.0
        !v = 0.0
        !w = 0.0
        !rot_rates = 0.0

        !! Check for sideslip angle
        !call json_get(j_initial, "trim.bank_angle[deg]", phi, found)
        !if (.not. found) then
            !if (trim_type == "shss") then
                !call json_get(j_initial, "trim.sideslip[deg]", beta, found)
                !if (.not. found) then
                    !write(*,*) "User must specify a bank or sideslip angle. Quitting..."
                    !stop
                !else
                    !beta = beta*PI/180.
                    !solve_bank = .true.
                !end if
            !else
                !write(*,*) "User must specify a bank angle. Quitting..."
                !stop
            !end if
        !else
            !phi = phi*PI/180.
            !solve_bank = .false.
        !end if

        !! Check for elevation angle
        !call json_get(j_initial, "trim.elevation_angle[deg]", theta, found)
        !if (.not. found) then
            !theta = 0.0
            !call json_get(j_initial, "trim.climb_angle[deg]", gamma, found)
            !if (.not. found) then
                !solve_elev = .false.
                !write(*,*) "User must specify a elevation or climb angle. Quitting..."
                !stop
            !else
                !gamma = gamma*PI/180.
                !solve_elev = .true.
            !end if
        !else
            !theta = theta*PI/180.
            !solve_elev = .false.
        !end if

        !write(*,*) "Initial theta [deg] = ", theta*180./PI
        !write(*,*) "Initial gamma [deg] = ", gamma*180./PI
        !write(*,*) "Initial phi[deg]    = ", phi*180./PI
        !write(*,*) "Initial beta[deg]   = ", beta*180./PI

        !do i = 1, max_iter

            !! Calculate velocities
            !u = V_mag*cos(alpha)*cos(beta)
            !v = V_mag*sin(beta)
            !w = V_mag*sin(alpha)*cos(beta)

            !if (solve_elev) then
                !! Calculate elevation angle
                !write(*,*) "Calculating Elevation angle: "
                !write(*,*) "phi: ", phi
                !write(*,*) "gamma: ", gamma
                !write(*,*) "v_mag: ", V_mag
                !write(*,*) "u: ", u
                !write(*,*) "v: ", v
                !write(*,*) "w: ", w
                !theta_1 = asin((u*V_mag*sin(gamma) + (v*sin(phi) + w*cos(phi))*sqrt(u**2 + (v*sin(phi) + w*cos(phi))**2 - &
                                                            !V_mag**2 * sin(gamma)**2))/(u**2 + (v*sin(phi) + w*cos(phi))**2))
                !theta_2 = asin((u*V_mag*sin(gamma) - (v*sin(phi) + w*cos(phi))*sqrt(u**2 + (v*sin(phi) + w*cos(phi))**2 - &
                                                            !V_mag**2 * sin(gamma)**2))/(u**2 + (v*sin(phi) + w*cos(phi))**2))
                !gamma_1 = asin((u*sin(theta_1) - (v*sin(phi) + w*cos(phi))*cos(theta_1))/V_mag)
                !gamma_2 = asin((u*sin(theta_2) - (v*sin(phi) + w*cos(phi))*cos(theta_2))/V_mag)

                !write(*,*) "    Theta 1: ", theta_1*180./PI
                !write(*,*) "    Gamma 1: ", gamma_1*180./PI
                !write(*,*) "    Theta 2: ", theta_2*180./PI
                !write(*,*) "    Gamma 2: ", gamma_2*180./PI

                !if (abs(gamma_1 - gamma) < 1.e-12) then
                    !theta = theta_1
                !else if (abs(gamma_2 - gamma) < 1.e-12) then
                    !theta = theta_2
                !else
                    !write(*,*) "Trim solver could not find correct elevation angle. Quitting..."
                    !stop
                !end if
                !write(*,*) "    Correct theta: ", theta
            !end if

            !! Calculate rotation rates
            !if (trim_type == "sct") then
                !rot_rates(1) = -sin(theta)
                !rot_rates(2) = sin(phi)*cos(theta)
                !rot_rates(3) = cos(phi)*cos(theta)

                !rot_rates = rot_rates*gravity*sin(phi)*cos(theta)/(u*cos(theta)*cos(phi) + w*sin(theta))

                !write(*,*) "Updating rotation rates for steady coordinated turn:"
                !write(*,'(A,ES20.12)') "  --> p [deg/s] = ", rot_rates(1)*180./PI
                !write(*,'(A,ES20.12)') "  --> q [deg/s] = ", rot_rates(2)*180./PI
                !write(*,'(A,ES20.12)') "  --> r [deg/s] = ", rot_rates(3)*180./PI
            !end if


            !write(*,*) "G defined as G = [alpha, beta, da, de, dr, throttle]"
            !write(*,'(A,6ES20.12)') " G = ", G
            !R = this%calc_R(V_mag, H, rot_rates, G, var, theta, psi, solve_bank)
            !write(*,'(A,6ES20.12)') " R = ", R

            !! Solve for trim state
            !!if (found_beta) then
                !!call newtons_solver(V_mag, H, euler, rot_rates, G, fd_step, relaxation, error, beta)
            !!else
                !!call newtons_solver(V_mag, H, euler, rot_rates, G, fd_step, relaxation, error)
            !!end if

            !! Set condition
            !if (solve_bank) then
                !G(2) = phi
                !var = beta
            !else
                !G(2) = beta
                !var = phi
            !end if

            !! Assemble Jacobian
            !do k = 1,6
                !write(*,*) "Calculating gradient relative to G(", k, ")"
                !G(k) = G(k) + fd_step
                !write(*,*) "    Positive Finite Difference Step"
                !write(*,'(A,6ES20.12)') "        G = ", G
                !R_1 = this%calc_R(V_mag, H, rot_rates, G, var, theta, psi, solve_bank)
                !write(*,'(A,6ES20.12)') "        R = ", R_1

                !G(k) = G(k) - 2*fd_step
                !write(*,*) "    Negative Finite Difference Step"
                !write(*,'(A,6ES20.12)') "        G = ", G
                !R_2 = this%calc_R(V_mag, H, rot_rates, G, var, theta, psi, solve_bank)
                !write(*,'(A,6ES20.12)') "        R = ", R_2

                !do l = 1,6
                    !J(l,k) = (R_1(l) - R_2(l))/(2*fd_step)
                !end do
                !G(k) = G(k) + fd_step
            !end do

            !write(*,*) "Jacobian J = "
            !do k = 1,6
                !write(*,'(6ES20.12)') J(k,:)
            !end do

            !! Calculate R
            !R = this%calc_R(V_mag, H, rot_rates, G, var, theta, psi, solve_bank)

            !call lu_solve(6, J, -R, dG)

            !! Update G
            !G = G + relaxation*dG
            !if (G(6) < 0.0) G(6) = 0.0

            !write(*,*) 
            !write(*,'(A,6ES20.12)') "Delta G: ", dG
            !write(*,'(A,6ES20.12)') "New G:   ", G
        
            !! Calculate error
            !error = maxval(abs(this%calc_R(V_mag, H, rot_rates, G, var, theta, psi, solve_bank)))


            !write(*,'(A,I4,A,ES20.12)') "Iteration: ", i, " Error: ", error

            !! Update alpha, beta, and controls
            !alpha = G(1)
            !if (solve_bank) then
                !phi = G(2)
            !else
                !beta = G(2)
            !end if
            !da = G(3)
            !de = G(4)
            !dr = G(5)
            !throttle = G(6)
        
            !if (error < tol) exit
        !end do 

        !write(*,'(A,ES20.12)') "Alpha (deg): ", alpha*180./PI
        !write(*,'(A,ES20.12)') "Beta  (deg): ", beta*180./PI
        !write(*,'(A,ES20.12)') "p     (deg/s): ", rot_rates(1)*180./PI
        !write(*,'(A,ES20.12)') "q     (deg/s): ", rot_rates(2)*180./PI
        !write(*,'(A,ES20.12)') "r     (deg/s): ", rot_rates(3)*180./PI
        !write(*,'(A,ES20.12)') "Phi   (deg): ", phi*180./PI
        !write(*,'(A,ES20.12)') "Theta (deg): ", theta*180./PI
        !write(*,'(A,ES20.12)') "da    (deg): ", da*180./PI
        !write(*,'(A,ES20.12)') "de    (deg): ", de*180./PI
        !write(*,'(A,ES20.12)') "dr    (deg): ", dr*180./PI
        !write(*,'(A,ES20.12)') "throttle    : ", throttle
        !write(*,'(11A20)') "alpha[deg]", "beta[deg]", "p[deg/s]", "q[deg/s]", "r[deg/s]", "phi[deg]", &
                             !"theta[deg]", "da[deg]", "de[deg]", "dr[deg]", "throttle"
        !write(*,'(11ES20.12)') alpha*180./PI, beta*180./PI, rot_rates(1)*180./PI, rot_rates(2)*180./PI, &
                    !rot_rates(3)*180./PI, phi*180./PI, theta*180./PI, da*180./PI, de*180./PI, dr*180./PI, throttle


        !! Set initial conditions
        !this%controls(1) = da
        !this%controls(2) = de
        !this%controls(3) = dr
        !this%controls(4) = throttle

        !this%states = 0.0

        !this%states(1) = V_mag*cos(alpha)*cos(beta)
        !this%states(2) = V_mag*sin(beta)
        !this%states(3) = V_mag*sin(alpha)*cos(beta)


        !this%states(4) = rot_rates(1)
        !this%states(5) = rot_rates(2)
        !this%states(6) = rot_rates(3)

        !this%states(9) = -H

        !this%states(10:13) = euler_to_quat((/phi, theta, psi/))
    !end subroutine aircraft_init_to_trim_old

    function aircraft_newtons_method (N, x, delta, gamma, tol) result(temp_x)

        implicit none
        
        integer, intent(in) :: N
        real, intent(in) :: x(N), delta, gamma, tol

        real :: R(N), R1(N), R2(N), J(N,N), temp_x(N), err
        real, allocatable :: dx(:)
        integer :: i, k

        temp_x = x

        R = calc_R_demo(N,temp_x)
        err = maxval(abs(R))

        do while (err > tol)

            do k = 1, N
                temp_x(k) = temp_x(k) + delta
                R1 = calc_R_demo(N,temp_x)
                temp_x(k) = temp_x(k) - 2*delta
                R2 = calc_R_demo(N,temp_x)
                do i = 1, N
                    J(i,k) = (R1(i) - R2(i))/(2*delta)
                end do
                temp_x(k) = temp_x(k) + delta
            end do
            
            J = -J

            call lu_solve(N, J, R, dx)
            
            temp_x = temp_x + gamma*dx

            R = calc_R_demo(N,temp_x)
            err = maxval(abs(R))
            write(*,*) "temp_x: ", temp_x, err
        end do

    end function aircraft_newtons_method

    function aircraft_calc_R(this, N, x, V, loc, euler) result(R)

        implicit none
        
        class(aircraft), intent(inout) :: this
        integer, intent(in) :: N
        real, intent(in) :: x(N), V, loc(3), euler(3)



        real :: y_temp(13), alpha, beta, dy_dt(13), R(6), gravity

        gravity = gravity_English(-loc(3))

        alpha = x(1)
        beta = x(2)

        y_temp(1) = V*cos(alpha)*cos(beta)
        y_temp(2) = V*sin(beta)
        y_temp(3) = V*sin(alpha)*cos(beta)

        ! sct
        y_temp(4) = -sin(euler(2))
        y_temp(5) = sin(euler(1))*cos(euler(2))
        y_temp(6) = cos(euler(1))*cos(euler(2))
        y_temp(4:6) = y_temp(4:6)*gravity*sin(euler(1))*cos(euler(2))&
                        /(y_temp(1)*cos(euler(2))*cos(euler(1)) + y_temp(3)*sin(euler(2)))
        
        y_temp(7:9) = loc

        y_temp(10:13) = euler_to_quat(euler)

        this%controls = x(3:6)

        dy_dt = this%diff_eq(y_temp)


    end function aircraft_calc_R

    function calc_R_demo(N, x) result(R)

        implicit none
        
        integer, intent(in) :: N
        real, intent(in) :: x(N)

        real :: R(N)

        R(1) = x(2) - x(1)**2 + 5
        R(2) = x(2) - x(1)**2 - (x(1) - 1)

    end function calc_R_demo

    function aircraft_calc_R(this, V, H, rot_rates, G, var, theta, psi, solve_bank) result(R)

        implicit none
        class(aircraft), intent(inout) :: this
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
        this%controls = G(3:6)

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
        dy_dt = this%diff_eq( y_temp)
        
        R = dy_dt(1:6)

    end function aircraft_calc_R

    function aircraft_runge_kutta(this, y_0, dt) result(y)

        implicit none

        class(aircraft), intent(inout) :: this
        real, intent(in) :: y_0(13), dt
        real :: y(13)

        real :: k1(13), k2(13), k3(13), k4(13)


        k1 = this%diff_eq(y_0)
        k2 = this%diff_eq(y_0 + k1*0.5*dt)
        k3 = this%diff_eq(y_0 + k2*0.5*dt)
        k4 = this%diff_eq(y_0 + k3*dt)

        y = y_0 + (dt*one_sixth)*(k1 + 2*k2 + 2*k3 + k4)

    end function aircraft_runge_kutta

    function aircraft_diff_eq(this, y) result(dy_dt)

        implicit none
        
        class(aircraft), intent(inout) :: this        
        real, intent(in) :: y(13)
        real :: dy_dt(13)

        real :: mass, I(3,3), F(3), M(3), g, I_inv(3,3), dummy(3), h(3), a_c

        if (verbose) then
            write(*,'(A,13ES20.12)') "y: ", y
        end if

        g = gravity_English(-y(9))


        call this%mass_inertia(y, mass, I)
        call this%pseudo_aero(y, F, M)
        call this%gyroscopic(y, h)

        dy_dt = 0.0

        ! Eq. 5.4.7
        dy_dt(7:9) = quat_dependent_to_base(y(1:3), y(10:13))

        ! Gravity relief
        a_c = (dy_dt(7)**2 + dy_dt(8)**2)/((r_e/0.3048) -y(9))

        ! Sim of Flight Eq. 5.4.5
        dy_dt(1:3) = (1.0/mass)*F
        dy_dt(1) = dy_dt(1) + (g-a_c)*(2*(y(11)*y(13) - y(12)*y(10))) + (y(6)*y(2) - y(5)*y(3))
        dy_dt(2) = dy_dt(2) + (g-a_c)*(2*(y(12)*y(13) + y(11)*y(10))) + (y(4)*y(3) - y(6)*y(1))
        dy_dt(3) = dy_dt(3) + (g-a_c)*(y(13)*y(13) + y(10)*y(10) - y(11)*y(11) - y(12)*y(12)) + (y(5)*y(1) - y(4)*y(2))

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

        write(*,*) "Gyroscope: ", h


        ! Eq. 5.4.8
        dy_dt(10) = 0.5*(- y(11)*y(4) - y(12)*y(5) - y(13)*y(6))
        dy_dt(11) = 0.5*(  y(10)*y(4) - y(13)*y(5) + y(12)*y(6))
        dy_dt(12) = 0.5*(  y(13)*y(4) + y(10)*y(5) - y(11)*y(6))
        dy_dt(13) = 0.5*(- y(12)*y(4) + y(11)*y(5) + y(10)*y(6))
        
        if (verbose) then
            write(*,'(A,13ES20.12)') "dy_dt: ", dy_dt
            write(*,*) "----------------------"
        end if

    end function aircraft_diff_eq

    subroutine aircraft_mass_inertia(this, y, mass, I)
    
        class(aircraft), intent(in) :: this
        real, intent(in) :: y(13)
        real, intent(out) :: mass, I(3,3)

        real :: E(3,3)
       
        mass = this%M/gravity_English(0.0)
        I = 0.0
        I(1,1) = this%Ixx
        I(2,2) = this%Iyy
        I(3,3) = this%Izz

        I(1,2) = -this%Ixy
        I(2,1) = -this%Ixy
        I(2,3) = -this%Iyz
        I(3,2) = -this%Iyz
        I(1,3) = -this%Ixz
        I(3,1) = -this%Ixz

        
    end subroutine aircraft_mass_inertia

    subroutine aircraft_gyroscopic(this, y, h)
    
        class(aircraft), intent(in) :: this
        real, intent(in) :: y(13)
        real, intent(out) :: h(3)
        
        h = this%h
        
    end subroutine aircraft_gyroscopic

    function aircraft_thrust(this, y, throttle) result(thrust)
   
        implicit none
        
        class(aircraft), intent(in) :: this
        real, intent(in) :: y(13), throttle
        
        real :: thrust(3)

        real :: Z, Temp, P, rho, a
        real :: Z_0, Temp_0, P_0, rho_0, a_0

        ! Get atmosphere
        call std_atm_English(-y(9), Z, temp, P, rho, a)
        call std_atm_English(0.0, Z_0, temp_0, P_0, rho_0, a_0)

        thrust = 0.0
        thrust(1) = throttle*this%T0*(rho/rho_0)**this%a
        !if (throttle < 0.0) thrust = 0.0

    end function aircraft_thrust

    !subroutine aircraft_landing_gear(this, t, y, F, M)

        !implicit none
        
        !class(aircraft), intent(in) :: this
        !real, intent(in) :: t, y(13)
        !real, intent(out) :: F(3), M(3)

        !integer :: N, i
        !real :: dummy_F(3)
        !real, allocatable :: g_f(:,:), v_b(:,:), v_f(:,:)

        !N = size(this%ks)

        !allocate(g_f(N, 3))
        !allocate(v_f(N, 3))
        !allocate(v_b(N, 3))

        !! Rotate into earth fixed coordinates
        !do i = 1, N
            !g_f(i,:) = quat_dependent_to_base(this%g_b(i,:), y(10:13)) + y(7:9)
            !v_b(i,1) = y(1) + y(5)*this%g_b(i,3) - y(6)*this%g_b(i,2)
            !v_b(i,2) = y(2) + y(6)*this%g_b(i,1) - y(4)*this%g_b(i,3)
            !v_b(i,3) = y(3) + y(4)*this%g_b(i,2) - y(5)*this%g_b(i,1)
            !v_f(i,:) = quat_dependent_to_base(v_b(i,:), y(10:13))
        !end do

        !F = 0.0
        !M = 0.0
    
        !! Determine spring forces and moments
        !do i = 1, N 
            !dummy_F = 0.0
            !! Only land on deck
            !if (abs(g_f(i,1)) < 400.0 .and. abs(g_f(i,2)) < 90.0) then 
                !if (g_f(i,3) > 0.0) then
                    !dummy_F(3) = - this%ks(i)*g_f(i,3) - this%kd(i)*v_f(i,3)
                    !!dummy_F(1) = -10.0*v_b(i,1) ! Braking force
                    !if (dummy_F(3) > 0.0) dummy_F = 0.0
                    !M(1) = M(1) + (this%g_b(i,2)*dummy_F(3) - this%g_b(i,3)*dummy_F(2))
                    !M(2) = M(2) + (this%g_b(i,3)*dummy_F(1) - this%g_b(i,1)*dummy_F(3))
                    !M(3) = M(3) + (this%g_b(i,1)*dummy_F(2) - this%g_b(i,2)*dummy_F(1))
                    !F = F + dummy_F
                !end if
            !end if
        !end do


    !end subroutine aircraft_landing_gear


    !subroutine aircraft_arresting_gear(this, t, y, F, M)

        !implicit none
        
        !class(aircraft), intent(inout) :: this
        !real, intent(in) :: t, y(13)
        !real, intent(out) :: F(3), M(3)

        !integer :: N, i
        !real :: dummy_F(3)
        !real :: th_f(3), v_b(3), v_f(3)

        !! Rotate into earth fixed coordinates
        !th_f = quat_dependent_to_base(this%th_b(:), y(10:13)) + y(7:9)
        !v_b(1) = y(1) + y(5)*this%th_b(3) - y(6)*this%th_b(2)
        !v_b(2) = y(2) + y(6)*this%th_b(1) - y(4)*this%th_b(3)
        !v_b(3) = y(3) + y(4)*this%th_b(2) - y(5)*this%th_b(1)
        !v_f = quat_dependent_to_base(v_b, y(10:13))

        !! Check if tail hook has caught cable
        !if (.not. this%ag_engaged) then
            !if (th_f(3) > -0.5) then
                !if (abs(th_f(2)) < this%ag_length) then
                    !if (abs(th_f(1) - this%ag_distance) < 0.01) then
                        !this%ag_engaged = .true.
                    !end if
                !end if
            !end if
        !end if

        !F = 0.0
        !M = 0.0

        !! Determine forces and moments
        !if (this%ag_engaged) then
            !F(1) = - this%ag_ks*(th_f(1) - this%ag_distance) - this%ag_kd*v_f(1)
            !M(1) = M(1) + (this%th_b(2)*dummy_F(3) - this%th_b(3)*dummy_F(2))
            !M(2) = M(2) + (this%th_b(3)*dummy_F(1) - this%th_b(1)*dummy_F(3))
            !M(3) = M(3) + (this%th_b(1)*dummy_F(2) - this%th_b(2)*dummy_F(1))
        !end if

    !end subroutine aircraft_arresting_gear


    !function aircraft_check_collision(this, t, y) result(crashed)
        !implicit none
        
        !class(aircraft), intent(in) :: this
        !real, intent(in) :: t, y(13)

        !logical :: crashed

        !integer :: N, i
        !real, allocatable :: r_f(:,:)

        !N = size(this%collision_points, 1)
        !allocate(r_f(N,3))

        !crashed = .false.


        !do i = 1, N
            !r_f(i,:) = quat_dependent_to_base(this%collision_points(i,:), y(10:13)) + y(7:9)
            !if (abs(r_f(i,1)) < 400.0 .and. abs(r_f(i,2)) < 90.0) then 
                !! Crash on deck
                !if (r_f(i,3) > 0.0) then
                    !crashed = .true.
                    !exit
                !end if
            !else
                !! Crash on water
                !if (r_f(i,3) > 60.0) then
                    !crashed = .true.
                    !exit
                !end if
            !end if
        !end do
    !end function aircraft_check_collision

    subroutine aircraft_pseudo_aero(this,y,F,M)

        implicit none
    
        class(aircraft), intent(inout) :: this
        real, intent(in) :: y(13)
        real, intent(out) :: F(3), M(3)

        real :: Z, Temp, P, rho, a, mu, Re
        real :: C_L1, C_L, C_D, C_S, C_ell, C_m, C_n
        real :: de, da, dr, throttle

        real :: alpha, beta, pbar, qbar, rbar, V, alphahat
        real :: S_alpha, C_alpha, S_beta, C_beta
        real :: thrust(3), F_g(3), M_g(3)
        real :: CLnewt, CDnewt, Cmnewt, pos, neg, sigma

        !print*, "State vector incoming: ", y
        da = this%controls(1)
        de = this%controls(2)
        dr = this%controls(3)
        throttle = this%controls(4)

        ! Get atmosphere
        call std_atm_English(-y(9), Z, temp, P, rho, a)
        
        V = sqrt(y(1)**2 + y(2)**2 + y(3)**2)

        pbar = (0.5/V)*y(4)*this%b
        qbar = (0.5/V)*y(5)*this%c
        rbar = (0.5/V)*y(6)*this%b

        alpha = atan2(y(3),y(1))
        alphahat = (this%c/(2.0*V))*y(5)

        ! traditional beta
        beta = asin(y(2)/V)
        ! Flank angle
        !beta = atan2(y(2),y(1))

        if (this%type == 'aircraft') then
            ! Calculate aerodynamic coefficients
            C_L1 = this%CL_0 + this%CL_alpha*alpha
            C_L = C_L1 + this%CL_qbar*qbar + this%CL_alphahat*alphahat + this%CL_de*de
            C_S = this%CS_beta*beta + (this%CS_pbar + this%CS_alpha_pbar*alpha)*pbar &
                    + this%CS_rbar*rbar + this%CS_da*da + this%CS_dr*dr
            C_D = this%CD_L0 + this%CD_CL1*C_L1 + this%CD_CL1_CL1*(C_L1**2) + this%CD_CS_CS*C_S**2 &
                    + (this%CD_qbar + this%CD_alpha_qbar*alpha)*qbar + this%CD_de*de + this%CD_alpha_de*alpha*de + this%CD_de_de*de**2
            C_ell = this%Cell_beta*beta + this%Cell_pbar*pbar + this%Cell_alpha_rbar*alpha*rbar + this%Cell_rbar*rbar &
                    + this%Cell_da*da + this%Cell_dr*dr
            C_m = this%Cm_0 + this%Cm_alpha*alpha + this%Cm_qbar*qbar + this%Cm_alphahat*alphahat + this%Cm_de*de
            C_n = this%Cn_beta*beta + this%Cn_pbar*pbar + this%Cn_alpha_pbar*alpha*pbar + this%Cn_rbar*rbar &
                    + this%Cn_da*da + this%Cn_alpha_da*alpha*da + this%Cn_dr*dr

            S_alpha = sin(alpha)
            C_alpha = cos(alpha)
            S_beta = sin(beta)
            C_beta = cos(beta)

            if (this%include_stall) then
                ! Compute stall CL
                CLnewt  = 2.0*sign(1.0,alpha)*S_alpha*S_alpha*C_alpha
                pos = exp( this%CL_stall%lambda_b*(alpha - this%CL_stall%alpha_0 + this%CL_stall%alpha_s))
                neg = exp(-this%CL_stall%lambda_b*(alpha - this%CL_stall%alpha_0 - this%CL_stall%alpha_s))
                sigma = (1.0 + neg + pos)/((1.0 + neg)*(1.0 + pos))
                C_L = (1.0 - sigma)*C_L + sigma*(CLnewt)
        
                ! Compute stall CD
                CDnewt  = 2.0*sin(abs(alpha))**3
                pos = exp( this%CD_stall%lambda_b*(alpha - this%CD_stall%alpha_0 + this%CD_stall%alpha_s))
                neg = exp(-this%CD_stall%lambda_b*(alpha - this%CD_stall%alpha_0 - this%CD_stall%alpha_s))
                sigma = (1.0 + neg + pos)/((1.0 + neg)*(1.0 + pos))
                C_D = (1.0 - sigma)*C_D + sigma*(CDnewt)

                ! Compute stall Cm
                Cmnewt  = this%Cm_stall%minval*sign(1.0,alpha)*S_alpha*S_alpha
                pos = exp( this%Cm_stall%lambda_b*(alpha - this%Cm_stall%alpha_0 + this%Cm_stall%alpha_s))
                neg = exp(-this%Cm_stall%lambda_b*(alpha - this%Cm_stall%alpha_0 - this%Cm_stall%alpha_s))
                sigma = (1.0 + neg + pos)/((1.0 + neg)*(1.0 + pos))
                C_m = (1.0 - sigma)*C_m + sigma*(Cmnewt)
            end if

            F(1) = 0.5*rho*V**2 * this%S_w * (C_L*S_alpha - C_S*C_alpha*S_beta - C_D*C_alpha*C_beta)
            F(2) = 0.5*rho*V**2 * this%S_w * (C_S*C_beta - C_D*S_beta)
            F(3) = 0.5*rho*V**2 * this%S_w * (-C_L*C_alpha - C_S*S_alpha*S_beta - C_D*S_alpha*C_beta)

            M(1) = 0.5*rho*V**2 * this%S_w * (this%b*(C_ell))
            M(2) = 0.5*rho*V**2 * this%S_w * (this%c*C_m)
            M(3) = 0.5*rho*V**2 * this%S_w * (this%b*(C_n))

        else if (this%type == 'arrow') then

            beta = atan2(y(2),y(1))

            C_L = this%CL_alpha*alpha
            C_S = -this%CL_alpha*beta
            C_D = this%CD_L0 + this%CD_CL1_CL1*C_L**2 + this%CD_CL1_CL1*C_S**2
            
            C_ell = this%Cell_0 + this%Cell_pbar*pbar
            C_m = this%Cm_alpha*alpha + this%Cm_qbar*qbar
            C_n = -this%Cm_alpha*beta + this%Cm_qbar*rbar

            beta = asin(y(2)/V)
            S_alpha = sin(alpha)
            C_alpha = cos(alpha)
            S_beta = sin(beta)
            C_beta = cos(beta)

            F(1) = 0.5*rho*V**2 * this%S_w * (C_L*S_alpha - C_S*C_alpha*S_beta - C_D*C_alpha*C_beta)
            F(2) = 0.5*rho*V**2 * this%S_w * (C_S*C_beta - C_D*S_beta)
            F(3) = 0.5*rho*V**2 * this%S_w * (-C_L*C_alpha - C_S*S_alpha*S_beta - C_D*S_alpha*C_beta)

            M(1) = 0.5*rho*V**2 * this%S_w * (this%b*(C_ell))
            M(2) = 0.5*rho*V**2 * this%S_w * (this%c*C_m)
            M(3) = 0.5*rho*V**2 * this%S_w * (this%b*(C_n))

        else if (this%type == 'sphere') then
            mu = sutherland_visc_English(temp)

            Re = 2.0*rho*norm2(y(1:3))*this%b/mu

            if (Re < 0.01) then
                C_D = 2405.0
            else if (Re >= 0.01 .and. Re <= 450000.) then
                C_D = 24.0/Re + 6.0/(1.0+sqrt(Re)) + 0.4
            else if (Re > 450000. .and. Re <= 560000) then
                C_D = 1.0e29*Re**(-5.211)
            else if (Re > 560000. .and. Re <= 14000000.) then
                C_D = -2.0e-23*Re**3 - 1.e-16*Re**2 + 9.0e-9*Re + 0.069
            else
                C_D = 0.12
            end if


            F = -0.5 * rho * V**2 * PI * this%b**2 * C_D * y(1:3)/V
            M = 0.0

        end if

        ! Apply CG shift to moments
        M(1) = M(1) + (this%CG_shift(2)*F(3) - this%CG_shift(3)*F(2))
        M(2) = M(2) + (this%CG_shift(3)*F(1) - this%CG_shift(1)*F(3))
        M(3) = M(3) + (this%CG_shift(1)*F(2) - this%CG_shift(2)*F(1))

        ! Add thrust influence
        thrust = this%thrust(y, throttle)
        F = F + thrust
        M(1) = M(1) + (this%t_location(2)*thrust(3) - this%t_location(3)*thrust(2))
        M(2) = M(2) + (this%t_location(3)*thrust(1) - this%t_location(1)*thrust(3))
        M(3) = M(3) + (this%t_location(1)*thrust(2) - this%t_location(2)*thrust(1))

        ! Add landing gear influence
        !call this%landing_gear(t, y, F_g, M_g)
        !F = F + F_g
        !M = M + M_g
        !call this%arresting_gear(t, y, F_g, M_g)
        !F = F + F_g
        !M = M + M_g
        
        if (verbose) then
            write(*,'(A,3ES20.12)') "F: ", F
            write(*,'(A,3ES20.12)') "M: ", M
        end if
        
    end subroutine aircraft_pseudo_aero


    function aircraft_tick_states(this, dt) result(new_states)

        implicit none
        
        class(aircraft), intent(inout) :: this
        real, intent(in) :: dt
        real :: new_states(13)


        if (this%run_physics) then
            if (verbose) then
                write(*,*) " State at beginning of RK4."
                write(*,'(A,13ES20.12)') "y: ", this%states
                write(*,*) "----------------------"
            end if

            this%states = this%runge_kutta(this%states, dt)
            if (verbose) then
                write(*,*) " State at end of RK4."
                write(*,'(A,13ES20.12)') "y: ", this%states
                write(*,*) "----------------------"
            end if
            call quat_norm(this%states(10:13))

        end if

    end function aircraft_tick_states


    subroutine aircraft_print_aero_table(this, V, H)

        implicit none
        class(aircraft), intent(inout) :: this
        real, intent(in) :: V, H
        integer :: i, iunit
        real :: alpha, beta, states(13), F(3), M(3), Ax, N, Y
        real :: ca, cb, sa, sb
        real :: CL, CD, Cm
        real :: Z,T,P,rho,a, mu, const
        real :: cntrl(4)

        call std_atm_English(H, Z, T, P, rho, a)
        const = (0.5*rho*V**2*this%S_w)

        cntrl = 0.0
        states = 0.0


        open(newunit=iunit, file='aero_table.csv', status='REPLACE')
        write(iunit,*) 'alpha[deg], CL, CD, Cm'
        do i = -180,180,1

            alpha = real(i)*PI/180.0
            beta = 0.0

            states(1) = V*cos(alpha)*cos(beta)
            states(2) = V*sin(beta)
            states(3) = V*sin(alpha)*cos(beta)
            states(9) = -H

            call this%pseudo_aero(states, F, M)
            Ax = -F(1)
            Y = F(2)
            N = -F(3)

            ca = cos(alpha)
            cb = cos(beta)
            sa = sin(alpha)
            sb = sin(beta)

            CL = N*ca - Ax*sa
            CD = Ax*ca*cb - Y*sb + N*sa*cb
            Cm = M(2)

            CL = CL/const
            CD = CD/const
            Cm = cm/(const*this%c)

            write(iunit, *) alpha*180./PI, ',', CL,',',CD,',',Cm
        end do
    end subroutine aircraft_print_aero_table
    
end module vehicle_m 

