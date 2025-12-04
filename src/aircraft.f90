module aircraft_m
    use goates_m
    use jsonx_m

    implicit none

    type stall_settings_t
        real :: alpha_0, alpha_s, lambda_b, minval        
    end type stall_settings_t
    
    type :: aircraft
        
        real :: M, Ixx, Iyy, Izz, Ixy, Iyz, Ixz ! Mass and Inertia matrix
        real, dimension(:), allocatable :: CG_shift ! CG shift from reference point (ft)
        real :: hx, hy, hz ! Gyroscopic components
        real :: T0, a ! Thrust parameters
        real, dimension(:), allocatable :: t_location
        real :: S_w, b, c ! Reference area, span, chord
        real :: CL_0, CL_alpha, CL_alphahat, CL_qbar, CL_de ! Lift coefficients
        real :: CS_beta, CS_pbar, CS_alpha_pbar, CS_rbar, CS_da, CS_dr ! Side force coefficients
        real :: CD_L0, CD_CL1, CD_CL1_CL1, CD_CS_CS, CD_qbar, CD_alpha_qbar, CD_de, CD_alpha_de, CD_de_de ! Drag coefficients
        real :: Cell_beta, Cell_pbar, Cell_alpha_rbar, Cell_rbar, Cell_da, Cell_dr ! Rolling moment coefficients
        real :: Cm_0, Cm_alpha, Cm_qbar, Cm_alphahat, Cm_de ! Pitching moment coefficients
        real :: Cn_beta, Cn_pbar, Cn_alpha_pbar, Cn_rbar, Cn_da, Cn_alpha_da, Cn_dr ! Yawing moment coefficients
        ! real :: rho0

        real, allocatable :: ks(:), kd(:) ! Landing gear spring and damping coefficients
        real, allocatable :: g_b(:,:) ! Landing gear body fixed coordinates
        real, allocatable :: th_b(:) ! tail hook location
        real :: ag_ks, ag_kd, ag_length, ag_distance ! arresting gear coefficients
        logical :: ag_engaged
        real, allocatable :: collision_points(:,:) ! Aircraft collision points
    
        type(stall_settings_t) :: CD_stall, CL_stall, Cm_stall
    
    contains

        procedure :: mass_inertia   => aircraft_mass_inertia
        procedure :: aerodynamics   => aircraft_aerodynamics
        procedure :: gyroscopic     => aircraft_gyroscopic
        procedure :: thrust         => aircraft_thrust
        procedure :: init           => aircraft_init
        procedure :: landing_gear   => aircraft_landing_gear
        procedure :: check_collision => aircraft_check_collision
        procedure :: arresting_gear => aircraft_arresting_gear
        
    
    end type aircraft
    
contains

    subroutine aircraft_init(this, settings)
    
        class(aircraft), intent(inout) :: this
        type(json_value), pointer, intent(in) :: settings

        type(json_value), pointer :: reference, coefficients, aerodynamics, p1, p2
        real, allocatable :: dummy_loc(:)
        integer :: N, cnt, i
        logical :: found

        ! Parse JSON settings
        call jsonx_get(settings, "aerodynamics", aerodynamics)
        call jsonx_get(aerodynamics, "coefficients", coefficients)
        call jsonx_get(aerodynamics, "reference", reference)

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
        call jsonx_get(settings, "mass.hx[slug-ft^2/s]", this%hx, 0.0)
        call jsonx_get(settings, "mass.hy[slug-ft^2/s]", this%hy, 0.0)
        call jsonx_get(settings, "mass.hz[slug-ft^2/s]", this%hz, 0.0)

        ! Thrust Properties
        call jsonx_get(settings, "thrust.T0[lbf]", this%T0, 0.0)
        call jsonx_get(settings, "thrust.a", this%a, 0.0)
        call jsonx_get(settings, "thrust.location[ft]", this%t_location, 0.0, 3)

        ! Reference properties
        call jsonx_get(reference, "area[ft^2]", this%S_w)
        call jsonx_get(reference, "longitudinal_length[ft]", this%c)
        call jsonx_get(reference, "lateral_length[ft]", this%b)
        call jsonx_get(reference, "relative_location[ft]", this%CG_shift, 0.0, 3)

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

        ! Landing Gear
        call jsonx_get(settings, "landing_gear", p1)
        N = json_value_count(p1)
        cnt= 0
        do i = 1, N
            call json_value_get(p1, i, p2)
            if (p2%name(1:1) == 'x') cycle
            cnt = cnt+1
        end do
        allocate(this%ks(cnt))
        allocate(this%kd(cnt))
        allocate(this%g_b(cnt, 3))
        cnt = 1
        do i = 1, N
            call json_value_get(p1, i, p2)
            if (p2%name(1:1) == 'x') cycle
            call jsonx_get(p2, "location[ft]", dummy_loc, 0.0, 3)
            this%g_b(cnt, :) = dummy_loc
            call jsonx_get(p2, "spring_constant[lb/ft]", this%ks(cnt))
            call jsonx_get(p2, "damping_constant[lb-s/ft]", this%kd(cnt))
            cnt = cnt+1
        end do

        ! Collision points
        call jsonx_get(settings, "collision_points", p1)
        N = json_value_count(p1)
        cnt = 0
        do i = 1, N
            call json_value_get(p1, i, p2)
            if (p2%name(1:1) == 'x') cycle
            cnt = cnt+1
        end do
        allocate(this%collision_points(cnt,3))
        cnt = 1
        do i = 1, N
            call json_value_get(p1, i, p2)
            if (p2%name(1:1) == 'x') cycle
            call jsonx_get(p1, p2%name, dummy_loc, 0.0, 3)
            this%collision_points(cnt, :) = dummy_loc
            cnt = cnt+1
        end do

        ! Arrestor Gear
        call jsonx_get(settings, "arresting_gear", p1)
        call jsonx_get(p1, "spring_constant[lb/ft]", this%ag_ks)
        call jsonx_get(p1, "damping_constant[lb-s/ft]", this%ag_kd)
        call jsonx_get(p1, "cable_length[ft]", this%ag_length)
        call jsonx_get(p1, "cable_distance[ft]", this%ag_distance)
        call jsonx_get(p1, "tail_hook[ft]", this%th_b)
        this%ag_engaged = .false.


    end subroutine aircraft_init

    subroutine aircraft_mass_inertia(this, t, y, mass, I)
    
        class(aircraft), intent(in) :: this
        real, intent(in) :: t, y(13)
        real, intent(out) :: mass, I(3,3)

        real :: E(3,3)
       
        !mass = this%M*0.3048/g_ssl ! Convert lbf to slugs
        !mass = this%M/32.1740485564304 ! Convert lbf to slugs
        mass = this%M/gravity_English(0.0)
        I(1,1) = this%Ixx
        I(2,2) = this%Iyy
        I(3,3) = this%Izz

        I(1,2) = -this%Ixy
        I(2,1) = -this%Ixy
        I(2,3) = -this%Iyz
        I(3,2) = -this%Iyz
        I(1,3) = -this%Ixz
        I(3,1) = -this%Ixz

        ! identity matrix
        E = 0.0
        E(1,1) = 1.0
        E(2,2) = 1.0
        E(3,3) = 1.0

        ! Apply CG shift to inertia matrix
        ! THIS EQUATION IS WRONG BUT IS USED TO CHECK BOOK EXAMPLES
        !I = I + mass * (dot_product(this%CG_shift, this%CG_shift)*E &
                        !- matmul(reshape(this%CG_shift, [3,1]), reshape(this%CG_shift, [1,3])))
        ! THIS IS THE CORRECT EQUATION
        !I = I - mass * (dot_product(this%CG_shift, this%CG_shift)*E &
                        !- matmul(reshape(this%CG_shift, [3,1]), reshape(this%CG_shift, [1,3])))
        
    end subroutine aircraft_mass_inertia

    subroutine aircraft_gyroscopic(this, t, y, h)
    
        class(aircraft), intent(in) :: this
        real, intent(in) :: t, y(13)
        real, intent(out) :: h(3)
        
        h(1) = this%hx
        h(2) = this%hy
        h(3) = this%hz
        
    end subroutine aircraft_gyroscopic

    function aircraft_thrust(this, t, y, throttle) result(thrust)
   
        implicit none
        
        class(aircraft), intent(in) :: this
        real, intent(in) :: t, y(13), throttle
        
        real :: thrust(3)

        real :: Z, Temp, P, rho, a
        real :: Z_0, Temp_0, P_0, rho_0, a_0

        ! Get atmosphere
        call std_atm_English(-y(9), Z, temp, P, rho, a)
        call std_atm_English(0.0, Z_0, temp_0, P_0, rho_0, a_0)

        thrust = 0.0
        thrust(1) = throttle*this%T0*(rho/rho_0)**this%a
        if (throttle < 0.0) thrust = 0.0

    end function aircraft_thrust

    subroutine aircraft_landing_gear(this, t, y, F, M)

        implicit none
        
        class(aircraft), intent(in) :: this
        real, intent(in) :: t, y(13)
        real, intent(out) :: F(3), M(3)

        integer :: N, i
        real :: dummy_F(3)
        real, allocatable :: g_f(:,:), v_b(:,:), v_f(:,:)

        N = size(this%ks)

        allocate(g_f(N, 3))
        allocate(v_f(N, 3))
        allocate(v_b(N, 3))

        ! Rotate into earth fixed coordinates
        do i = 1, N
            g_f(i,:) = quat_dependent_to_base(this%g_b(i,:), y(10:13)) + y(7:9)
            v_b(i,1) = y(1) + y(5)*this%g_b(i,3) - y(6)*this%g_b(i,2)
            v_b(i,2) = y(2) + y(6)*this%g_b(i,1) - y(4)*this%g_b(i,3)
            v_b(i,3) = y(3) + y(4)*this%g_b(i,2) - y(5)*this%g_b(i,1)
            v_f(i,:) = quat_dependent_to_base(v_b(i,:), y(10:13))
        end do

        F = 0.0
        M = 0.0
    
        ! Determine spring forces and moments
        do i = 1, N 
            dummy_F = 0.0
            if (g_f(i,3) > 0.0) then
                dummy_F(3) = - this%ks(i)*g_f(i,3) - this%kd(i)*v_f(i,3)
                !dummy_F(1) = -10.0*v_b(i,1) ! Braking force
                if (dummy_F(3) > 0.0) dummy_F = 0.0
                M(1) = M(1) + (this%g_b(i,2)*dummy_F(3) - this%g_b(i,3)*dummy_F(2))
                M(2) = M(2) + (this%g_b(i,3)*dummy_F(1) - this%g_b(i,1)*dummy_F(3))
                M(3) = M(3) + (this%g_b(i,1)*dummy_F(2) - this%g_b(i,2)*dummy_F(1))
                F = F + dummy_F
            end if
        end do


    end subroutine aircraft_landing_gear


    subroutine aircraft_arresting_gear(this, t, y, F, M)

        implicit none
        
        class(aircraft), intent(inout) :: this
        real, intent(in) :: t, y(13)
        real, intent(out) :: F(3), M(3)

        integer :: N, i
        real :: dummy_F(3)
        real :: th_f(3), v_b(3), v_f(3)

        ! Rotate into earth fixed coordinates
        th_f = quat_dependent_to_base(this%th_b(:), y(10:13)) + y(7:9)
        v_b(1) = y(1) + y(5)*this%th_b(3) - y(6)*this%th_b(2)
        v_b(2) = y(2) + y(6)*this%th_b(1) - y(4)*this%th_b(3)
        v_b(3) = y(3) + y(4)*this%th_b(2) - y(5)*this%th_b(1)
        v_f = quat_dependent_to_base(v_b, y(10:13))

        ! Check if tail hook has caught cable
        if (.not. this%ag_engaged) then
            if (th_f(3) > -0.5) then
                if (abs(th_f(2)) < this%ag_length) then
                    if (abs(th_f(1) - this%ag_distance) < 0.01) then
                        this%ag_engaged = .true.
                    end if
                end if
            end if
        end if

        F = 0.0
        M = 0.0

        ! Determine forces and moments
        if (this%ag_engaged) then
            F(1) = - this%ag_ks*(th_f(1) - this%ag_distance) - this%ag_kd*v_f(1)
            M(1) = M(1) + (this%th_b(2)*dummy_F(3) - this%th_b(3)*dummy_F(2))
            M(2) = M(2) + (this%th_b(3)*dummy_F(1) - this%th_b(1)*dummy_F(3))
            M(3) = M(3) + (this%th_b(1)*dummy_F(2) - this%th_b(2)*dummy_F(1))
        end if

    end subroutine aircraft_arresting_gear


    function aircraft_check_collision(this, t, y) result(crashed)
        implicit none
        
        class(aircraft), intent(in) :: this
        real, intent(in) :: t, y(13)

        logical :: crashed

        integer :: N, i
        real, allocatable :: r_f(:,:)

        N = size(this%collision_points, 1)
        allocate(r_f(N,3))

        crashed = .false.

        do i = 1, N
            r_f(i,:) = quat_dependent_to_base(this%collision_points(i,:), y(10:13)) + y(7:9)
            if (r_f(i,3) > 0.0) then
                crashed = .true.
                exit
            end if
        end do
    end function aircraft_check_collision

    subroutine aircraft_aerodynamics(this, t, y, F, M, controls)

        implicit none
    
        class(aircraft), intent(inout) :: this
        real, intent(in) :: t, y(13)
        real, intent(out) :: F(3), M(3), controls(4)

        real :: Z, Temp, P, rho, a
        real :: C_L1, C_L, C_D, C_S, C_ell, C_m, C_n
        real :: de, da, dr, throttle

        real :: alpha, beta, pbar, qbar, rbar, V, alphahat
        real :: S_alpha, C_alpha, S_beta, C_beta
        real :: thrust(3), F_g(3), M_g(3)

        !print*, "State vector incoming: ", y
        da = controls(1)
        de = controls(2)
        dr = controls(3)
        throttle = controls(4)

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

        ! sompute stall CL
        !CLnewt  = 2.0*sign(alpha)*S_alpha*S_alpha*C_alpha
        !pos = exp()
        !neg = exp()
        !sigma = (1.0 + neg + pos)/((1.0 + neg)*(1.0 + pos))
        !C_L = (1.0 - sigma)*C_L + sigma*(CLnewt)

        F(1) = 0.5*rho*V**2 * this%S_w * (C_L*S_alpha - C_S*C_alpha*S_beta - C_D*C_alpha*C_beta)
        F(2) = 0.5*rho*V**2 * this%S_w * (C_S*C_beta - C_D*S_beta)
        F(3) = 0.5*rho*V**2 * this%S_w * (-C_L*C_alpha - C_S*S_alpha*S_beta - C_D*S_alpha*C_beta)

        M(1) = 0.5*rho*V**2 * this%S_w * (this%b*(C_ell))
        M(2) = 0.5*rho*V**2 * this%S_w * (this%c*C_m)
        M(3) = 0.5*rho*V**2 * this%S_w * (this%b*(C_n))

        ! Apply CG shift to moments
        M(1) = M(1) + (this%CG_shift(2)*F(3) - this%CG_shift(3)*F(2))
        M(2) = M(2) + (this%CG_shift(3)*F(1) - this%CG_shift(1)*F(3))
        M(3) = M(3) + (this%CG_shift(1)*F(2) - this%CG_shift(2)*F(1))

        ! Add thrust influence
        thrust = this%thrust(t, y, throttle)
        F = F + thrust
        M(1) = M(1) + (this%t_location(2)*thrust(3) - this%t_location(3)*thrust(2))
        M(2) = M(2) + (this%t_location(3)*thrust(1) - this%t_location(1)*thrust(3))
        M(3) = M(3) + (this%t_location(1)*thrust(2) - this%t_location(2)*thrust(1))

        ! Add landing gear influence
        call this%landing_gear(t, y, F_g, M_g)
        F = F + F_g
        M = M + M_g
        call this%arresting_gear(t, y, F_g, M_g)
        F = F + F_g
        M = M + M_g
        
        if (verbose) then
            write(*,'(A,3ES20.12)') "F: ", F
            write(*,'(A,3ES20.12)') "M: ", M
        end if
        
    end subroutine aircraft_aerodynamics
    
end module aircraft_m

