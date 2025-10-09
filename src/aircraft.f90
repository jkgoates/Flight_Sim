module aircraft_m
    use goates_m
    use jsonx_m

    implicit none

    
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
    
    
    contains

        procedure :: mass_inertia   => aircraft_mass_inertia
        procedure :: aerodynamics   => aircraft_aerodynamics
        procedure :: gyroscopic     => aircraft_gyroscopic
        procedure :: thrust         => aircraft_thrust
        procedure :: init           => aircraft_init
        
    
    end type aircraft
    
contains

    subroutine aircraft_init(this, settings)
    
        class(aircraft), intent(inout) :: this
        type(json_value), pointer, intent(in) :: settings

        type(json_value), pointer :: reference, coefficients, aerodynamics
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

    subroutine aircraft_aerodynamics(this, t, y, F, M, controls)

        implicit none
    
        class(aircraft), intent(in) :: this
        real, intent(in) :: t, y(13)
        real, intent(out) :: F(3), M(3), controls(4)

        real :: Z, Temp, P, rho, a
        real :: C_L1, C_L, C_D, C_S, C_ell, C_m, C_n
        real :: de, da, dr, throttle

        real :: alpha, beta, pbar, qbar, rbar, V, alphahat
        real :: S_alpha, C_alpha, S_beta, C_beta
        real :: thrust(3)

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
        
        if (verbose) then
            write(*,'(A,3ES20.12)') "F: ", F
            write(*,'(A,3ES20.12)') "M: ", M
        end if
        
    end subroutine aircraft_aerodynamics
    
end module aircraft_m

