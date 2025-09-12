module aircraft_m
    use goates_m
    use jsonx_m

    implicit none

    
    type :: aircraft
        
        real :: M, Ixx, Iyy, Izz, Ixy, Iyz, Ixz ! Mass and Inertia matrix
        real :: S_w, b, c ! Reference area, span, chord
        real :: CL_0, CL_alpha, CL_alphahat, CL_qbar, CL_de ! Lift coefficients
        real :: CS_beta, CS_pbar, CS_alpha_pbar, CS_rbar, CS_da, CS_dr ! Side force coefficients
        real :: CD_L0, CD_CL1, CD_CL1_CL1, CD_CS_CS, CD_qbar, CD_alpha_qbar, CD_de, CD_alpha_de, CD_de_de ! Drag coefficients
        real :: Cell_beta, Cell_pbar, Cell_alpha_rbar, Cell_rbar, Cell_da, Cell_dr ! Rolling moment coefficients
        real :: Cm_0, Cm_alpha, Cm_qbar, Cm_alphahat, Cm_de ! Pitching moment coefficients
        real :: Cn_beta, Cn_pbar, Cn_alpha_pbar, Cn_rbar, Cn_da, Cn_alpha_da, Cn_dr ! Yawing moment coefficients
    
    
    contains

        procedure :: mass_inertia => aircraft_mass_inertia
        procedure :: aerodynamics => aircraft_aerodynamics
        procedure :: init => aircraft_init
        
    
    end type aircraft
    
contains

    subroutine aircraft_init(this, settings)
    
        class(aircraft), intent(inout) :: this
        type(json_value), pointer, intent(in) :: settings

        type(json_value), pointer :: reference, coefficients

        ! Parse JSON settings
        call jsonx_get(settings, "aerodynamics.reference", reference)
        call jsonx_get(settings, "aerodynamics.coefficients", coefficients)

        ! Mass properties
        call jsonx_get(settings, "mass.weight[lbf]", this%M)
        call jsonx_get(settings, "mass.Ixx[slug-ft^2]", this%Ixx, 0.0)
        call jsonx_get(settings, "mass.Iyy[slug-ft^2]", this%Iyy, 0.0)
        call jsonx_get(settings, "mass.Izz[slug-ft^2]", this%Izz, 0.0)
        call jsonx_get(settings, "mass.Ixy[slug-ft^2]", this%Ixy, 0.0)
        call jsonx_get(settings, "mass.Iyz[slug-ft^2]", this%Iyz, 0.0)
        call jsonx_get(settings, "mass.Ixz[slug-ft^2]", this%Ixz, 0.0)

        ! Reference properties
        call jsonx_get(reference, "area[ft^2]", this%S_w)
        call jsonx_get(reference, "longitudinal_length[ft]", this%c)
        call jsonx_get(reference, "lateral_length[ft]", this%b)

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
        real, intent(out) :: mass
        real, intent(out) :: I(3,3)
        
        mass = 0.0
        
        I = 0.0
        
    end subroutine aircraft_mass_inertia

    subroutine aircraft_aerodynamics(this, t, y, F, M)
    
        class(aircraft), intent(in) :: this
        real, intent(in) :: t, y(13)
        real, intent(out) :: F(3), M(3)
        
        real :: Z, Temp, P, rho, a
        real :: C_L1, C_L, C_D, C_S, C_ell, C_m, C_n
        real :: de, da, dr

        real :: alpha, beta, pbar, qbar, rbar, V, alphahat
        real :: S_alpha, C_alpha, S_beta, C_beta

        !print*, "State vector incoming: ", y
        de = 0.0
        da = 0.0
        dr = 0.0

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
        C_D = this%CD_L0 + this%CD_CL1*C_L1 + this%CD_CL1_CL1*C_L1**2 + this%CD_CS_CS*C_S**2 &
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

        F(1) = 0.5*rho*V**2 * this%S_w * (C_L*S_alpha - C_S*C_alpha*S_beta -C_D*C_alpha*C_beta)
        F(2) = 0.5*rho*V**2 * this%S_w * (C_S*C_beta + C_D*S_beta)
        F(3) = 0.5*rho*V**2 * this%S_w * (-C_L*C_alpha - C_S*S_alpha*S_beta - C_D*S_alpha*C_beta)

        M(1) = 0.5*rho*V**2 * this%S_w * (this%b*(C_ell*C_alpha*C_beta - C_n*S_alpha) - this%c*C_m*C_alpha*S_beta)
        M(2) = 0.5*rho*V**2 * this%S_w * (this%b*C_ell*S_beta + this%c*C_m*C_beta)
        M(3) = 0.5*rho*V**2 * this%S_w * (this%b*(C_ell*S_alpha*C_beta + C_n*C_alpha) - this%c*C_m*S_alpha*S_beta)
        
        write(*,*) "F: ", F
        write(*,*) "M: ", M
        
    end subroutine aircraft_aerodynamics
    
end module aircraft_m

