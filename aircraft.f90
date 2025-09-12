module aircraft_m
    use goates_m
    use jsonx_m
    implicit none

    
    type :: aircraft
        
        real :: M, I(3,3) ! Mass and Inertia matrix
        real :: S_w, b, c ! Reference area, span, chord
        real :: CL_0, CL_alpha, CL_alphahat, CL_qbar, CL_de ! Lift coefficients
        real :: CS_beta, CS_pbar, CS_alpha_pbar, CS_rbar, CS_da, CS_dr ! Side force coefficients
        real :: CD_L0, CD_CL1, CD_CL1_CL1, CD_CS_CS, CD_qbar, CD_alpha_qbar, CD_de, CD_alpha_de, CD_de_de ! Drag coefficients
        real :: Cell_beta, Cell_pbar, Cell_alpha_rbar, Cell_rbar, Cell_da, Cell_dr ! Rolling moment coefficients
        real :: Cm_0, Cm_alpha, Cm_qbar, Cm_alphahat, Cm_de ! Pitching moment coefficients
        real :: Cn_beta, Cn_pbar, Cn_alpha_rbar, Cn_rbar, Cn_da, Cn_alpha_da, Cn_dr ! Yawing moment coefficients
    
    
    contains

        procedure :: mass_inertia => aircraft_mass_inertia
        procedure :: aerodynamics => aircraft_aerodynamics
        
    
    end type aircraft
    
contains

    subroutine aircraft_init(this, settings)
    
        class(aircraft), intent(inout) :: this
        type(json_value), pointer, intent(in) :: settings


        
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
        
        F = 0.0
        M = 0.0
        
    end subroutine aircraft_aerodynamics
    
end module aircraft_m

