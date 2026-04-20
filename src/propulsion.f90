module propulsion_m
    use goates_m
    use atmosphere_m
    use jsonx_m
    implicit none
    
    type propulsion_t
        character(len=:), allocatable :: name, type, units
        real, allocatable :: location(:), orientation_eul(:)
        real :: orientation_quat(4)
        integer :: control_ID


        ! For type = T=f(v)
        real, allocatable :: T_coeffs(:)
        real :: Ta

        ! For type = propeller_polynomial
        real :: diameter, Ixx
        integer :: rotation_delta
        real, allocatable :: CT_J(:), CP_J(:), CNa_J(:), Cnna_J(:)

    end type propulsion_t

contains

    subroutine propulsion_init(this, j_propulsion)

        implicit none

        type(propulsion_t), intent(inout) :: this
        type(json_value), pointer :: j_propulsion
        integer :: i
        character(len=:), allocatable :: temp

        this%name = j_propulsion%name
        write(*,*) "       Initilizing Propulsion : ", trim(this%name)

        call jsonx_get(j_propulsion, 'location[ft]', this%location, 0.0, 3)
        call jsonx_get(j_propulsion, 'orientation[deg]', this%orientation_eul, 0.0, 3)
        this%orientation_eul = this%orientation_eul*pi/180.0
        this%orientation_quat = euler_to_quat(this%orientation_eul)

        call jsonx_get(j_propulsion, 'type', this%type)
        select case (this%type)
        case("T=f(V)")
            call jsonx_get(j_propulsion, "T_coefficients[lbf]", this%T_coeffs, 0.0)
            call jsonx_get(j_propulsion, "Ta", this%Ta)
            this%rotation_delta = 1
            this%Ixx = 0.0
        case("propeller_polynomial")
            call jsonx_get(j_propulsion, 'diameter[ft]', this%diameter)
            call jsonx_get(j_propulsion, 'Ixx[slug-ft^2]', this%Ixx)
            call jsonx_get(j_propulsion, 'rotation', temp)
            this%rotation_delta = 1
            if(trim(temp) == "LH") this%rotation_delta = -1
            call jsonx_get(j_propulsion, 'CT(J)', this%CT_J, 0.0)
            call jsonx_get(j_propulsion, 'CPb(J)', this%CP_J, 0.0)
            call jsonx_get(j_propulsion, 'CN,alpha(J)', this%CNa_J, 0.0)
            call jsonx_get(j_propulsion, 'Cn,alpha(J)', this%Cnna_J, 0.0)
        end select
    end subroutine propulsion_init


    function propulsion_get_FMh(this, states, tau) result(ans)
        implicit none
        type(propulsion_t), intent(in) :: this
        real, intent(in) :: states(21), tau
        real :: ans(9)
        real :: Vc(3), Vc_mag, uc(3), vN(3), vN_mag, uN(3)
        real :: Fc(3), Mc(3)
        real :: Thrust, Normal, Torque, Yaw, hxx, Torque_s, Torque_m
        real :: Z_dum, T_dum, P_dum, rho, rho0, a_dum, mu_dum, dyp
        real :: Hz, omega, J, alpha_c, N_r, N_m, error, I_m, eta_c, B, C, E_m, E_b, I_b


        Vc = quat_base_to_dependent(states(1:3)+cross_product(states(4:6), this%location), this%orientation_quat)
        Vc_mag = norm2(Vc)
        uc = Vc/Vc_mag
        if (Vc_mag < 1e-12) uc = [1.0, 0.0, 0.0]
        vN = -[0.0, uc(2), uc(3)]
        vN_mag = norm2(vN)
        uN = vN/vN_mag
        if (vN_mag < 1e-12) uN = [0.0, 0.0, 1.0]
        alpha_c = acos(uc(1))
        !write(*,*) "alpha_c: ", alpha_c*180.0/pi

        call std_atm_English(0.0, Z_dum, T_dum, P_dum, rho0, a_dum)
        call std_atm_English(-states(9), Z_dum, T_dum, P_dum, rho, a_dum)

        select case (this%type)
            case("T=f(V)")
                Thrust = tau*calc_polynomial(this%T_coeffs, Vc_mag)*(rho/rho0)**this%Ta
                Normal = 0.0
                Torque = 0.0
                Yaw = 0.0
                hxx = 0.0
            case("propeller_polynomial")
                select case(this%units)
                case("RPM")
                    Hz = tau/60.0
                    omega = Hz*2*pi
                case("thrust[lbf]")
                    omega = -Vc_mag*this%CT_J(1) + sqrt(Vc_mag**2*this%CT_J(2)**2 - 4.0*this%CT_J(1)*(Vc_mag**2*this%CT_J(3) - tau/rho/this%diameter**2))
                    omega = omega*PI/this%diameter/this%CT_J(1)
                    Hz = omega/2.0/pi
                case("torque[ft-lbf]")
                    omega = -Vc_mag*this%CP_J(2) + sqrt(Vc_mag**2*this%CP_J(2)**2 - 4.0*this%CP_J(1)*(Vc_mag**2*this%CP_J(3) - tau*2.0*pi/rho/this%diameter**3))
                    omega = omega*PI/this%diameter/this%CP_J(1)
                    Hz = omega/2.0/pi
                case("throttle")
                    N_r = 1000
                    error = 1.0
                    do while(abs(error) > 1e-12)
                        write(*,*) "N_r [rpm]: ", N_r
                        Hz = N_r/60.0
                        omega = Hz*2*pi
                        J = 2.0*PI*Vc_mag/omega/this%diameter
                        Torque = calc_polynomial(this%CP_J, J) * rho *(Hz**3)*(this%diameter**5) / omega
                        Torque_s = Torque
                        Torque_m = Torque_s ! Gear ratio of 1
                        I_m = Torque_m*450/7.04319971369755 + 0.3
                        eta_c = 1-0.078*(1-tau)
                        B = 2*I_m*0.01 - tau*22.2 + tau**2*I_m*0.025/eta_c
                        C = I_m**2*0.01**2 - tau*22.2*I_m*0.01
                        E_m = 0.5*(-B + sqrt(B**2 - 4*C))
                        N_m = 450*(E_m - I_m*0.12)
                        write(*,*) "N_s [rpm]: ", N_m
                        error = N_m - N_r
                        write(*,*) "error: ", error
                        N_r = N_r + 0.9*(N_m - N_r)
                    end do

                    E_b = 0.5*(22.2 + sqrt(22.2**2 - 4*E_m*I_m*0.025/eta_c))
                    I_b = E_m*I_m/(eta_c*E_b)
                    
                    write(*,*) "E_m: ", E_m, "volts"
                    write(*,*) "I_m: ", I_m, "amps"
                    write(*,*) "E_b: ", E_b, "volts"
                    write(*,*) "I_b: ", I_b, "amps"
                end select

                J = 2.0*PI*Vc_mag/omega/this%diameter
                !write(*,*) "Omega : ", omega
                !write(*,*) "J : ", J

                Thrust = calc_polynomial(this%CT_J, J) * rho *(Hz**2)*(this%diameter**4)
                Torque = calc_polynomial(this%CP_J, J) * rho *(Hz**3)*(this%diameter**5) / omega
                Normal = calc_polynomial(this%CNa_J, J) * rho *(Hz**2)*(this%diameter**4)*alpha_c
                Yaw = calc_polynomial(this%Cnna_J, J) * rho *(Hz**2)*(this%diameter**5)*alpha_c
                hxx = this%rotation_delta*this%Ixx*omega
                if ((omega < 1e-12) .or. (omega .ne. omega)) then
                    Thrust = 0.0
                    Torque = 0.0
                    Normal = 0.0
                    Yaw = 0.0
                    hxx = 0.0
                end if
                write(*,*) " Thrust: ", Thrust
                write(*,*) " Torque: ", Torque
                write(*,*) " Normal: ", Normal
                write(*,*) " Yaw: ", Yaw

        end select

        Fc = [Thrust, 0.0, 0.0] + Normal*uN
        Mc = -real(this%rotation_delta)*([Torque, 0.0, 0.0] + Yaw*uN)
        !write(*,*) "tau: ", tau
        !write(*,*) " Fc: ", Fc
        !write(*,*) " Mc: ", Mc
        !write(*,*) " hxx: ", hxx

        ans(1:3) = quat_dependent_to_base(Fc, this%orientation_quat)
        ans(4:6) = quat_dependent_to_base(Mc, this%orientation_quat) + cross_product(this%location, ans(1:3))
        ans(7:9) = quat_dependent_to_base([hxx, 0.0, 0.0], this%orientation_quat)

        !write(*,*) " Fb: ", ans(1:3)
        !write(*,*) " Mb: ", ans(4:6)
        !write(*,*) " hxxb: ", ans(7:9)

    end function propulsion_get_FMh

    function calc_polynomial(coeffs, var) result(ans)
        implicit none
        real, intent(in) :: coeffs(:), var
        real :: ans
        integer :: i

        ans = 0.0
        do i = 1, size(coeffs)
            ans = ans + coeffs(i)*var**(i-1)
        end do
    end function calc_polynomial

end module propulsion_m
