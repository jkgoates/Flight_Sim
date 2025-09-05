module goates_m

    implicit none

    ! CONSTANTS
    real,parameter :: PI = 3.141592653589793
    real,parameter :: PI2 = PI*0.5
    real,parameter :: PI4 = PI*0.25
    real,parameter :: g_ssl = 9.80665 ! Standard Gravity at sea level 
    real,parameter :: r_ez = 6356766.0 ! US standard atmosphere earth radius in meters
    real,parameter :: R_gas = 287.0528
    real,parameter :: gamma = 1.4

    ! CONVERSION FACTORS
    real,parameter :: Pa_to_lbf_ft2 = 1.0/47.880258
    real,parameter :: kg_m3_to_slug_ft3 = 1.0/515.379   

    ! ATMOSPHERE TABLES
    real, parameter :: Z_i(8)= [0., 11000., 20000., 32000., 47000., 52000., 61000., 79000.]
    real, parameter :: Z_i1(8) = [11000., 20000., 32000., 47000., 52000., 61000., 79000., 90000.]
    real, parameter :: T_i(8)= [288.150, 216.650, 216.650, 228.650, 270.650, 270.650, 252.650, 180.650]
    real, parameter :: Tp_i(8) = [-6.5, 0.0, 1.0, 2.8, 0.0, -2.0, -4.0, 0.0]
    real, parameter :: P_i(8)= [1.01325e5, 2.26320318222212e4, 5.47487352827083e3, 8.68014769086723e2, &
                1.10905588989225e2, 5.90005242789244e1, 1.82099249050177e1, 1.03770045489203]

    ! OTHER STUFF
    
contains
    
    !!! COORDINATE TRANSFORM !!!

    function quat_mult(q1, q2) result(q3)

        implicit none
        
        real, intent(in) :: q1(4), q2(4)
        real :: q3(4)

        q3(1) = q1(1)*q2(1) - q1(2)*q2(2) - q1(3)*q2(3) - q1(4)*q2(4)
        q3(2) = q1(1)*q2(2) + q1(2)*q2(1) + q1(3)*q2(4) - q1(4)*q2(3)
        q3(3) = q1(1)*q2(3) - q1(2)*q2(4) + q1(3)*q2(1) + q1(4)*q2(2)
        q3(4) = q1(1)*q2(4) + q1(2)*q2(3) - q1(3)*q2(2) + q1(4)*q2(1)

    end function quat_mult


    function quat_base_to_dependent(vec1, q1) result(vec2)

        implicit none

        real, intent(in) :: q1(4), vec1(3)
        real :: vec2(3)

        real :: t(4)

        t(1) = -vec1(1)*q1(2) - vec1(2)*q1(3) - vec1(3)*q1(4)
        t(2) = vec1(1)*q1(1) + vec1(2)*q1(4) - vec1(3)*q1(3)
        t(3) = -vec1(1)*q1(4) + vec1(2)*q1(1) + vec1(3)*q1(2)
        t(4) = vec1(1)*q1(3) - vec1(2)*q1(2) + vec1(3)*q1(1)

        vec2(1) = q1(1)*t(2) - q1(2)*t(1) - q1(3)*t(4) + q1(4)*t(3)
        vec2(2) = q1(1)*t(3) + q1(2)*t(4) - q1(3)*t(1) - q1(4)*t(2)
        vec2(3) = q1(1)*t(4) - q1(2)*t(3) + q1(3)*t(2) - q1(4)*t(1)

    end function quat_base_to_dependent


    function quat_dependent_to_base(vec1, q1) result(vec2)

        implicit none
        
        real, intent(in) :: vec1(3), q1(4)
        real :: vec2(3)

        real :: t(4)

        t(1) = vec1(1)*q1(2) + vec1(2)*q1(3) + vec1(3)*q1(4)
        t(2) = vec1(1)*q1(1) - vec1(2)*q1(4) + vec1(3)*q1(3)
        t(3) = -vec1(1)*q1(4) + vec1(2)*q1(1) - vec1(3)*q1(2)
        t(4) = -vec1(1)*q1(3) + vec1(2)*q1(2) + vec1(3)*q1(1)

        vec2(1) = q1(1)*t(2) + q1(2)*t(1) + q1(3)*t(4) - q1(4)*t(3)
        vec2(2) = q1(1)*t(3) - q1(2)*t(4) + q1(3)*t(1) + q1(4)*t(2)
        vec2(3) = q1(1)*t(4) + q1(2)*t(3) - q1(3)*t(2) + q1(4)*t(1)
        
    end function quat_dependent_to_base


    subroutine quat_norm(q1)

        implicit none
        
        real, intent(inout) :: q1(4)

        q1 = q1/sqrt(q1(1)**2 + q1(2)**2 + q1(3)**2 + q1(4)**2)

    end subroutine quat_norm


    function euler_to_quat(e1) result(q1)
        implicit none
        
        real, intent(in) :: e1(3)
        real :: q1(4)

        !q1(1) = cos(e1(1)*0.5)*cos(e1(2)*0.5)*cos(e1(3)*0.5) + sin(e1(1)*0.5)*sin(e1(2)*0.5)*sin(e1(3)*0.5)
        !q1(2) = sin(e1(1)*0.5)*cos(e1(2)*0.5)*cos(e1(3)*0.5) - cos(e1(1)*0.5)*sin(e1(2)*0.5)*sin(e1(3)*0.5)
        !q1(3) = cos(e1(1)*0.5)*sin(e1(2)*0.5)*cos(e1(3)*0.5) + sin(e1(1)*0.5)*cos(e1(2)*0.5)*sin(e1(3)*0.5)
        !q1(4) = cos(e1(1)*0.5)*cos(e1(2)*0.5)*sin(e1(3)*0.5) - sin(e1(1)*0.5)*sin(e1(2)*0.5)*cos(e1(3)*0.5)

        !real :: s(3), c(3)

        !s(1) = sin(e1(1)*0.5)
        !s(2) = sin(e1(2)*0.5)
        !s(3) = sin(e1(3)*0.5)

        !c(1) = cos(e1(1)*0.5)
        !c(2) = cos(e1(2)*0.5)
        !c(3) = cos(e1(3)*0.5)

        !q1(1) = c(1)*c(2)*c(3) + s(1)*s(2)*s(3)    
        !q1(2) = s(1)*c(2)*c(3) - c(1)*s(2)*s(3)    
        !q1(3) = c(1)*s(2)*c(3) + s(1)*c(2)*s(3)    
        !q1(4) = c(1)*c(2)*s(3) - s(1)*s(2)*c(3)    

        real :: cphi, ctheta, cpsi, sphi, stheta, spsi

        sphi = sin(e1(1)*0.5)
        stheta = sin(e1(2)*0.5)
        spsi = sin(e1(3)*0.5)

        cphi = cos(e1(1)*0.5)
        ctheta = cos(e1(2)*0.5)
        cpsi = cos(e1(3)*0.5)

        q1(1) = cphi*ctheta*cpsi + sphi*stheta*spsi    
        q1(2) = sphi*ctheta*cpsi - cphi*stheta*spsi    
        q1(3) = cphi*stheta*cpsi + sphi*ctheta*spsi    
        q1(4) = cphi*ctheta*spsi - sphi*stheta*cpsi    
        
    end function euler_to_quat


    function quat_to_euler(q1) result(e1)
        implicit none
        
        real, intent(in) :: q1(4)
        real :: e1(3)

        if (abs(q1(1)*q1(3) - q1(2)*q1(4) - 0.5) < 1e-12) then
            
            e1(1) = 2*asin(q1(2)/cos(PI4))
            e1(2) = PI2
            e1(3) = 0.0

        else if (abs(q1(1)*q1(3) - q1(2)*q1(4) + 0.5) < 1e-12) then

            e1(1) = 2*asin(q1(2)/cos(PI4))
            e1(2) = -PI2
            e1(3) = 0.0

        else

            e1(1) = atan2(2*(q1(1)*q1(2) + q1(3)*q1(4)), (q1(1)**2 + q1(4)**2 - q1(2)**2 - q1(3)**2))
            e1(2) = asin(2*(q1(1)*q1(3) - q1(2)*q1(4)))
            e1(3) = atan2(2*(q1(1)*q1(4) + q1(2)*q1(3)), (q1(1)**2 + q1(2)**2 - q1(3)**2 - q1(4)**2))

        end if
    
        
    end function quat_to_euler

    !!! STD ATMOSPHERE !!!

    function gravity_SI(h) result(g)
        implicit none
        
        real, intent(in) :: h 
        real :: g

        g = g_ssl*(r_ez/(r_ez + h))**2
        
    end function gravity_SI

    function gravity_English(h_ft) result(g)
        implicit none
        
        real, intent(in) :: h_ft
        real :: g

        g = meters_to_feet(gravity_SI(feet_to_meters(h_ft)))
        
    end function gravity_English


    subroutine std_atm_SI(h, Z, T, P, rho, a)

        implicit none

        real, intent(in) :: h
        real, intent(inout) :: Z, T, P, rho, a
        real :: h_p

        integer :: i

        ! THIS WAS MOVED TO THE TOP OF THE MODULE
        !real :: Z_i(8), Z_i1(8), T_i(8), Tp_i(8), P_i(8)
        ! TABLE DECLARATION DO NOT EDIT
        !Z_i = [0., 11000., 20000., 32000., 47000., 52000., 61000., 79000.]
        !Z_i1 = [11000., 20000., 32000., 47000., 52000., 61000., 79000., 90000.]
        !T_i = [288.150, 216.650, 216.650, 228.650, 270.650, 270.650, 252.650, 180.650]
        !Tp_i = [-6.5, 0.0, 1.0, 2.8, 0.0, -2.0, -4.0, 0.0]
        !P_i = [1.01325e5, 2.26320318222212e4, 5.47487352827083e3, 8.68014769086723e2, &
                !1.10905588989225e2, 5.90005242789244e1, 1.82099249050177e1, 1.03770045489203]

        if (h < 0.0) then
            h_p = 0.0
        else
            h_p = h
        end if

        ! Calculate geopotential altitude
        Z = r_ez*h_p/(r_ez + h_p)

        ! Calculate pressure
        do i = 1, 8

            if (Z >= Z_i(i) .and. Z < Z_i1(i)) then

                T = T_i(i) + 0.001*Tp_i(i)*(Z - Z_i(i))

                if (Tp_i(i) == 0) then
                    P = P_i(i)*exp(-g_ssl*(Z - Z_i(i))/(R_gas * T_i(i)))
                else
                    P = P_i(i)*(T/T_i(i))**(-g_ssl/(R_gas * 0.001*Tp_i(i)))
                end if 

                rho = P/(R_gas*T)

                a = sqrt(gamma*R_gas*T)

                exit
            else
                continue
            end if

        end do

    end subroutine std_atm_SI


    subroutine std_atm_English(h, Z, T, P, rho, a)

        implicit none
        
        real, intent(in) :: h
        real, intent(inout) :: Z, T, P, rho, a

        real :: h_m
        
        h_m = feet_to_meters(h)

        call std_atm_SI(h_m, Z, T, P, rho, a)

        Z = meters_to_feet(Z)
        T = kelvin_to_rankine(T)
        P = P*1.0/47.880258
        rho = rho*kg_m3_to_slug_ft3
        a = meters_to_feet(a)

    end subroutine std_atm_English

    function sutherland_visc_SI(T) result(mu)

        implicit none
        
        real, intent(in) :: T
        real :: mu

        mu = 1.716e-05*((273.15+110.4)/(T + 110.4))*(T/273.15)**(3./2.)

    end function sutherland_visc_SI

    function sutherland_visc_English(T) result(mu)

        implicit none
        
        real, intent(in) :: T
        real :: mu, T_K

        T_K = rankine_to_kelvin(T)

        mu = sutherland_visc_SI(T_K)*1.0/47.880258

    end function sutherland_visc_English


    !!! UNIT CONVERSIONS !!!

    function feet_to_meters(l_ft) result(l_m)

        implicit none
        
        real, intent(in) :: l_ft
        real :: l_m

        l_m = 0.3048*l_ft

    end function feet_to_meters

    function meters_to_feet(l_m) result(l_ft)

        implicit none
        
        real, intent(in) :: l_m
        real :: l_ft

        l_ft = 3.280839895013123*l_m


    end function meters_to_feet

    function kelvin_to_rankine(T) result(T_R)

        implicit none
        
        real, intent(in) :: T
        real :: T_R

        T_R = 1.8*T

    end function kelvin_to_rankine

    function rankine_to_kelvin(T) result(T_K)
        implicit none
        
        real, intent(in) :: T
        real :: T_K
    
        T_K = T*1.0/1.8

    end function rankine_to_kelvin

    function y_English_to_SI(y_eng) result(y_SI)

        implicit none

        real, intent(in) :: y_eng(13)
        real :: y_SI(13)

        y_SI(1:3) = y_eng(1:3)*0.3048
        y_SI(4:6) = y_eng(4:6)
        y_SI(7:9) = y_eng(7:9)*0.3048
        y_SI(10:13) = y_eng(10:13)

    end function y_English_to_SI

    function y_SI_to_English(y_SI) result(y_eng)

        implicit none

        real, intent(in) :: y_SI(13)
        real :: y_eng(13)

        y_eng(1:3) = y_SI(1:3)*1.0/0.3048
        y_eng(4:6) = y_SI(4:6)
        y_eng(7:9) = y_SI(7:9)*1.0/0.3048
        y_eng(10:13) = y_SI(10:13)

    end function y_SI_to_English

end module goates_m


!program test_gravity

    !use helper_m

    !implicit none

    !integer :: i, j
    !real :: g

    !write(*,*) "SI units"
    !do i = 0, 100000, 5000
        !g = gravity_SI(real(i))

        !write(*,*) "Altitude:", i, " meters. Gravity", g, "m/s2"

    !end do

    !write(*,*)
    !write(*,*) "ENGLISH UNITS"

    !do j = 0, 200000, 10000
        !g = gravity_English(real(j))

        !write(*,*) "Altitude:", j, " meters. Gravity", g, "ft/s2"

    !end do
    
!end program test_gravity